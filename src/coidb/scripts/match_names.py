#!/usr/bin/env python

from multiprocessing import get_context
from pygbif import species
import polars as pl
from multiprocessing import Pool
from argparse import ArgumentParser
from tqdm import tqdm
import sys


def check_alternatives(alternatives, confidence):
    for a in alternatives:
        conf = a["diagnostics"]["confidence"]
        if conf >= confidence:
            return True
    return False


def gbif_match(
    value,
    rank="species",
    ranks=None,
    strict=True,
    checklist_key=None,
):
    """
    Matches species names/bin URIs to Catalog of Life. Only returns a taxonomy
    if the matching is exact.
    """
    if ranks is None:
        ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    res = species.name_backbone(
        scientificName=value,
        checklistKey=checklist_key,
        strict=strict,
        taxonRank=rank,
        verbose=True,
    )
    d = {}
    taxres = {"name": value}
    for rank in ranks:
        taxres[rank] = None
    if res["diagnostics"]["matchType"] != "EXACT":
        return pl.DataFrame(taxres)
    confidence = res["diagnostics"]["confidence"]
    if (
        "alternatives" in res["diagnostics"].keys()
        and len(res["diagnostics"]["alternatives"]) > 0
    ):
        alternatives = res["diagnostics"]["alternatives"]
        equal_best = check_alternatives(alternatives, confidence)
        if equal_best:
            return pl.DataFrame(taxres)
    for item in res["classification"]:
        rank = item["rank"].lower()
        name = item["name"]
        d[rank] = name
    for rank in ranks:
        try:
            taxres[rank] = d[rank]
        except KeyError:
            continue
    return pl.DataFrame(taxres)


def get_unique(f, col="bin_uri"):
    """
    Return unique values for the column.
    """
    info = pl.scan_csv(f, separator="\t")
    v = (
        (info.filter(~pl.col(col).str.contains(r"_X+$")).select(col).unique())
        .collect()
        .to_series()
        .to_list()
    )
    return v


def collapse_nulls(df, partition_rank="genus", ranks=None):
    """
    Replace null values for higher ranks if all non-null values are the same.
    Example:
    "Animalia" "Arthropoda" "Arachnida"	"Trombidiformes" "Unionicolidae" "Unionicola" "Unionicola crassipes"
    "Animalia" "Arthropoda" null	    "Trombidiformes" "Unionicolidae" "Unionicola" "Unionicola figuralis"
    "Animalia"	null        null         null	         null            "Unionicola" "Unionicola trapezidens"

    Result:
    "Animalia" "Arthropoda" "Arachnida"	"Trombidiformes" "Unionicolidae" "Unionicola" "Unionicola crassipes"
    "Animalia" "Arthropoda" "Arachnida" "Trombidiformes" "Unionicolidae" "Unionicola" "Unionicola figuralis"
    "Animalia" "Arthropoda" "Arachnida"	"Trombidiformes" "Unionicolidae" "Unionicola" "Unionicola trapezidens"
    """
    if ranks is None:
        ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    # Create a unique dataframe with the least amount of null values
    uniq = (
        df.drop_nulls(partition_rank)
        .select(ranks[0 : ranks.index(partition_rank) + 1])
        .unique()
        .with_columns(nulls=pl.sum_horizontal(pl.col("*").is_null()))
        .sort("nulls")
        .head(1)
        .drop("nulls")
    )
    # require that the kingdom has only one unique non-null label
    if df.drop_nulls("kingdom").unique("kingdom").height > 1:
        return df
    # iterate the ranks
    for rank in ranks[0 : ranks.index(partition_rank)]:
        # if there are more than one unique non-null taxlabel for the rank,
        # return the original df
        if df.drop_nulls(rank).n_unique(rank) > 1:
            return df
    # return the collapsed df
    return (
        df.drop(ranks[0 : ranks.index(partition_rank)])
        .join(uniq, on=partition_rank)
        .select(["name"] + ranks)
    )


def refine_worker(arg):
    df, partition_rank, ranks = arg
    return collapse_nulls(df=df, partition_rank=partition_rank, ranks=ranks)


def worker(arg):
    """
    Helper function allowing more than 1 argument to be passed.
    """
    value, rank, ranks, strict, checklist_key = arg
    return gbif_match(
        value=value, rank=rank, ranks=ranks, strict=strict, checklist_key=checklist_key
    )


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--infile",
        type=str,
        help="TSV infile with taxonomy for processids",
        required=True,
    )
    parser.add_argument(
        "-c", "--col", type=str, help="Column to match by", default="species"
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        help="TSV outfile with taxonomic ranks matched to backbone",
        required=True,
    )
    parser.add_argument(
        "--unrefined_out", type=str, help="Write un-refined dataframe to this outfile"
    )
    parser.add_argument(
        "-p", dest="cpus", type=int, default=1, help="Number of cpus to use"
    )
    parser.add_argument(
        "--refine_only",
        action="store_true",
        help="Skip matching names in input and only perform refinement",
    )
    parser.add_argument(
        "-r",
        "--ranks",
        nargs="+",
        help="Ranks to write taxonomic information for",
        default=["kingdom", "phylum", "class", "order", "family", "genus", "species"],
    )
    parser.add_argument("--strict", action="store_true", help="Use strict matching")
    parser.add_argument(
        "--checklist_key",
        type=str,
        help="Checklist key to use for matching",
    )
    args = parser.parse_args()
    if not args.refine_only:
        sys.stderr.write(f"Reading unique values for {args.col} from {args.infile}\n")
        unique_ids = get_unique(args.infile, args.col)
        sys.stderr.write(f"{len(unique_ids)} unique values loaded\n")
        with get_context("spawn").Pool(args.cpus) as p:
            df_list = list(
                tqdm(
                    p.imap_unordered(
                        worker,
                        (
                            (
                                value,
                                args.col,
                                args.ranks,
                                args.strict,
                                args.checklist_key,
                            )
                            for value in unique_ids
                        ),
                    ),
                    unit=f" {args.col}",
                    leave=False,
                    desc=f"Matching {args.col} names to GBIF",
                    total=len(unique_ids),
                    ncols=120,
                )
            )
        matched_df = pl.concat(df_list, how="vertical_relaxed").select(
            ["name"] + args.ranks
        )
        if args.unrefined_out:
            sys.stderr.write(f"Writing matched table to {args.unrefined_out}\n")
            matched_df.write_csv(args.unrefined_out, separator="\t")
    else:
        sys.stderr.write(f"Only performing refinement of taxa in {args.infile}\n")
        matched_df = pl.read_csv(args.infile, separator="\t")
        matched_df.columns = ["name"] + matched_df.columns[1:]
        matched_df = matched_df.select(["name"] + args.ranks)
    orig_null_counts = matched_df.drop("name").null_count()
    sys.stderr.write("Missing values per rank:\n")
    sys.stderr.write(f"{str(orig_null_counts)}\n")
    partition_rank = args.ranks[args.ranks.index(args.col) - 1]
    partitioned = matched_df.partition_by(partition_rank)
    sys.stderr.write(f"Refining lineages per {partition_rank}\n")
    with get_context("spawn").Pool(args.cpus) as p:
        refined_list = list(
            tqdm(
                p.imap(
                    refine_worker,
                    ((df, partition_rank, args.ranks) for df in partitioned),
                ),
                unit=f" {partition_rank}",
                leave=False,
                desc="Refining",
                total=len(partitioned),
                ncols=120,
            )
        )
    refined_df = pl.concat(refined_list, how="diagonal_relaxed")
    refined_null_counts = refined_df.drop("name").null_count()
    sys.stderr.write("Missing values per rank after refinement:\n")
    sys.stderr.write(f"{str(refined_null_counts)}\n")
    refined_df.write_csv(args.outfile, separator="\t")
