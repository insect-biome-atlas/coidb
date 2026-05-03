#!/usr/bin/env python

from pygbif import species
import polars as pl
from multiprocessing import Pool
from argparse import ArgumentParser
from tqdm import tqdm
import sys


def gbif_match(value):
    """
    Matches species names/bin URIs to Catalog of Life. Only returns a taxonomy
    if the matching is exact.
    """
    ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    res = species.name_backbone(
        scientificName=value,
        checklistKey="7ddf754f-d193-4cc9-b351-99906754a03b",
        strict=True,
        verbose=True,
    )
    d = {}
    taxres = {"name": value}
    for rank in ranks:
        taxres[rank] = None
    if (
        "alternatives" in res["diagnostics"].keys()
        and len(res["diagnostics"]["alternatives"]) > 0
    ):
        return taxres
    if res["diagnostics"]["matchType"] != "EXACT":
        return taxres
    for item in res["classification"]:
        rank = item["rank"].lower()
        name = item["name"]
        d[rank] = name
    for rank in ranks:
        try:
            taxres[rank] = d[rank]
        except KeyError:
            continue
    return taxres


def get_true_name(name):
    res = species.name_backbone(
        scientificName=name,
        strict=True,
        checklistKey="7ddf754f-d193-4cc9-b351-99906754a03b",
        verbose=True,
    )
    if res["diagnostics"]["matchType"] != "EXACT":
        return name
    classification_df = pl.DataFrame(res["classification"])
    rank = res["usage"]["rank"]
    true_name = classification_df.filter(pl.col("rank") == rank).item(0, 1)
    return true_name


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


def collapse_nulls(df):
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
    ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    # Create a unique dataframe with the least amount of null values
    uniq = (
        df.drop_nulls("genus")
        .select(ranks[0 : ranks.index("genus") + 1])
        .unique()
        .with_columns(nulls=pl.sum_horizontal(pl.col("*").is_null()))
        .sort("nulls")
        .head(1)
        .drop("nulls")
    )
    collapse = True
    for rank in ranks[0 : ranks.index("genus")]:
        if df.drop_nulls(rank).n_unique(rank) > 1:
            collapse = False
    if collapse and df.unique("kingdom").height == 1:
        return (
            df.drop(ranks[0 : ranks.index("genus")])
            .join(uniq, on="genus")
            .select(["name"] + ranks)
        )
    else:
        return df


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
        "-c", "--col", type=str, help="Column to match by", default="bin_uri"
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
    args = parser.parse_args()
    if not args.refine_only:
        sys.stderr.write(f"Reading unique values for {args.col} from {args.infile}\n")
        unique_ids = get_unique(args.infile, args.col)
        sys.stderr.write(f"{len(unique_ids)} unique values loaded\n")
        with Pool(args.cpus) as p:
            matches = pl.DataFrame(
                list(
                    tqdm(
                        p.imap_unordered(gbif_match, unique_ids),
                        total=len(unique_ids),
                        unit=f" {args.col}",
                        ncols=120,
                        leave=False,
                        desc=f"Matching {args.col} to Catalog of Life",
                    )
                )
            )
        matched_df = matches.select(["name"] + args.ranks)
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
    partitioned = matched_df.partition_by("genus")
    sys.stderr.write("Refining lineages per genera\n")
    with Pool(args.cpus) as p:
        refined_df = pl.concat(
            tqdm(
                p.imap_unordered(collapse_nulls, partitioned, chunksize=1),
                total=len(partitioned),
                unit=f" genera",
                ncols=120,
                leave=False,
                desc="Refining",
            ),
            how="diagonal_relaxed",
        )
    refined_null_counts = refined_df.drop("name").null_count()
    sys.stderr.write("Missing values per rank after refinement:\n")
    sys.stderr.write(f"{str(refined_null_counts)}\n")
    refined_df.write_csv(args.outfile, separator="\t")
