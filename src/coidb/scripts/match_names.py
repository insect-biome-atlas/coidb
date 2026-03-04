#!/usr/bin/env python

from pygbif import species
import polars as pl
from multiprocessing import Pool
from argparse import ArgumentParser
from tqdm import tqdm
import sys


def col_match(
    value, ranks=["kingdom", "phylum", "class", "order", "family", "genus", "species"]
):
    res = species.name_backbone(
        scientificName=value, checklistKey="7ddf754f-d193-4cc9-b351-99906754a03b"
    )
    d = {}
    taxres = {"name": value}
    for rank in ranks:
        taxres[rank] = None
    if not "classification" in res.keys():
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


def refine(v, d=None):
    res = species.name_backbone(
        scientificName=v, checklistKey="7ddf754f-d193-4cc9-b351-99906754a03b"
    )
    try:
        return d[v], d
    except KeyError:
        pass
    refined = v
    if not "classification" in res.keys():
        return refined, d
    classification_df = pl.DataFrame(res["classification"])
    rank = res["usage"]["rank"]
    refined = classification_df.filter(pl.col("rank") == rank).item(0, 1)
    d[v] = refined
    return refined, d


def refine_df(df):
    columns = df.columns
    refined_data = {}
    for col in columns:
        refined_data[col] = []
    d = {}
    for row in df.sort(columns[1:]).iter_rows():
        refined_data["name"].append(row[0])
        for i, value in enumerate(row[1:], start=1):
            col = columns[i]
            refined, d = refine(value, d)
            refined_data[col].append(refined)
    return pl.DataFrame(refined_data)


def get_unique(f, col="bin_uri"):
    info = pl.scan_csv(f, separator="\t")
    v = (
        (info.filter(~pl.col(col).str.contains(r"_X+$")).select(col).unique())
        .collect()
        .to_series()
        .to_list()
    )
    return v


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
        "--refine_partition",
        type=str,
        help="When refining names, partition data by this rank (default: family). Selecting a lower rank and giving more cpus can speed up runs.",
        default="family",
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
                        p.imap_unordered(col_match, unique_ids),
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
            sys.stderr.write(f"Writing unrefined table to {args.unrefined_out}\n")
            matched_df.write_csv(args.unrefined_out, separator="\t")
    else:
        sys.stderr.write(f"Only performing refinement of taxa in {args.infile}\n")
        matched_df = pl.read_csv(args.infile, separator="\t")
        matched_df.columns = ["name"] + matched_df.columns[1:]
        matched_df = matched_df.select(["name"] + args.ranks)
    partitioned = matched_df.partition_by(args.refine_partition)
    sys.stderr.write("Refining taxa names\n")
    with Pool(args.cpus) as p:
        refined_df = pl.concat(
            tqdm(
                p.imap_unordered(refine_df, partitioned),
                total=len(partitioned),
                unit=f" {args.refine_partition} partitions",
                ncols=120,
                leave=False,
                desc="Refining",
            )
        )
    refined_df.write_csv(args.outfile, separator="\t")
