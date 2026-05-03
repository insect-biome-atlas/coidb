#!/usr/bin/env python

import polars as pl
from argparse import ArgumentParser
import sys


def generate_complete(df, ranks=None):
    if ranks is None:
        ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    for rank in ranks:
        df = df.filter(~pl.col(rank).str.contains(r"_X+$"))
    return df


def generate_unique(df, ranks=None):
    if ranks is None:
        ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    return df.unique(ranks)


def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_taxfile", help="TSV file with taxonomic labels")
    parser.add_argument(
        "-m",
        "--matched",
        help="TSV file with taxonomic labels and a 'name' column corresponding to the original species name",
    )
    parser.add_argument(
        "--filter_strategy",
        choices=["complete", "species-complete"],
        default="complete",
        help="Strategy to filter the matched data. With 'complete' (default) only matched species without missing taxlabels are used. With 'species-complete' missing taxlabels are allowed for ranks higher than species",
    )
    parser.add_argument(
        "-o", "--output_taxfile", help="TSV output file with consolidated naming"
    )
    args = parser.parse_args()
    ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    sys.stderr.write(f"Reading matched results from {args.matched}\n")
    matched_df = pl.scan_csv(args.matched, separator="\t")
    matched_unique = generate_unique(df=matched_df, ranks=ranks).collect()
    sys.stderr.write(f"{matched_unique.height} unique taxa in matched results\n")
    if args.filter_strategy == "complete":
        _ranks = ranks
    else:
        _ranks = ["species"]
    matched_complete = generate_complete(df=matched_unique, ranks=_ranks)
    sys.stderr.write(
        f"{matched_complete.height} taxa after removing taxa with missing data for {" ".join(_ranks)}\n"
    )
    sys.stderr.write(f"Reading BOLD data from {args.input_taxfile}\n")
    input_df = pl.scan_csv(args.input_taxfile, separator="\t")
    sys.stderr.write(f"Consolidating names and writing to {args.output_taxfile}\n")
    (
        input_df.drop(["kingdom", "phylum", "class", "order", "family", "genus"])
        .join(matched_complete.lazy(), left_on="species", right_on="name")
        .drop("species")
        .rename({"species_right": "species"})
        .select(["processid"] + ranks + ["bin_uri", "seq"])
    ).sink_csv(args.output_taxfile, separator="\t")
