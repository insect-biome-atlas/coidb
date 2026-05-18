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
    parser.add_argument(
        "-i", "--input_taxfile", required=True, help="TSV file with taxonomic labels"
    )
    parser.add_argument(
        "-m",
        "--matched",
        required=True,
        help="TSV file with taxonomic labels and a 'name' column corresponding to the original species name",
    )
    parser.add_argument(
        "-o",
        "--output_taxfile",
        required=True,
        help="TSV output file with consolidated naming",
    )
    args = parser.parse_args()
    ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    sys.stderr.write(f"Reading matched results from {args.matched}\n")
    matched_df = pl.scan_csv(args.matched, separator="\t")
    # generate unique rows
    matched_unique = generate_unique(df=matched_df, ranks=ranks)
    # filter to only species-matched results
    matched_unique = matched_unique.filter(~pl.col("species").str.contains(r"_X+$"))
    # sys.stderr.write(f"{matched_unique.height} unique taxa in matched results\n")
    # sys.stderr.write(f"Reading BOLD data from {args.input_taxfile}\n")
    input_df = pl.scan_csv(args.input_taxfile, separator="\t")
    input_unique = generate_unique(df=input_df, ranks=ranks)
    # identify potential errors where species in Arthropoda are assigned to a
    # different phylum. Examples include 'Murphyana rayi' (BOLD:AAI8355) which
    # is matched to phylum Mollusca in Catalogue of Life.
    # the code below identifies cases where
    # 1) the original and matched phyla differ
    # 2) the original phylum is Arthropoda and the matched phylum is not unassigned
    not_allowed = (
        input_unique.join(matched_unique, left_on="species", right_on="name")
        .filter(
            (pl.col("phylum") != pl.col("phylum_right"))
            & (~pl.col("phylum_right").str.contains(r"_X+$"))
            & (pl.col("phylum") == "Arthropoda")
        )
        .select(["kingdom", "phylum", "class", "order", "family", "genus", "species"])
    )
    fixed = not_allowed.join(matched_unique, on="species").select(
        ["name", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
    )
    untouched = matched_unique.join(not_allowed, on="species", how="anti")
    matched_unique = pl.concat([fixed, untouched])
    sys.stderr.write(f"Consolidating names and writing to {args.output_taxfile}\n")
    (
        input_df.drop(["kingdom", "phylum", "class", "order", "family", "genus"])
        .join(matched_unique, left_on="species", right_on="name")
        .drop("species")
        .rename({"species_right": "species"})
        .select(["processid"] + ranks + ["bin_uri", "seq"])
    ).sink_csv(args.output_taxfile, separator="\t")
