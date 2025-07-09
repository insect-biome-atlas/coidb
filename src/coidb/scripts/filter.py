#!/usr/bin/env python

import polars as pl
from argparse import ArgumentParser
import string


def filter_tsv(infile, outfile, min_len=0):
    letters = set(list(string.ascii_uppercase) + ["-"])
    DNA = set(["A", "C", "T", "G"])
    non_DNA = list(letters.difference(DNA))
    tsv = pl.scan_csv(
        infile,
        has_header=True,
        separator="\t",
        schema_overrides={"nuc_basecount": int},
        ignore_errors=True,
        low_memory=True,
    )
    (
        tsv.
        # select columns
        select(
            [
                "processid",
                "bin_uri",
                "kingdom",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
                "nuc",
                "nuc_basecount",
                "marker_code",
            ]
        )
        .
        # filter by gene type
        filter(
            (pl.col("marker_code") == "COI-5P")
            & (pl.col("bin_uri").str.contains("BOLD:[A-Z0-9]+"))
            & (pl.col("nuc_basecount") >= min_len)
        )
        .
        # translate seq to uppercase
        with_columns(seq=pl.col("nuc").str.to_uppercase())
        .
        # remove left gap characters
        with_columns(seq=pl.col("seq").str.strip_chars_start("-"))
        .
        # remove right gap characters
        with_columns(seq=pl.col("seq").str.strip_chars_end("-"))
        .
        # remove remaining seqs with non DNA characters
        filter(~pl.col("seq").str.contains_any(non_DNA))
        .
        # select and order columns
        select(
            [
                "processid",
                "kingdom",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
                "bin_uri",
                "seq",
            ]
        )
    ).sink_csv(outfile, separator="\t")


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "-i", "--infile", type=str, help="Input TSV file", required=True
    )
    parser.add_argument(
        "-o", "--outfile", type=str, help="Output TSV file", required=True
    )
    parser.add_argument(
        "-l",
        "--min_len",
        type=int,
        help="Minimum length of sequences to include",
        default=0,
    )
    args = parser.parse_args()
    filter_tsv(args.infile, args.outfile, args.min_len)
