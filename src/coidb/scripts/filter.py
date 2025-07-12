#!/usr/bin/env python

import polars as pl
from argparse import ArgumentParser
import string
import sys
import os
from coidb import get_header, extract_columns


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
    # proper_header defines the columns required for filtering
    proper_header = [
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
    infile = args.infile
    header = get_header(args.infile)
    # If the header does not match proper_header, get indices of the required
    # columns and use extract_non_proper to write these to a temporary file
    if header != proper_header:
        missing_cols = ",".join(list(set(proper_header).difference(header)))
        if len(missing_cols) > 0:
            missing_cols = ",".join(list(set(proper_header).difference(header)))
            sys.exit(
                f"Not all required columns found in {args.infile}. Missing columns: {missing_cols}"
            )
        indices = [header.index(x) for x in proper_header]
        infile = extract_columns(args.infile, indices)
        sys.stderr.write(
            f"Wrote columns {','.join([str(x) for x in indices])} from {args.infile} to {infile}\n"
        )
    sys.stderr.write(f"Reading from {infile}\n")
    filter_tsv(infile, args.outfile, args.min_len)
    if infile != args.infile:
        sys.stderr.write(f"Removing {infile}\n")
        os.remove(infile)
