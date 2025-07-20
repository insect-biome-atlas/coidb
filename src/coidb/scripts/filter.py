#!/usr/bin/env python

import polars as pl
from argparse import ArgumentParser
import string
import sys
import os
from coidb import get_header, extract_columns
import time


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
    # select columns
    select_cols = tsv.select(
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
    # filter rows by marker_code, length and BOLD BIN assignment
    # records assigned to Bacteria or Archaea are kept even if they do not
    # have a proper BIN assigned
    filter_rows = select_cols.filter(
        (pl.col("marker_code") == "COI-5P")
        & (pl.col("nuc_basecount") >= min_len)
        & (
            (pl.col("bin_uri").str.contains("BOLD:[A-Z0-9]+"))
            | (
                (pl.col("bin_uri") == "None")
                & (pl.col("kingdom").is_in(["Bacteria", "Archaea"]))
            )
        )
    )
    transform_seq = (
        # translate seq to uppercase
        filter_rows.with_columns(seq=pl.col("nuc").str.to_uppercase())
        # remove left gap characters
        .with_columns(seq=pl.col("seq").str.strip_chars_start("-"))
        # remove right gap characters
        .with_columns(seq=pl.col("seq").str.strip_chars_end("-"))
    )
    filter_seqs = (
        # remove remaining seqs with non DNA characters
        transform_seq.filter(~pl.col("seq").str.contains_any(non_DNA))
    )
    order_cols = (
        # select and order columns
        filter_seqs.select(
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
    )
    order_cols.with_columns(
        bin_uri=pl.when(pl.col("bin_uri") == "None")
        .then("processid")
        .otherwise("bin_uri")
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
    sys.stderr.write(f"Prefiltering {args.infile}\n")
    prefilter_start = time.time()
    tempfile = extract_columns(args.infile, indices)
    prefilter_time = round(time.time() - prefilter_start)
    sys.stderr.write(
        f"Wrote columns {','.join([str(x) for x in indices])} from {args.infile} to temporary file {tempfile} in {prefilter_time} s\n"
    )
    sys.stderr.write(f"Running filter on {tempfile} and writing to {args.outfile}\n")
    filter_tsv(tempfile, args.outfile, args.min_len)
    sys.stderr.write(f"Removing {tempfile}\n")
    os.remove(tempfile)
