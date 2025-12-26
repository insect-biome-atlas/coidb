#!/usr/bin/env python

from argparse import ArgumentParser
import polars as pl


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--infile",
        type=str,
        help="Input TSV file with 'bin_uri' column",
        required=True,
    )
    parser.add_argument(
        "-o", "--outfile", type=str, help="Output fasta file", required=True
    )
    parser.add_argument(
        "--single", action="store_true", help="Only output single-record BOLD bins"
    )
    parser.add_argument(
        "--multi", action="store_true", help="Only output multi-record BOLD bins"
    )
    parser.add_argument(
        "--id_col",
        type=str,
        default="processid",
        help="Name of column to use for fasta sequence id (default: 'processid')",
    )
    parser.add_argument(
        "--seq_col",
        type=str,
        default="seq",
        help="Name of column containing sequences (default: 'seq')",
    )
    parser.add_argument(
        "--low_memory",
        action="store_true",
        help="Read from file with polars 'low_memory' to reduce memory pressure at the expense of performance.",
    )
    args = parser.parse_args()
    df = pl.scan_csv(args.infile, separator="\t", low_memory=args.low_memory).select(
        "bin_uri", args.id_col, args.seq_col
    )
    if args.single:
        single_record_bins = (
            (df.group_by("bin_uri").len().filter(pl.col("len") == 1).select("bin_uri"))
            .collect(engine="streaming")
            .to_series()
            .to_list()
        )
        df = df.filter(pl.col("bin_uri").is_in(single_record_bins))
    elif args.multi:
        multi_record_bins = (
            df.filter(pl.col("bin_uri").str.starts_with("BOLD:"))
            .group_by("bin_uri")
            .len()
            .filter(pl.col("len") > 1)
            .select("bin_uri")
            .collect(engine="streaming")
            .to_series()
            .to_list()
        )
        df = df.filter(pl.col("bin_uri").is_in(multi_record_bins))
    df.with_columns(s=">" + pl.col(args.id_col) + "\n" + pl.col(args.seq_col)).select(
        "s"
    ).sink_csv(args.outfile, include_header=False, quote_style="never")
