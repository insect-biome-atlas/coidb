#!/usr/bin/env python

from argparse import ArgumentParser
import polars as pl
from shutil import rmtree
from coidb.scripts.bold2fasta import multi_record_bins, specific_bins
from itertools import batched
import os


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
        "-o",
        "--outdir",
        type=str,
        help="Write batches into this output directory",
        required=True,
    )
    parser.add_argument(
        "-b",
        "--batch_size",
        type=int,
        help="Number of BINs in each batch",
        default=50000,
    )
    args = parser.parse_args()
    if os.path.exists(args.outdir):
        rmtree(args.outdir)
    os.makedirs(args.outdir, exist_ok=True)
    df = pl.scan_csv(args.infile, separator="\t", low_memory=True).select(
        "bin_uri", "processid", "seq"
    )
    df = multi_record_bins(df)
    bins = df.select("bin_uri").collect().unique().to_series().to_list()
    for i, batch in enumerate(batched(bins, args.batch_size), 1):
        outfile = os.path.join(args.outdir, f"batch_{i}.fasta")
        _df = specific_bins(df, batch)
        _df.with_columns(
            s=">"
            + pl.col("processid")
            + "|"
            + "bin_uri:"
            + pl.col("bin_uri")
            + "\n"
            + pl.col("seq")
        ).select("s").sink_csv(outfile, include_header=False, quote_style="never")
