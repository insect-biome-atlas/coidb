#!/usr/bin/env python

from argparse import ArgumentParser
import polars as pl


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "-i", "--infile", type=str, help="Input TSV file", required=True
    )
    parser.add_argument(
        "-o", "--outfile", type=str, help="Output TSV file", required=True
    )
    args = parser.parse_args()

    tsv = pl.scan_csv(
        args.infile, has_header=True, separator="\t", null_values=["None"]
    )
    columns = tsv.collect_schema().names()
    ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    for rank in ranks:
        if rank in columns:
            if rank == "kingdom":
                # if kingdom is null, set to "unassigned"
                tsv = tsv.with_columns(
                    pl.when(pl.col(rank).is_null())
                    .then(pl.col(rank).fill_null("unassigned"))
                    .otherwise(pl.col(rank))
                    .alias(rank)
                )
            # if phylum is null, set to kingdom + "_X"
            elif rank == "phylum":
                tsv = tsv.with_columns(
                    pl.when(pl.col(rank).is_null())
                    .then(pl.col("kingdom") + "_X")
                    .otherwise(pl.col(rank))
                    .alias(rank)
                )
            else:
                parent = ranks[ranks.index(rank) - 1]
                # if rank is class, order, family, genus or species
                # and parent rank does not end in "_X", set to parent + "_X"
                # otherwise set to parent + "X"
                tsv = tsv.with_columns(
                    pl.when(
                        (pl.col(rank).is_null())
                        & (~pl.col(parent).str.contains(r"_X+$"))
                    )
                    .then(pl.col(parent) + "_X")
                    .otherwise(pl.col(rank))
                    .alias(rank)
                )
                tsv = tsv.with_columns(
                    pl.when(
                        (pl.col(rank).is_null())
                        & (pl.col(parent).str.contains(r"_X+$"))
                    )
                    .then(pl.col(parent) + "X")
                    .otherwise(pl.col(rank))
                    .alias(rank)
                )
    tsv.sink_csv(args.outfile, separator="\t")
