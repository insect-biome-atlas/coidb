#!/usr/bin/env python

from argparse import ArgumentParser
import polars as pl


def main(infile, outfile):
    tsv = pl.scan_csv(infile, has_header=True, separator="\t", null_values=["None"])
    (
        tsv.with_columns(
            # if kingdom is null, set to "unassigned"
            pl.when(pl.col("kingdom").is_null())
            .then(pl.col("kingdom").fill_null("unassigned"))
            .otherwise(pl.col("kingdom"))
            .alias("kingdom")
        )
        .with_columns(
            # if phylum is null, set to kingdom + "_X"
            pl.when(pl.col("phylum").is_null())
            .then(pl.col("kingdom") + "_X")
            .otherwise(pl.col("phylum"))
            .alias("phylum")
        )
        .with_columns(
            # if class is null, and phylum does not end in "_X", set to ph_fixed + "_X"
            pl.when(
                (pl.col("class").is_null()) & (~pl.col("phylum").str.contains(r"_X+$"))
            )
            .then(pl.col("phylum") + "_X")
            .otherwise(pl.col("class"))
            .alias("class")
        )
        .with_columns(
            # if class is null, and phylum ends in "_X", set to phylum + "X"
            pl.when(
                (pl.col("class").is_null()) & (pl.col("phylum").str.contains(r"_X+$"))
            )
            .then(pl.col("phylum") + "X")
            .otherwise(pl.col("class"))
            .alias("class")
        )
        .with_columns(
            # if order is null and class does not end in "_X", add "_X"
            pl.when(
                (pl.col("order").is_null()) & (~pl.col("class").str.contains(r"_X+$"))
            )
            .then(pl.col("class") + "_X")
            .otherwise(pl.col("order"))
            .alias("order")
        )
        .with_columns(
            # if order is null and class ends with "_X", add "X"
            pl.when(
                (pl.col("order").is_null()) & (pl.col("class").str.contains(r"_X+$"))
            )
            .then(pl.col("class") + "X")
            .otherwise(pl.col("order"))
            .alias("order")
        )
        .with_columns(
            # if family is null and order does not end in "_X", add "_X"
            pl.when(
                (pl.col("family").is_null()) & (~pl.col("order").str.contains(r"_X+$"))
            )
            .then(pl.col("order") + "_X")
            .otherwise(pl.col("family"))
            .alias("family")
        )
        .with_columns(
            # if family is null and order ends in "_X", add "X"
            pl.when(
                (pl.col("family").is_null()) & (pl.col("order").str.contains(r"_X+$"))
            )
            .then(pl.col("order") + "X")
            .otherwise(pl.col("family"))
            .alias("family")
        )
        .with_columns(
            # if genus is null and family does not end in "_X", add "_X"
            pl.when(
                (pl.col("genus").is_null()) & (~pl.col("family").str.contains(r"_X+$"))
            )
            .then(pl.col("family") + "_X")
            .otherwise(pl.col("genus"))
            .alias("genus")
        )
        .with_columns(
            # if genus is null and family ends in "_X", add "X"
            pl.when(
                (pl.col("genus").is_null()) & (pl.col("family").str.contains(r"_X+$"))
            )
            .then(pl.col("family") + "X")
            .otherwise(pl.col("genus"))
            .alias("genus")
        )
        .with_columns(
            # if species is null and ge_fixed does not end in "_X", add "_X"
            pl.when(
                (pl.col("species").is_null()) & (~pl.col("genus").str.contains(r"_X+$"))
            )
            .then(pl.col("genus") + "_X")
            .otherwise(pl.col("species"))
            .alias("species")
        )
        .with_columns(
            # if species is null and ge_fixed ends in "_X", add "X"
            pl.when(
                (pl.col("species").is_null()) & (pl.col("genus").str.contains(r"_X+$"))
            )
            .then(pl.col("genus") + "X")
            .otherwise(pl.col("species"))
            .alias("species")
        )
    ).sink_csv(outfile, separator="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-i", "--infile", type=str, help="Input TSV file", required=True
    )
    parser.add_argument(
        "-o", "--outfile", type=str, help="Output TSV file", required=True
    )
    args = parser.parse_args()
    main(infile=args.infile, outfile=args.outfile)
