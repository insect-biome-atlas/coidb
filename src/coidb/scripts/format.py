#!/usr/bin/env python

from argparse import ArgumentParser
import polars as pl
import gzip as gz
from coidb import series_to_fasta, read_records


def format_dada2(df, outfile_toGenus, outfile_toSpecies, outfile_assignSpecies):
    """
    Generate files compatible with DADA2
    """
    df.with_columns(
        s=">"
        + pl.col("kingdom")
        + ";"
        + pl.col("phylum")
        + ";"
        + pl.col("class")
        + ";"
        + pl.col("order")
        + ";"
        + pl.col("family")
        + ";"
        + pl.col("genus")
        + ";\n"
        + pl.col("seq")
    ).select("s").sink_csv(outfile_toGenus, include_header=False, quote_style="never")
    # Generate a fasta up to the species level
    df.with_columns(
        s=">"
        + pl.col("kingdom")
        + ";"
        + pl.col("phylum")
        + ";"
        + pl.col("class")
        + ";"
        + pl.col("order")
        + ";"
        + pl.col("family")
        + ";"
        + pl.col("genus")
        + ";"
        + pl.col("species")
        + ";\n"
        + pl.col("seq")
    ).select("s").sink_csv(outfile_toSpecies, include_header=False, quote_style="never")
    # Generate a fasta with only species
    df.with_columns(
        s=">" + pl.col("processid") + " " + pl.col("species") + "\n" + pl.col("seq")
    ).select("s").sink_csv(
        outfile_assignSpecies, include_header=False, quote_style="never"
    )


def format_sintax(df, outfile):
    """
    Output a fasta file compatible with SINTAX
    Here the BOLD BIN is set as 'strain' ('t:' prefix)
    """
    df.with_columns(
        s=">"
        + pl.col("processid")
        + ";tax=k:"
        + pl.col("kingdom")
        + ",p:"
        + pl.col("phylum")
        + ",c:"
        + pl.col("class")
        + ",o:"
        + pl.col("order")
        + ",f:"
        + pl.col("family")
        + ",g:"
        + pl.col("genus")
        + ",s:"
        + pl.col("species")
        + ",t:"
        + pl.col("bin_uri")
        + "\n"
        + pl.col("seq")
    ).select("s").sink_csv(outfile, include_header=False, quote_style="never")


def format_qiime2(df, outfile):
    """
    Output a tab separated taxonomy file with columns 'Taxon' and 'Feature ID'
    which can be used with QIIME2
    """
    df.with_columns(
        [
            ("k__" + pl.col("kingdom")).alias("kingdom"),
            ("p__" + pl.col("phylum")).alias("phylum"),
            ("c__" + pl.col("class")).alias("class"),
            ("o__" + pl.col("order")).alias("order"),
            ("f__" + pl.col("family")).alias("family"),
            ("g__" + pl.col("genus")).alias("genus"),
            ("s__" + pl.col("species")).alias("species"),
            ("t__" + pl.col("bin_uri")).alias("bin_uri"),
        ]
        # concatenate rank columns
    ).with_columns(
        pl.concat_str(
            [
                "kingdom",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
                "bin_uri",
            ],
            separator="; ",
        ).alias("Taxon")
    ).select(
        ["processid", "Taxon"]
    ).rename(
        {"processid": "Feature ID"}
    ).sink_csv(
        outfile, separator="\t"
    )


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "--tsv",
        type=str,
        required=True,
        help="TSV file with taxonomic info and sequences",
    )
    parser.add_argument(
        "--fasta", type=str, required=False, help="Fasta file of sequences"
    )
    parser.add_argument(
        "--consensus",
        type=str,
        required=False,
        help="TSV file with consensus taxonomy for BOLD BINs",
    )
    subparsers = parser.add_subparsers(dest="format")
    sintax_parser = subparsers.add_parser("sintax", help="SINTAX format options")
    sintax_parser.add_argument(
        "--outfile", type=str, required=True, help="SINTAX formatted fasta"
    )
    dada2_parser = subparsers.add_parser("dada2", help="DADA2 format options")
    dada2_parser.add_argument(
        "--outfile_toGenus",
        type=str,
        required=True,
        help="Output file for assignments up to genus",
    )
    dada2_parser.add_argument(
        "--outfile_toSpecies",
        type=str,
        required=True,
        help="Output file for assignments up to species",
    )
    dada2_parser.add_argument(
        "--outfile_assignSpecies",
        type=str,
        required=True,
        help="Output file for assignSpecies",
    )
    qiime2_parser = subparsers.add_parser("qiime2", help="QIIME2 format options")
    qiime2_parser.add_argument(
        "--qiime2_outfile",
        type=str,
        required=True,
        help="TSV file with taxonomic info for use with QIIME2",
    )
    args = parser.parse_args()
    if args.fasta:
        records = read_records(args.fasta)
    if args.consensus:
        # Create a dataframe of only the records from the clustered fasta file
        # and only keep the processid, bin_uri and seq columns
        df = (
            pl.scan_csv(args.tsv, separator="\t")
            .filter(pl.col("processid").is_in(records))
            .select(["processid", "bin_uri", "seq"])
        )
        # Read consensus taxonomy for BOLD BINs
        consensus = pl.scan_csv(args.consensus, separator="\t")
        # Merge dataframe with consensus, adding BOLD BIN taxonomy to records
        df = consensus.join(df, on="bin_uri")
    else:
        df = pl.scan_csv(args.tsv, separator="\t")
    # Write the requested format
    if args.format == "sintax":
        format_sintax(df, args.outfile)
    elif args.format == "dada2":
        format_dada2(
            df,
            outfile_toGenus=args.outfile_toGenus,
            outfile_toSpecies=args.outfile_toSpecies,
            outfile_assignSpecies=args.outfile_assignSpecies,
        )
    elif args.format == "qiime2":
        format_qiime2(df, args.qiime2_outfile)
