#!/usr/bin/env python

from argparse import ArgumentParser
import polars as pl
import gzip as gz
from coidb import series_to_fasta


def format_dada2(df, outfile_toGenus, outfile_toSpecies, outfile_assignSpecies):
    """
    Generate files compatible with DADA2
    """
    # Generate a fasta up to the genus level
    togenus = (
        ">"
        + df["kingdom"]
        + ";"
        + df["phylum"]
        + ";"
        + df["class"]
        + ";"
        + df["order"]
        + ";"
        + df["family"]
        + ";"
        + df["genus"]
        + ";"
        + "\n"
        + df["seq"]
    )
    if outfile_toGenus.endswith(".gz"):
        compress = True
    else:
        compress = False
    series_to_fasta(togenus, outfile_toGenus, compress=compress)
    # Generate a fasta up to the species level
    tospecies = (
        ">"
        + df["kingdom"]
        + ";"
        + df["phylum"]
        + ";"
        + df["class"]
        + ";"
        + df["order"]
        + ";"
        + df["family"]
        + ";"
        + df["genus"]
        + ";"
        + df["species"]
        + ";"
        + "\n"
        + df["seq"]
    )
    if outfile_toSpecies.endswith(".gz"):
        compress = True
    else:
        compress = False
    series_to_fasta(tospecies, outfile_toSpecies, compress=compress)
    # Generate a fasta with only species
    assignspecies = ">" + df["processid"] + " " + df["species"] + "\n" + df["seq"]
    if outfile_assignSpecies.endswith(".gz"):
        compress = True
    else:
        compress = False
    series_to_fasta(assignspecies, outfile_assignSpecies, compress=compress)


def format_sintax(df, outfile):
    """
    Output a fasta file compatible with SINTAX
    Here the BOLD BIN is set as 'strain' ('t:' prefix)
    """
    fasta_series = (
        ">"
        + df["processid"]
        + ";tax=k:"
        + df["kingdom"]
        + ",p:"
        + df["phylum"]
        + ",c:"
        + df["class"]
        + ",o:"
        + df["order"]
        + ",f:"
        + df["family"]
        + ",g:"
        + df["genus"]
        + ",s:"
        + df["species"]
        + ",t:"
        + df["bin_uri"]
        + "\n"
        + df["seq"]
    )
    if outfile.endswith(".gz"):
        compress = True
    else:
        compress = False
    series_to_fasta(fasta_series, outfile, compress=compress)


def read_records(f):
    """
    Read records from fasta file
    """
    records = []
    if f.endswith(".gz"):
        open_fn = gz.open
        mode = "rt"
    else:
        open_fn = open
        mode = "r"
    with open_fn(f, mode) as fhin:
        for line in fhin:
            line = line.rstrip()
            if line.startswith(">"):
                records.append(line.lstrip(">"))
    return records


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
            .collect()
        )
        # Read consensus taxonomy for BOLD BINs
        consensus = pl.read_csv(args.consensus, separator="\t")
        # Merge dataframe with consensus, adding BOLD BIN taxonomy to records
        df = consensus.join(df, on="bin_uri")
    else:
        df = pl.scan_csv(args.tsv, separator="\t").collect()
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
