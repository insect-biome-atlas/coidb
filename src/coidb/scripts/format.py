#!/usr/bin/env python

from argparse import ArgumentParser
import polars as pl
import gzip as gz
import sys


def read_records(f):
    """
    Read records from fasta file
    """
    d = {"processid": [], "bin_uri": []}
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
                record, bin_uri = line.lstrip(">").split(" ")
                bin_uri = bin_uri.lstrip("bin_uri:")
                d["processid"].append(record)
                d["bin_uri"].append(bin_uri)
    return pl.DataFrame(d).lazy()


def generate_id_file(ids, outfile):
    with open(outfile, "w") as fhout:
        for i in ids:
            fhout.write(f"{i}\n")


def generate_kv_file(df, outfile, format="sintax"):
    if format == "sintax":
        df = df.with_columns(
            value=";tax=k:"
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
        )
    elif format == "dada2.toGenus":
        df = df.with_columns(
            value=pl.col("kingdom")
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
        )
    elif format == "dada2.toSpecies":
        df = df.with_columns(
            value=pl.col("kingdom")
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
            + ";"
        )
    elif format == "dada2.addSpecies":
        df = df.with_columns(value=pl.col("species"))
    df = df.with_columns(
        key="bin_uri:" + pl.col("bin_uri"),
    )
    df.select("key", "value").sink_csv(outfile, separator="\t", include_header=False)


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
        "--fasta",
        type=str,
        required=True,
        help="Fasta file of sequences. Required when --format=qiime2",
    )
    parser.add_argument(
        "--consensus",
        type=str,
        required=True,
        help="TSV file with consensus taxonomy for BOLD BINs",
    )
    parser.add_argument(
        "--format",
        type=str,
        default="sintax",
        help="Output format. Choose from 'sintax', 'dada2 toGenus', 'dada2 toSpecies', 'dada2 addSpecies' or 'qiime2'",
        choices=[
            "sintax",
            "dada2.toGenus",
            "dada2.toSpecies",
            "dada2.addSpecies",
            "qiime2",
        ],
    )
    parser.add_argument("-o", "--outfile", type=str, required=True, help="Output file")
    parser.add_argument("--idfile", type=str, help="Fasta identifiers output file")
    args = parser.parse_args()
    # Ensure idfile argument passed if format != qiime2
    if args.format != "qiime2" and args.idfile is None:
        sys.exit("Argument --idfile required when format!=qiime2\n")
    # Read consensus taxonomy for BOLD BINs
    consensus = pl.scan_csv(args.consensus, separator="\t")
    sys.stderr.write(f"Reading fasta headers from {args.fasta}\n")
    records = read_records(args.fasta)
    df = consensus.join(records, on="bin_uri")
    sys.stderr.write(
        f"{df.collect().height}/{records.collect().height} records found in consensus\n"
    )
    if args.format == "qiime2":
        sys.stderr.write("Generating QIIME2 format file\n")
        sys.stderr.write(f"Writing QIIME2 format file to {args.outfile}\n")
        format_qiime2(df, args.outfile)
    else:
        sys.stderr.write(f"Writing {args.format} key-value file to {args.outfile}\n")
        generate_kv_file(df=consensus, outfile=args.outfile, format=args.format)
        generate_id_file(
            ids=df.select("processid").unique().collect().to_series().to_list(),
            outfile=args.idfile,
        )
