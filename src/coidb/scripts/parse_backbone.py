#!/usr/bin/env python
from argparse import ArgumentParser
import polars as pl
from tempfile import NamedTemporaryFile
import sys
import os
from coidb import get_header
from coidb import extract_columns


def parse_backbone(infile, proper_header):
    """
    Parses the GBIF backbone into a dataframe of BOLD bins (first column) and
    taxonomic labels for ranks from kingdom->species.
    """
    tsv = pl.scan_csv(infile, separator="\t", low_memory=True)
    backbone = (tsv.select(proper_header)).collect()
    # Create dataframe with BOLD bins
    backbone_bins = backbone.filter(
        pl.col("scientificName").str.starts_with("BOLD:")
    ).rename({"scientificName": "bin_uri"})
    # Create taxonID->species dataframe
    species_df = (
        backbone.filter(
            (
                pl.col("taxonID").is_in(
                    backbone_bins.select("parentNameUsageID")
                    .unique()
                    .to_series()
                    .to_list()
                )
            )
            & (pl.col("taxonRank") == "species")
        )
        .select(["taxonID", "canonicalName"])
        .rename({"canonicalName": "species"})
    )
    # Add species names for BINs with species assignments
    bin_df = backbone_bins.join(
        species_df, left_on="parentNameUsageID", right_on="taxonID", how="left"
    )
    return bin_df.select(
        ["bin_uri", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
    )


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "infile",
        type=str,
        help="Taxon.tsv file with columns 'taxonID', 'parentNameUsageID', 'scientificName', 'canonicalName', 'taxonRank', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus'",
    )
    parser.add_argument(
        "outfile",
        type=str,
        help="Output file with BOLD_BIN in first column followed by taxonomic ranks (kingdom..species)",
    )
    args = parser.parse_args()
    proper_header = [
        "taxonID",
        "parentNameUsageID",
        "scientificName",
        "canonicalName",
        "taxonRank",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
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
    backbone = parse_backbone(infile, proper_header=proper_header)
    sys.stderr.write(f"Writing parsed backbone to {args.outfile}\n")
    backbone.write_csv(args.outfile, separator="\t")
    if infile != args.infile:
        sys.stderr.write(f"Removing {infile}\n")
        os.remove(infile)
