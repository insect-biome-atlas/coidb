#!/usr/bin/env python

import polars as pl
from tqdm import tqdm
import gzip as gz
from argparse import ArgumentParser
import os


def count_bin_clusters(f):
    """
    Count sequences and BINs in the clustered fasta file

    :param f: Fasta file output from clustering
    """
    bins = {}
    seqids = []
    with gz.open(f, "rt") as fhin:
        for line in tqdm(fhin, desc=f"Reading {f}", unit=" lines"):
            if line.startswith(">"):
                line = line.rstrip()
                bin_uri = line.split(" ")[1].replace("bin_uri:", "")
                seqid = line.split(" ")[0].lstrip(">")
                seqids.append(seqid)
                try:
                    bins[bin_uri] += 1
                except KeyError:
                    bins[bin_uri] = 1
            continue
    bins_df = pl.LazyFrame(data={"bin_uri": bins.keys(), "n": bins.values()})
    mean_clusters_per_bins = (
        bins_df.filter(pl.col("bin_uri").str.starts_with("BOLD:"))
        .select("n")
        .mean()
        .collect()
        .item(0, 0)
    )
    median_clusters_per_bins = (
        bins_df.filter(pl.col("bin_uri").str.starts_with("BOLD:"))
        .select("n")
        .median()
        .collect()
        .item(0, 0)
    )
    max_clusters_per_bins = (
        bins_df.filter(pl.col("bin_uri").str.starts_with("BOLD:"))
        .select("n")
        .max()
        .collect()
        .item(0, 0)
    )
    return (
        seqids,
        bins_df,
        mean_clusters_per_bins,
        median_clusters_per_bins,
        max_clusters_per_bins,
    )


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "--fasta", required=True, type=str, help="COIDB clustered fasta.gz file"
    )
    parser.add_argument(
        "--consensus", required=True, type=str, help="Consensus taxonomy TSV file"
    )
    parser.add_argument(
        "--general_stats_out",
        required=True,
        type=str,
        help="Output file for general stats",
    )
    parser.add_argument(
        "--taxa_stats_out",
        required=True,
        type=str,
        help="Counts of bins and sequences in kingdoms/phyla",
    )
    args = parser.parse_args()
    fasta = args.fasta
    consensus = args.consensus
    cons_type = os.path.splitext(os.path.basename(consensus))[0].replace(".tsv", "")
    (
        seqids,
        bins_df,
        mean_clusters_per_bins,
        median_clusters_per_bins,
        max_clusters_per_bins,
    ) = count_bin_clusters(fasta)
    consensus_joined_df = pl.scan_csv(consensus, separator="\t").join(
        bins_df, on="bin_uri"
    )
    # calculate stats on sequences per bin
    bold_bin_df = consensus_joined_df.filter(pl.col("bin_uri").str.starts_with("BOLD:"))
    mean_seqs_per_bin = bold_bin_df.select("n").mean().collect().item(0, 0)
    median_seqs_per_bin = bold_bin_df.select("n").median().collect().item(0, 0)
    min_seqs_per_bin = bold_bin_df.select("n").min().collect().item(0, 0)
    max_seqs_per_bin = bold_bin_df.select("n").max().collect().item(0, 0)
    # calculate total sequences
    total_seqs = consensus_joined_df.select("n").sum().collect().item(0, 0)
    # calculate total bins
    total_bins = bold_bin_df.collect().height
    # calculate total non-BOLD-bins
    total_nonbins = (
        consensus_joined_df.filter(~pl.col("bin_uri").str.starts_with("BOLD:"))
        .collect()
        .height
    )
    seqs_per_kingdom = (
        consensus_joined_df.group_by("kingdom")
        .agg(pl.sum("n"))
        .rename({"n": "n_seqs"})
        .collect()
        .sort("kingdom")
    )
    bins_per_kingdom = (
        bold_bin_df.group_by("kingdom")
        .len()
        .rename({"len": "n_bins"})
        .sort("kingdom")
        .collect()
    )
    seqs_per_phylum = (
        consensus_joined_df.group_by("phylum")
        .agg(pl.sum("n"))
        .rename({"n": "n_seqs"})
        .collect()
        .sort("phylum")
    )
    bins_per_phyla = (
        bold_bin_df.group_by("phylum")
        .len()
        .rename({"len": "n_bins"})
        .collect()
        .sort("phylum")
    )
    # calculate total unique species
    total_species = consensus_joined_df.select("species").unique().collect().height
    # calculate total unique species assigned to BOLD BINs
    total_bin_species = bold_bin_df.select("species").unique().collect().height
    total_nonbin_species = bold_bin_df.select("species").unique().collect().height
    # filter to ambiguous species
    ambig_species = consensus_joined_df.filter(
        (pl.col("species").str.contains(r"_X+$"))
        & (~pl.col("species").str.starts_with("unresolved"))
    )
    # calculate sequences in ambiguous species
    ambig_species_seqs = ambig_species.select("n").sum().collect().item(0, 0)
    # filter to ambiguous species assigned to BOLD BINs
    ambig_bin_species = ambig_species.filter(pl.col("bin_uri").str.starts_with("BOLD:"))
    # calculate sequences in ambiguous species assigned to BOLD BINs
    ambig_bin_species_seqs = ambig_bin_species.select("n").sum().collect().item(0, 0)
    # filter to unresolved species
    unresolved_species = consensus_joined_df.filter(
        (pl.col("species").str.starts_with("unresolved"))
        & (~pl.col("species").str.contains(r"_X+$"))
    )
    # calculate sequences in unresolved species
    unresolved_species_seqs = unresolved_species.select("n").sum().collect().item(0, 0)
    # filter to unresolved species assigned to BOLD BINs
    unresolved_bin_species = unresolved_species.filter(
        pl.col("bin_uri").str.starts_with("BOLD:")
    )
    # calculate sequences in unresolved species assigned to BOLD BINs
    unresolved_bin_species_seqs = (
        unresolved_bin_species.select("n").sum().collect().item(0, 0)
    )
    unresolved_ambig_species = consensus_joined_df.filter(
        (pl.col("species").str.starts_with("unresolved"))
        & (pl.col("species").str.contains(r"_X+$"))
    )
    unresolved_ambig_species_seqs = (
        unresolved_ambig_species.select("n").sum().collect().item(0, 0)
    )
    unresolved_ambig_bin_species = unresolved_ambig_species.filter(
        pl.col("bin_uri").str.starts_with("BOLD:")
    )
    unresolved_ambig_bin_species_seqs = (
        unresolved_ambig_bin_species.select("n").sum().collect().item(0, 0)
    )
    general_stats = pl.DataFrame(
        data={
            "type": [cons_type],
            "total_seqs": [total_seqs],
            "total_bins": [total_bins],
            "mean_seqs_per_bin": [mean_seqs_per_bin],
            "median_seqs_per_bin": [median_seqs_per_bin],
            "min_seqs_per_bin": [min_seqs_per_bin],
            "max_seqs_per_bin": [max_seqs_per_bin],
            "total_non-bins": [total_nonbins],
            "total_species": [total_species],
            "total_bin_species": [total_bin_species],
            "total_nonbin_species": [total_nonbin_species],
            "ambiguous_species": [
                ambig_species.select("species").unique().collect().height
            ],
            "seqs_in_ambiguous_species": [ambig_species_seqs],
            "ambiguous_bin_species": [
                ambig_bin_species.select("species").unique().collect().height
            ],
            "seqs_in_ambiguous_bin_species": [ambig_bin_species_seqs],
            "unresolved_species": [
                unresolved_species.select("species").unique().collect().height
            ],
            "seqs_in_unresolved_species": [unresolved_species_seqs],
            "unresolved_bin_species": [
                unresolved_bin_species.select("species").unique().collect().height
            ],
            "seqs_in_unresolved_bin_species": [unresolved_bin_species_seqs],
            "unresolved_ambiguous_species": [
                unresolved_ambig_species.select("species").unique().collect().height
            ],
            "seqs_in_unresolved_ambiguous_species": [unresolved_ambig_species_seqs],
            "unresolved_ambiguous_bin_species": [
                unresolved_ambig_bin_species.select("species").unique().collect().height
            ],
            "seqs_in_unresolved_ambiguous_bin_species": [
                unresolved_ambig_bin_species_seqs
            ],
        }
    )
    phylum_counts = (
        bins_per_phyla.join(seqs_per_phylum, on="phylum", how="full", coalesce=True)
        .with_columns(rank=pl.lit("phylum"))
        .fill_null(0)
        .rename({"phylum": "taxa"})
    )
    kingdom_counts = (
        bins_per_kingdom.join(seqs_per_kingdom, on="kingdom", how="full", coalesce=True)
        .with_columns(rank=pl.lit("kingdom"))
        .fill_null(0)
        .rename({"kingdom": "taxa"})
    )
    pl.concat([kingdom_counts, phylum_counts]).write_csv(
        args.taxa_stats_out, separator="\t"
    )
    general_stats.write_csv(args.general_stats_out, separator="\t")
