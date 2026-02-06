#!/usr/bin/env python

import polars as pl
from tqdm import tqdm
import gzip as gz
from argparse import ArgumentParser
import sys


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
    parser.add_argument("--tsv", type=str, help="COIDB info TSV file")
    parser.add_argument("--fasta", type=str, help="COIDB clustered fasta.gz file")
    parser.add_argument("--consensus", type=str, help="Consensus taxonomy TSV file")
    args = parser.parse_args()
    tsv = args.tsv
    fasta = args.fasta
    consensus = args.consensus
    (
        seqids,
        bins_df,
        mean_clusters_per_bins,
        median_clusters_per_bins,
        max_clusters_per_bins,
    ) = count_bin_clusters(fasta)
    consensus_df = pl.scan_csv(consensus, separator="\t")
    with gz.open(tsv, "rt") as fhin:
        seqs = -1
        for line in tqdm(fhin, desc=f"Reading {tsv}", unit=" lines"):
            seqs += 1
    total_seqs = seqs
    total_bins = (
        bins_df.filter(pl.col("bin_uri").str.starts_with("BOLD:"))
        .select(pl.len())
        .collect()
        .item()
    )
    seqs_per_kingdom = (
        consensus_df.join(bins_df, on="bin_uri")
        .group_by("kingdom")
        .agg(pl.sum("n"))
        .rename({"n": "n_seqs"})
        .sort("n_seqs", descending=True)
        .collect()
    )
    bins_per_kingdom = (
        consensus_df.filter(pl.col("bin_uri").str.starts_with("BOLD:"))
        .group_by("kingdom")
        .len()
        .sort("len", descending=True)
        .rename({"len": "n_bins"})
        .collect()
    )
    bins_per_phyla = (
        consensus_df.filter(pl.col("bin_uri").str.starts_with("BOLD:"))
        .group_by("phylum")
        .len()
        .sort("len", descending=True)
        .rename({"len": "n_bins"})
        .collect()
    )
    sys.stdout.write(f"Total records: {total_seqs}\n")
    sys.stdout.write(f"Total number of BINs: {total_bins}\n")
    sys.stdout.write(f"Total sequences (clustered): {len(seqids)}\n")
    sys.stdout.write(f"Clustered sequences per BIN:\n")
    sys.stdout.write(f"    mean: {round(mean_clusters_per_bins)}\n")
    sys.stdout.write(f"    median: {round(median_clusters_per_bins)}\n")
    sys.stdout.write(f"    max: {max_clusters_per_bins}\n")
    sys.stdout.write("Total sequences (clustered) per kingdom:\n")
    seqs_per_kingdom.write_csv(sys.stdout, separator="\t")
    sys.stdout.write("\n")
    sys.stdout.write("Total BINs per kingdom:\n")
    bins_per_kingdom.write_csv(sys.stdout, separator="\t")
    sys.stdout.write("\n")
    sys.stdout.write("Total BINs per phylum:\n")
    bins_per_phyla.write_csv(sys.stdout, separator="\t")
    sys.stdout.write("\n")
    total_species = consensus_df.select("species").unique().collect().height
    total_bin_species = (
        consensus_df.filter(pl.col("bin_uri").str.starts_with("BOLD:"))
        .select("species")
        .unique()
        .collect()
        .height
    )
    ambig_species = consensus_df.filter(
        (pl.col("species").str.contains(r"_X+$"))
        & (~pl.col("species").str.starts_with("unresolved"))
    )
    ambig_bin_species = ambig_species.filter(pl.col("bin_uri").str.starts_with("BOLD:"))
    unresolved_species = consensus_df.filter(
        (pl.col("species").str.starts_with("unresolved"))
        & (~pl.col("species").str.contains(r"_X+$"))
    )
    unresolved_bin_species = unresolved_species.filter(
        pl.col("bin_uri").str.starts_with("BOLD:")
    )
    unresolved_ambig_species = consensus_df.filter(
        (pl.col("species").str.starts_with("unresolved"))
        & (pl.col("species").str.contains(r"_X+$"))
    )
    unresolved_ambig_bin_species = unresolved_ambig_species.filter(
        pl.col("bin_uri").str.starts_with("BOLD:")
    )
    sys.stdout.write(f"Total species: {total_species}\n")
    sys.stdout.write(f"Total BIN species: {total_bin_species}\n")
    sys.stdout.write(
        f"Ambiguous species: {ambig_species.select("species").unique().collect().height}\n"
    )
    sys.stdout.write(
        f"Ambiguous BIN species: {ambig_bin_species.select("species").unique().collect().height}\n"
    )
    sys.stdout.write(
        f"Unresolved species: {unresolved_species.select("species").unique().collect().height}\n"
    )
    sys.stdout.write(
        f"Unresolved BIN species: {unresolved_bin_species.select("species").unique().collect().height}\n"
    )
    sys.stdout.write(
        f"Unresolved and ambiguous species: {unresolved_ambig_species.select("species").unique().collect().height}\n"
    )
    sys.stdout.write(
        f"Unresolved and ambiguous BIN species: {unresolved_ambig_bin_species.select("species").unique().collect().height}\n"
    )
