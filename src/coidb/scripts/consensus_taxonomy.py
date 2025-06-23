#!/usr/bin/env python

from multiprocessing import Pool
import polars as pl
import sys
import time
from argparse import ArgumentParser
import re

r = re.compile(r"_X+$")


def cons(df, group_by, threshold):
    _df = (
        df.group_by(group_by)
        .agg(pl.col("n").sum())
        .with_columns((pl.col("n") / pl.col("n").sum() * 100).alias("percent"))
        .filter(pl.col("percent") >= threshold)
    )
    return _df


def calc_unassigned(row):
    return sum([True if r.search(x) else False for x in row])


def calculate_consensus(
    df,
    ranks=["kingdom", "phylum", "class", "order", "family", "genus", "species"],
    threshold=80,
    method="full",
):
    """
    Calculate consensus taxonomy for a dataframe. The dataframe is expected to
    have the following format:

    kingdom	    phylum	        class	    order	        family	        genus	    species	        bin_uri	        n	bin_rows
    "Animalia"	"Arthropoda"	"Insecta"	"Lepidoptera"	"Oecophoridae"	"Garrha"	"Garrha carnea"	"BOLD:AGS2783"	41	3
    "Animalia"	"Arthropoda"	"Insecta"	"Lepidoptera"	"Oecophoridae"	"Garrha"	"Garrha_X"	    "BOLD:AGS2783"	7	3
    "Animalia"	"Arthropoda"	"Insecta"	"Lepidoptera"	"Oecophoridae"	"Garrha"	"Garrha sp."	"BOLD:AGS2783"	2	3

    The function iterates over the length of the ranks list, the first time this
    list contains all ranks: ['kingdom', 'phylum', 'class', 'order', 'family',
    'genus', 'species'] and the number of records for each lineage using these
    ranks is calculated by summing the 'n' column. A percentage is calculated
    and the dataframe is filtered to only include rows with at least
    <threshold>%. If the filtered dataframe has only 1 row, this consensus
    taxonomy is assigned for the BOLD BIN, otherwise the for loop continues by
    taking away one rank, i.e.: ['kingdom', 'phylum', 'class', 'order',
    'family', 'genus'] then attemps to calculate the consensus again. If the
    consensus taxonomy could only be assigned at ranks above species, the
    remaining ranks are given the 'unresolved.' prefix.
    """
    # Store the BOLD BIN id
    bin_uri = df["bin_uri"].unique()[0]
    # Iterate over the list of ranks, successively calculating a consensus using
    # a shorter list
    for i in range(len(ranks), 0, -1):
        # Group by the list of ranks, sum up the 'n' column and calculate a
        # percentage. Then filter the percent column using the threshold.
        if method == "full":
            group_by = ranks[:i]
        elif method == "rank":
            group_by = ranks[i - 1]
        _df = cons(df, group_by, threshold)
        # If the resulting dataframe contains only 1 row, i.e. if only one
        # lineage is above the threshold, use this lineage for the BOLD BIN.
        if _df.shape[0] == 1:
            break
    # Get the most resolved rank
    last_known_rank = [x for x in ranks if x in _df.columns][-1]
    # Get the name of the taxa of the most resolved rank
    last_known = _df[last_known_rank].item()
    # Get any unresolved ranks
    unresolved_ranks = ranks[ranks.index(last_known_rank) + 1 :]
    resolved_ranks = [x for x in ranks if not x in unresolved_ranks]
    # Create a dictionary of the unresolved ranks, all labeled with the most
    # resolved taxonomic name and the 'unresolved.' prefix
    d = dict(
        zip(unresolved_ranks, [f"unresolved.{last_known}"] * len(unresolved_ranks))
    )
    if method == "rank":
        # Get rows matching last_known taxa at last_known_rank
        rows_at_rank = df.filter(pl.col(last_known_rank) == last_known)
        if rows_at_rank.shape[0] > 1:
            # Calculate number of unassigned taxlabels per row
            un = (
                rows_at_rank.select(resolved_ranks)
                .map_rows(calc_unassigned)
                .rename({"map": "unassigned"})
            )
            rows_at_rank = pl.concat([rows_at_rank, un], how="horizontal")
            # Sort by number of unassigned labels in ascending order, then by
            # number of records in descending order and take the first row
            _df = (
                rows_at_rank.sort(["unassigned", "n"], descending=[False, True])
                .head(1)
                .select(resolved_ranks)
            )
        else:
            _df = rows_at_rank.select(resolved_ranks)
    # Add the BOLD BIN id in the 'bin_uri' column
    d["bin_uri"] = bin_uri
    # Create a dataframe by concatenating the resolved dataframe with the
    # unresolved dictionary.
    consensus_df = pl.concat(
        [
            _df.select(resolved_ranks),
            pl.DataFrame(d).select(unresolved_ranks + ["bin_uri"]),
        ],
        how="horizontal",
    )
    return consensus_df


def worker(arg):
    """
    Helper function allowing more than 1 argument to be passed to the
    calculate_consensus function.
    """
    df, ranks, threshold, method = arg
    return calculate_consensus(df=df, ranks=ranks, threshold=threshold, method=method)


def load_taxonomy(infile):
    """
    Load data with lazy API, group by taxonomic ranks + BOLD BIN and count
    unique lineages per BIN. Then sort dataframe by number of records with each
    lineage per BIN.
    """
    bin_taxonomy = (
        pl.scan_csv(infile, separator="\t")
        .group_by(
            [
                "kingdom",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
                "bin_uri",
            ]
        )
        .len("n")
        .sort(["bin_uri", "n"], descending=True)
    ).collect()
    sys.stderr.write(f"Read {bin_taxonomy["bin_uri"].n_unique()} BOLD BINs\n")
    # Create a dataframe for BOLD BINs with more than 1 unique lineage
    bin_taxonomy_ambig = bin_taxonomy.join(
        bin_taxonomy.group_by("bin_uri").len("bin_rows"), on="bin_uri"
    ).filter(pl.col("bin_rows") > 1)
    sys.stderr.write(
        f"Found {bin_taxonomy_ambig["bin_uri"].n_unique()} BOLD BINs with non-unique lineages\n"
    )
    # Create a dataframe for BOLD BINs with only 1 lineage
    bin_taxonomy_unambig = (
        bin_taxonomy.join(
            bin_taxonomy.group_by("bin_uri").len("bin_rows"), on="bin_uri"
        )
        .filter(pl.col("bin_rows") < 2)
        .drop(["n", "bin_rows"])
    )
    sys.stderr.write(
        f"Found {bin_taxonomy_unambig["bin_uri"].n_unique()} BOLD BINs with unique lineages\n"
    )
    return bin_taxonomy, bin_taxonomy_ambig, bin_taxonomy_unambig


def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str, help="Input TSV file")
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        help="Output file with consensus taxonomies",
        required=True,
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=int,
        default=80,
        help="Consensus threshold (in %%) for assigning a taxonomic label at a rank",
    )
    parser.add_argument(
        "-r",
        "--ranks",
        nargs="+",
        type=str,
        choices=["kingdom", "phylum", "class", "order", "family", "genus", "species"],
        default=["kingdom", "phylum", "class", "order", "family", "genus", "species"],
        help="Ranks to include in output",
    )
    parser.add_argument(
        "-m",
        "--method",
        type=str,
        choices=["rank", "full"],
        default="full",
        help="Method to calculate consensus. If 'full' the number of records is summed to the full lineage which takes into account if higher level ranks differ for each BOLD BIN. If 'rank' the records are summed per rank irrespective of higher ranks, then the full lineage is set to the most frequent lineage with fewest unassigned taxlabels.",
    )
    parser.add_argument(
        "-p", dest="cpus", type=int, default=1, help="Number of cpus to use"
    )
    args = parser.parse_args()

    sys.stderr.write(f"Loading data from {args.infile}\n")
    bin_taxonomy, bin_taxonomy_ambig, bin_taxonomy_unambig = load_taxonomy(args.infile)
    # Partition the ambiguous BOLD BINs into a list of dataframes
    dataframes = bin_taxonomy_ambig.partition_by("bin_uri")
    # Apply the consensus function to each dataframe in the dataframes list
    # Use the worker function to supply the threshold as an argument
    sys.stderr.write(
        f"Calculating consensus taxonomies for non-unique BINS using {args.cpus} cpus\n"
    )
    with Pool(args.cpus) as p:
        consensus_list = p.map(
            worker, ((df, args.ranks, args.threshold, args.method) for df in dataframes)
        )
    sys.stderr.write(f"Concatenating results\n")
    # Concatenate the dataframes in the consensus_list + the unambiguous dataframe
    consensus_df = pl.concat(
        consensus_list + [bin_taxonomy_unambig.select(args.ranks + ["bin_uri"])]
    )
    consensus_df.write_csv(args.outfile, separator="\t")
