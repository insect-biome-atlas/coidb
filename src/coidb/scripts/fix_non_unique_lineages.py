#!/usr/bin/env python

from argparse import ArgumentParser
import polars as pl
import sys
import re


def find_non_unique(df, ranks):
    """
    This function loops through the ranks and for each rank creates a new
    dataframe with a 'lineage' column which contains all the taxlabels of the
    current + parent ranks concatenated. It then groups by the current rank and
    counts number of unique 'lineages' found for the rank. If there's more than
    1 unique lineage this means that parent ranks have conflicting taxlabels.
    Taxa with conflicting parent ranks are saved to a list in a dictionary and
    returned.

    Parameters
    ----------
    df : polars DataFrame
        Data frame with taxonomic ranks as columns
    ranks : list
        List of ranks for which to search for non-unique lineages

    Returns
    -------
    non_unique: dict
        Dictionary with ranks as keys and a list of taxa names as values
    """
    non_unique = {}
    for rank in ranks:
        sys.stderr.write(f"{rank}...\n")
        i = ranks.index(rank) + 1
        _ranks = ranks[0:i]
        parent_rank = ranks[ranks.index(rank) - 1]
        q = (
            df.filter(~pl.col(parent_rank).str.contains(r"_X+$"))
            .with_columns(lineage=pl.concat_str(_ranks, separator=";"))
            .select(rank, "lineage")
        )
        non_unique[rank] = (
            q.group_by(rank)
            .n_unique()
            .filter(pl.col("lineage") > 1)
            .select(rank)
            .collect()
            .to_series()
            .to_list()
        )
    return non_unique


def find_bins_to_remove(df, group_ranks, id_col="bin_uri"):
    rows = df.select([id_col] + group_ranks).with_columns(
        assigned_ranks=len(group_ranks)
        - pl.sum_horizontal(pl.col((group_ranks)).str.contains(r"_X+$"))
    )
    return (
        rows.filter(pl.col("assigned_ranks") == 0).select(id_col).to_series().to_list()
    )


def fix_non_unique_lineages(df, non_unique, ranks, id_col="bin_uri", remove=False):
    """

    This function takes as input a dataframe of rows with non-unique lineages
    and prefixes the taxonomic labels with the parent rank in order to make the
    lineages unique.

    For example, the genus 'Achlya' can have the lineage strings:
    Animalia;Arthropoda;Insecta;Lepidoptera;Drepanidae;Achlya
    Protista;Oomycota;Peronosporea;Saprolegniales;Saprolegniaceae;Achlya

    This will be modified to:
    Animalia;Arthropoda;Insecta;Lepidoptera;Drepanidae;Drepanidae_Achlya
    Protista;Oomycota;Peronosporea;Saprolegniales;Saprolegniaceae;Saprolegniaceae_Achlya

    Parameters
    ----------
    df : polars DataFrame
        DataFrame of rows with non-unique lineages
    non_unique : dict
        Dictionary with ranks as keys and taxa with non-unique parent lineages
        as values
    ranks : list
        List of taxonomic ranks to iterate
    id_col : string
        Column name to use for removal of non-unique lineages
    remove : boolean
        If True, the function will attempt to fix the non-unique lineages
        by first removing ambiguous assignments
    """
    bins_to_remove = []
    for rank in ranks[1:]:
        try:
            taxa = non_unique[rank]
        except KeyError:
            continue
        for t in taxa:
            group_ranks = ranks[: ranks.index(rank)]
            parent_rank = ranks[ranks.index(rank) - 1]
            _df = df.filter(pl.col(rank) == t)
            if remove:
                _bins_to_remove = find_bins_to_remove(_df, group_ranks, id_col)
                if (
                    _df.filter(~pl.col(id_col).is_in(_bins_to_remove))
                    .unique(group_ranks)
                    .height
                    == 1
                ):
                    bins_to_remove += _bins_to_remove
                    sys.stderr.write(
                        f"Removing {len(_bins_to_remove)} features for {rank}:{t}\n"
                    )
                    continue
            sys.stderr.write(f"Prefixing {rank}:{t} with {parent_rank}\n")
            df = df.with_columns(
                pl.when(pl.col(rank) == t)
                .then(pl.concat_str([parent_rank, rank], separator="_"))
                .otherwise(pl.col(rank))
                .alias(rank)
            )
    return df.filter(~pl.col(id_col).is_in(bins_to_remove))


def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str, help="Input TSV file")
    parser.add_argument("-o", "--outfile", type=str, help="Output TSV file")
    parser.add_argument(
        "-r",
        "--ranks",
        nargs="+",
        default=["kingdom", "phylum", "class", "order", "family", "genus", "species"],
    )
    parser.add_argument(
        "--id_col",
        type=str,
        help="Column name containing identifier to remove if running with --remove",
        default="bin_uri",
    )
    parser.add_argument(
        "--remove",
        action="store_true",
        help="Attempt to remove rows with only missing information for higher ranks in order to make lineages unique.",
    )
    args = parser.parse_args()
    df = pl.scan_csv(args.infile, separator="\t")
    index_key = df.collect_schema().names()[0]
    sys.stderr.write(f"Finding non-unique lineages in {args.infile}\n")
    non_unique = find_non_unique(df, args.ranks)
    # generate a list of Lazy Frames containing all non-unique rows
    dup_list = []
    for rank in args.ranks:
        dup_list.append(df.filter(pl.col(rank).is_in(non_unique[rank])))
    # concatenate and collect into a single DataFrame
    dups = pl.concat(dup_list).collect()
    sys.stderr.write(
        "Non-unique taxa per rank:"
        + " ".join([f"{rank}:{len(non_unique[rank])}" for rank in args.ranks])
        + "\n"
    )
    sys.stderr.write(f"{dups.height} records with non-unique lineages\n")
    sys.stderr.write("Fixing non-unique lineages\n")
    # generate a dataframe with unique lineages
    unique_df = fix_non_unique_lineages(
        dups, non_unique, args.ranks, args.id_col, args.remove
    )
    pl.concat(
        [
            df.filter(
                ~pl.col(index_key).is_in(
                    unique_df.select(index_key).to_series().to_list()
                )
            ),
            unique_df.lazy(),
        ]
    ).sink_csv(args.outfile, separator="\t", engine="streaming")
