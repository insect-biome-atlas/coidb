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
    """
    non_unique = {}
    for rank in ranks:
        sys.stderr.write(f"{rank}...\n")
        i = ranks.index(rank) + 1
        _ranks = ranks[0:i]
        q = df.with_columns(lineage=pl.concat_str(_ranks, separator=";")).select(
            rank, "lineage"
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


def check_parent_ranks(_df, group_ranks, regex):
    """
    This function checks the taxonomic labels of <group_ranks> for each row and
    counts number of ranks that end in '_X' (any number of Xs) using a regex
    pattern. If all ranks match the regex, the BOLD bin is added to the
    bins_to_remove list which is finally returned.
    """
    bins_to_remove = []
    rows = _df.select(["bin_uri"] + group_ranks)
    for row in rows.iter_rows():
        unassigned = sum([1 if regex.match(x) else 0 for x in row[2:]])
        if unassigned == len(group_ranks) - 1:
            bins_to_remove.append(row[0])
    return bins_to_remove


def fix_non_unique_lineages(df, non_unique, ranks):
    """
    This function iterates the duplicated ranks/names and attempts to
    identify BINs that can be removed in order to make the
    dataframe unique for parent lineages. If BINs cannot be removed, the
    taxa are instead prefixed with the parent rank

    As an example, the genus Aphaenogaster can be present for BINs like this:
    kingdom  phylum     class       order         family        genus
    Animalia Animalia_X Animalia_XX Animalia_XXX  Animalia_XXXX Aphaenogaster
    Animalia Arthropoda Insecta 	Hymenoptera   Formicidae 	Aphaenogaster

    This function will identify BINs assigned according to the first row, and mark
    them for removal, while keeping BINs assigned as in the second row.

    If removal of BINs assigned as in the first row is not enough to generate a
    unique lineage, then the conflicting taxlabels are prefixed with their parent taxa.
    """
    regex = re.compile(".+(_[X]+)$")
    bins_to_remove = []
    for rank in ranks:
        try:
            taxa = non_unique[rank]
        except KeyError:
            continue
        for t in taxa:
            group_ranks = ranks[: ranks.index(rank)]
            parent_rank = ranks[ranks.index(rank) - 1]
            _df = df.filter(pl.col(rank) == t)
            _bins_to_remove = check_parent_ranks(_df, group_ranks, regex)
            if (
                _df.filter(~pl.col("bin_uri").is_in(_bins_to_remove))
                .group_by(group_ranks)
                .len()
                .height
                == 1
            ):
                bins_to_remove += _bins_to_remove
                sys.stderr.write(
                    f"Removing {len(_bins_to_remove)} BINs for {rank}:{t}\n"
                )
            else:
                sys.stderr.write(f"Prefixing {rank}:{t} with {parent_rank}\n")
                df = df.with_columns(
                    pl.when(pl.col(rank) == t)
                    .then(pl.concat_str([parent_rank, rank], separator="_"))
                    .otherwise(pl.col(rank))
                    .alias(rank)
                )
    return df.filter(~pl.col("bin_uri").is_in(bins_to_remove))


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
    args = parser.parse_args()
    df = pl.scan_csv(args.infile, separator="\t")
    sys.stderr.write(f"Finding non-unique lineages in {args.infile}\n")
    non_unique = find_non_unique(df, args.ranks)
    dups = df.filter(
        (pl.col("kingdom").is_in(non_unique["kingdom"]))
        | (pl.col("phylum").is_in(non_unique["phylum"]))
        | (pl.col("class").is_in(non_unique["class"]))
        | (pl.col("order").is_in(non_unique["order"]))
        | (pl.col("family").is_in(non_unique["family"]))
        | (pl.col("genus").is_in(non_unique["genus"]))
        | (pl.col("species").is_in(non_unique["species"]))
    ).collect()
    sys.stderr.write(
        "Non-unique taxa per rank:"
        + " ".join([f"{rank}:{len(non_unique[rank])}" for rank in args.ranks])
        + "\n"
    )
    sys.stderr.write(f"{dups.height} records with non-unique lineages\n")
    sys.stderr.write("Fixing non-unique lineages\n")
    unique_df = fix_non_unique_lineages(dups, non_unique, args.ranks)
    pl.concat(
        [
            df.filter(
                ~pl.col("processid").is_in(
                    unique_df.select("processid").to_series().to_list()
                )
            ),
            unique_df.lazy(),
        ]
    ).sink_csv(args.outfile, separator="\t", engine="streaming")
