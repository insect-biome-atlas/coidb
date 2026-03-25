#!/usr/bin/env python

import polars as pl
from argparse import ArgumentParser
import sys


def closest(lst, K):
    """
    Returns the closest value to K from lst
    """
    return min(lst, key=lambda x: abs(x - K))


def polars_to_string(df, header=""):
    s = str(df)
    items = s.split("\n")
    if header:
        items = [header] + items[1:]
    return "\n".join(items)


def make_lineage_col(
    df, ranks=["kingdom", "phylum", "class", "order", "family", "genus", "species"]
):
    return df.insert_column(
        0,
        df.fill_null("NA")
        .with_columns(lineage=pl.concat_str(ranks, separator=";"))
        .select("lineage")
        .to_series(),
    )


def fix_reptilia(df):
    """
    Set correct class and phylum for species assigned to Reptilia
    """
    cols = df.collect_schema().names()
    return (
        df.with_columns(
            # species with class assignment of Squamata, Testudines, Crocodyla and Sphenodontia
            # are incorrect as these correspond to orders within class Reptilia
            pl.when(
                pl.col("class").is_in(
                    ["Squamata", "Testudines", "Crocodylia", "Sphenodontia"]
                )
            )
            # for these, leave kingdom, family, genus and species as is
            # but set order to class, set class to 'Reptilia' and phylum to 'Chordata'
            .then(
                pl.struct(
                    kingdom="kingdom",
                    phylum=pl.lit("Chordata"),
                    Class=pl.lit("Reptilia"),
                    order=pl.col("class"),
                    family="family",
                    genus="genus",
                    species="species",
                )
            )
            # for all other species keep assignments as is
            .otherwise(pl.struct(cols)).struct.unnest()
        )
        # because we can't use the class as a variable in the struct above
        # set Class column to class value if null
        .with_columns(
            Class=pl.when(pl.col("Class").is_not_null())
            .then(pl.col("Class"))
            .otherwise(pl.col("class"))
        )
        # then drop class column and rename Class to class
        .drop("class")
        .rename({"Class": "class"})
        .select(cols)
    )


def fix_symbionts(
    df, ranks=["kingdom", "phylum", "class", "order", "family", "genus", "species"]
):
    """
    Fixes potential symbiont cases by prefixing genus and species labels
    """
    return df.with_columns(
        pl.when(
            (pl.col("kingdom") == "Bacteria")
            & (
                ~pl.col("kingdom_right").is_in(
                    ["Bacteria", "Archaea", "incertae sedis"]
                )
            )
            & (pl.col("species_right").is_not_null())
        )
        .then(
            pl.struct(
                kingdom_right=pl.col("kingdom"),
                phylum_right=pl.col("phylum"),
                class_right=pl.col("class"),
                order_right=pl.col("order"),
                family_right=pl.col("family"),
                genus_right=pl.concat_str(
                    [pl.lit("symbiont of "), pl.col("genus_right")]
                ),
                species_right=pl.concat_str(
                    [pl.lit("symbiont of "), pl.col("species_right")]
                ),
            )
        )
        .otherwise(pl.struct([f"{rank}_right" for rank in ranks]))
        .struct.unnest()
    )


def rename_protista(df):
    """
    Renames the Protista kingdom based on phylum assignment
    """
    return df.with_columns(
        kingdom=pl.when(
            pl.col("phylum").is_in(
                [
                    "Ochrophyta",
                    "Bacillariophyta",
                    "Heterokontophyta",
                    "Haptophyta",
                    "Cryptophyta",
                    "Ciliophora",
                    "Pyrrophycophyta",
                    "Apicomplexa",
                ]
            )
        )
        .then(pl.lit("Chromista"))
        .when(pl.col("phylum").is_in(["Rhodophyta", "Glaucophyta"]))
        .then(pl.lit("Plantae"))
        .when(pl.col("phylum").is_in(["Euglenida", "Amoebozoa"]))
        .then(pl.lit("Protozoa"))
        .otherwise(pl.col("kingdom"))
    )


def gap_filling(
    partially_matched,
    ranks=["kingdom", "phylum", "class", "order", "family", "genus", "species"],
):
    """
    # for each row,
    # if a taxon name is missing in GBIF but not in BOLD:
    #   if the next higher rank that is not empty in either GBIF or BOLD has the same content:
    #       import it from BOLD
    #   else:
    #       leave the higher taxon name empty.
    # if a taxon name is missing in BOLD but not in GBIF:
    #
    # if a taxon name is present in GBIF and in BOLD:
    #   use the GBIF taxon name
    """
    cons = {}
    errors = {}
    for rank in ranks:
        cons[rank] = []
        cons[f"{rank}_right"] = []
    # for each row in partially_matched
    for row in partially_matched.iter_rows():
        index = row[0]
        # cons["i"].append(index)
        # take the lineage from BOLD
        lineage = [
            row[i] for i in [partially_matched.columns.index(rank) for rank in ranks]
        ]
        # take the matched taxonomy from GBIF
        matched = [
            row[i]
            for i in [
                partially_matched.columns.index(f"{rank}_right") for rank in ranks
            ]
        ]
        # extract all indices with missing data for GBIF and BOLD
        matched_missing = [i for i, val in enumerate(matched) if str(val) == "None"]
        lineage_missing = [i for i, val in enumerate(lineage) if str(val) == "None"]
        # extract all indices with data for GBIF and BOLD
        matched_present = [i for i, val in enumerate(matched) if str(val) != "None"]
        lineage_present = [i for i, val in enumerate(lineage) if str(val) != "None"]
        # the default will be to use BOLD taxonomic labels to fill the holes in the GBIF data
        use_bold = True
        # for each index with missing data in GBIF
        for i in matched_missing:
            # find the nearest higer rank with non-missing taxonomic labels in both BOLD and GBIF
            closest_assigned_parent = closest(
                [
                    x
                    for x in set(matched_present).intersection(lineage_present)
                    if x < 3
                ],
                i,
            )
            # if the nearest higher rank label is not the same in both GBIF and BOLD
            # set use_bold = False
            if matched[closest_assigned_parent] != lineage[closest_assigned_parent]:
                error = f"{ranks[closest_assigned_parent]}:{matched[closest_assigned_parent]} != {lineage[closest_assigned_parent]}"
                try:
                    errors[error].append(row)
                except KeyError:
                    errors[error] = [row]
                use_bold = False
                break
        # if all the assigned parents have the same labels
        # use BOLD to fill in the gaps
        if use_bold:
            for i in matched_missing:
                matched[i] = lineage[i]
        for i, rank in enumerate(ranks):
            cons[rank].append(lineage[i])
            cons[f"{rank}_right"].append(matched[i])
    cons_df = pl.DataFrame(cons).select(ranks + [f"{rank}_right" for rank in ranks])
    errors_df = pl.DataFrame(
        data={
            "error": errors.keys(),
            "count": [len(errors[key]) for key in errors.keys()],
        }
    ).sort("count", descending=True)
    return cons_df, errors_df


def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_taxfile", help="TSV file with taxonomic labels")
    parser.add_argument(
        "-m",
        "--matched",
        help="TSV file with taxonomic labels and a 'name' column corresponding to the original species name",
    )
    args = parser.parse_args()
    ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    sys.stderr.write(f"Reading GBIF matched results from {args.matched}\n")
    matched_df = pl.scan_csv(args.matched, separator="\t")
    sys.stderr.write("Fixing Reptilia class\n")
    matched_df = fix_reptilia(matched_df)
    sys.stderr.write(f"Reading BOLD data from {args.input_taxfile}\n")
    input_df = pl.scan_csv(args.input_taxfile, separator="\t").with_columns(
        pl.col(pl.String).replace("None", None)
    )
    sys.stderr.write("Fixing Protista kingdom\n")
    input_df = rename_protista(input_df)
    sys.stderr.write("Creating unique lineage dataframe\n")
    rank_unique = input_df.group_by(ranks).first().select(ranks).collect()
    sys.stderr.write("BEFORE MATCHING\n")
    x = sys.stderr.write(
        polars_to_string(
            rank_unique.select(pl.all().n_unique()), header="unique taxa per rank:"
        )
        + "\n"
    )
    x = sys.stderr.write(
        polars_to_string(rank_unique.null_count(), header="unassigned per rank:") + "\n"
    )
    # make dataframe with BOLD + GBIF by joining with matched dataframe
    # for rows with assigned species
    to_match = rank_unique.join(
        matched_df.filter(pl.col("species").is_not_null()).collect(),
        left_on="species",
        right_on="name",
    )
    sys.stderr.write("Fixing labels for symbionts\n")
    to_match = fix_symbionts(to_match)
    # Get rows which are fully matched (contain no null values in the matched ranks)
    fully_matched = to_match.drop_nulls([f"{r}_right" for r in ranks])
    sys.stderr.write(
        f"Matched {fully_matched.select("species").n_unique()} species names in full to {fully_matched.select("species_right").n_unique()} species names\n"
    )
    # Get rows for remaining species
    partially_matched = to_match.filter(
        (
            ~pl.col("species_right").is_in(
                fully_matched.select("species_right").to_series().to_list()
            )
        )
    )
    sys.stderr.write(
        f"{partially_matched.select("species").n_unique()} species names matched with partial taxonomy\n"
    )
    sys.stderr.write("Attempting to fill gaps in taxonomy\n")
    cons_df, errors_df = gap_filling(partially_matched)
    consolidated = pl.concat([cons_df, fully_matched])
