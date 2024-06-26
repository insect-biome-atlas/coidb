#!/usr/bin/env python

import pandas as pd
import numpy as np
from tqdm import tqdm
import sys
import shutil
import os
import urllib.request
import datetime
import json
import re


def logg(f):
    """
    Decorator for dataframe processing
    :param f:
    :return:
    """

    def wrapper(dataf, *args, **kwargs):
        """
        This wrapper outputs stats on running dataframe processing functions

        :param dataf: input dataframe
        :param args: arguments
        :param kwargs: keyword arguments
        :return: processed dataframes
        """
        # Get rows before processing
        rows_before = dataf.shape[0]
        # Get start time
        tic = datetime.datetime.now()
        # Perform the processing
        result = f(dataf, *args, **kwargs)
        # Get end time
        toc = datetime.datetime.now()
        # Get rows after processing
        rows_after = result.shape[0]
        sys.stderr.write(
            f"{toc-tic} seconds to {f.__name__}, {rows_before-rows_after} rows removed, {rows_after} rows remaining\n"
        )
        return result

    return wrapper


def api_match_species(bin_id):
    """
    This function uses the GBIF API to match a BIN id to a species name.
    It's only used as a fallback in case the add_species function does
    not work, for instance due to improper taxonomic info added for a BIN.

    :param bin_id: The BOLD bin id to match
    :return: Species name or NaN
    """
    try:
        with urllib.request.urlopen(
            f"https://api.gbif.org/v1/species/match?name={bin_id}"
        ) as response:
            json_text = response.read()
        response_dict = json.loads(json_text)
        return response_dict["species"]
    except:
        return np.nan


def add_species(species_bins, bin_tax_df, parent_df):
    """
    This function goes through putative species BINs, i.e. BINs that are
    classified down to genus level and which may have species annotations also.
    These have to be identified by looking up the name of the parent for each BIN.

    :param species_bins: list of BOLD BIN ids that may have species assignments
    :param bin_tax_df: the taxonomic dataframe
    :param parent_df: dataframe with parents to the species bins
    :return: the taxonmomic dataframe with species names added
    """
    for bold_id in tqdm(species_bins, unit=" BINs", desc="adding species"):
        parent = bin_tax_df.loc[bold_id, "parentNameUsageID"]
        if len(parent_df.loc[parent_df.taxonID == parent]) == 0:
            continue
        try:
            # Try to get the species name from the parent backbone
            bin_tax_df.loc[bold_id, "species"] = parent_df.loc[
                parent_df.taxonID == parent, "canonicalName"
            ].values[0]
        except IndexError:
            # If that fails for some reason, do one attempt with the API
            bin_tax_df.loc[bold_id, "species"] = api_match_species(bold_id)
    return bin_tax_df


def find_non_unique_lineages(dataf, ranks):
    """
    This function identifies rank/taxa combinations where parent lineages
    are not unique.

    :param dataf: record dataframe
    :param ranks: ranks to use
    :return: list of rank:name items
    """
    bin_tax_profiles = dataf.reset_index().groupby(ranks).size().reset_index()
    tax_strings = {}
    dups = []
    # Iterate middle ranks
    sys.stderr.write("Looking for non-unique lineages in ranks\n")
    for rank in ranks[1:-1]:
        tax_strings[rank] = {}
        for taxa in tqdm(
            bin_tax_profiles[rank].unique(),
            total=len(bin_tax_profiles[rank].unique()),
            desc=f"{rank}: ",
        ):
            try:
                tax_strings[rank][taxa]
            except KeyError:
                tax_strings[rank][taxa] = []
            # Get all rows with taxa
            rows = bin_tax_profiles.loc[bin_tax_profiles[rank] == taxa]
            parent_ranks = ranks[0 : ranks.index(rank)]
            for row in rows.iterrows():
                tax_strings[rank][taxa].append(
                    ";".join(row[1][p] for p in parent_ranks)
                )
            if len(set(tax_strings[rank][taxa])) > 1:
                dups.append(f"{rank}:{taxa}")
    sys.stderr.write(
        f"Found {len(dups)} rank/taxa combinations with non-unique lineages\n"
    )
    return dups


def check_uniqueness(df, bin_df, group_ranks, rank, name):
    """
    Check whether uniqueness can be achieved by removing bins

    :param df:
    :param bin_df:
    :param group_ranks:
    :return:
    """
    regex = re.compile(".+(_[X]+)$")
    _bins_to_remove = []
    for row in df.groupby(group_ranks).size().reset_index().iterrows():
        # Check whether parent ranks have '_X' in them
        unassigned = sum(
            [1 if x else 0 for x in [regex.match(row[1][rank]) for rank in group_ranks]]
        )
        # if all parent ranks (up to kingdom) are unassigned, mark this bin for removal
        if unassigned == len(group_ranks) - 1:
            # Mark BINs with lineage to be deleted
            _ = df.copy()
            for r in group_ranks:
                _ = _.loc[_[r] == row[1][r]]
            _bins_to_remove += list(_.loc[_[rank] == name].bold_id.values)
    return _bins_to_remove


def prefix_taxa(dataf, d, rank, name, group_ranks, parent_rank, child_ranks):
    """
    Prefixes non-unique taxa labels with the parent rank, e.g.
    Plantae 	Rhodophyta 	Florideophyceae 	Gigartinales 	Acrotylaceae 	Acrotylus 	Acrotylus_X
    Animalia 	Arthropoda 	Insecta 	        Orthoptera 	    Acrididae 	    Acrotylus 	Acrotylus_X

    becomes
    Plantae 	Rhodophyta 	Florideophyceae 	Gigartinales 	Acrotylaceae 	Acrotylaceae_Acrotylus 	Acrotylaceae_Acrotylus_X
    Animalia 	Arthropoda 	Insecta 	        Orthoptera 	    Acrididae 	    Acrididae_Acrotylus 	Acrididae_Acrotylus_X


    :param dataf:
    :param d:
    :param rank:
    :param name:
    :param group_ranks:
    :param parent_rank:
    :param child_ranks:
    :return:
    """
    for key, row_dict in d.items():
        parent = row_dict[parent_rank]
        for child_rank in child_ranks:
            row_dict[child_rank] = row_dict[child_rank].replace(
                name, f"{parent}_{name}"
            )
        d[key] = row_dict
    dataf = pd.concat([dataf.drop(d.keys()), pd.DataFrame(d).T])
    return dataf


def clean_up_non_unique_lineages(dataf, dups, ranks):
    """
    This function iterates the duplicated ranks/names and attempts to
    identify BINs that can be removed  in order to make the
    dataframe unique for parent lineages. If BINs cannot be removed, the
    taxa are instead prefixed with the parent rank

    As an example, the genus Aphaenogaster can be present for BINs like this:
    Animalia 	Animalia_X 	Animalia_XX 	Animalia_XXX 	Animalia_XXXX 	Aphaenogaster
    Animalia 	Arthropoda 	Insecta 	    Hymenoptera 	Formicidae 	    Aphaenogaster

    This function will identify BINs assigned according to the first row, and mark
    them for removal, while keeping BINs assigned as in the second row.

    :param dataf: Record dataframe
    :param dups: Rank/names causing non-uniqueness in dataframe
    :param ranks: Ranks used in dataframe
    :return: Return dataframe cleaned of BINs identified for removal
    """
    bins_to_remove = []
    log = {}
    # items in dups have the format 'rank:taxname', e.g. 'genus:Peyssonnelia'
    for item in dups:
        rank, name = item.split(":")
        log[name] = {"rank": rank}
        group_ranks = ranks[: ranks.index(rank)]
        parent_rank = ranks[ranks.index(rank) - 1]
        child_ranks = ranks[ranks.index(rank) :]
        _df = dataf.loc[dataf[rank] == name].copy()
        _bins_to_remove = check_uniqueness(_df, dataf, group_ranks, rank, name)
        # Test if lineages are unique after removing BINs
        if (
            _df.loc[~_df.bold_id.isin(_bins_to_remove)]
            .groupby(group_ranks)
            .size()
            .reset_index()
            .shape[0]
            == 1
        ):
            # if removing these BINs was enough to generate a single lineage
            # add bins to list to remove
            bins_to_remove += _bins_to_remove
            log[name]["decision"] = "removing BINs with unassigned parent ranks"
            log[name]["bins_removed"] = ",".join(_bins_to_remove)
        else:
            # if not, prefix the name with the parent rank
            d = _df.to_dict(orient="index")
            dataf = prefix_taxa(
                dataf, d, rank, name, group_ranks, parent_rank, child_ranks
            )
            log[name]["decision"] = "prefixing with parent rank name"
            log[name]["bins_removed"] = ""
    return dataf.loc[~dataf["bold_id"].isin(bins_to_remove)], log


def fill_unassigned(
    df,
    bins,
    ranks=["kingdom", "phylum", "class", "order", "family", "genus", "species"],
):
    """
    This function fills in 'the blanks' for each row in a dataframe
    For example:
    BOLD:ACQ8069 	Animalia 	Mollusca 	Gastropoda 	NaN 	Hermaeidae 	Mourgona 	NaN
    BOLD:AAN1572 	Animalia 	Mollusca 	Gastropoda 	NaN 	Hermaeidae 	NaN 	    NaN

    becomes:
    BOLD:ACQ8069 	Animalia 	Mollusca 	Gastropoda 	Gastropoda_X 	Hermaeidae 	Mourgona 	Mourgona_X
    BOLD:AAN1572 	Animalia 	Mollusca 	Gastropoda 	Gastropoda_X 	Hermaeidae 	Hermaeidae_X 	Hermaeidae_XX

    :param df: Dataframe to fill
    :param bins: Bins with NaN ranks
    :param ranks: Ranks to iterate through
    :return: A dataframe with filled ranks
    """
    d = {}
    # Drop the bins to check into separate dataframe
    others = df.drop(bins).loc[:, ranks]
    # Extract bins to check for iteration
    df = df.loc[bins, ranks]
    for bold_bin in tqdm(
        bins, unit=" BINs", total=len(bins), desc="filling unassigned ranks"
    ):
        unknowns = 0
        row = df.loc[bold_bin].to_dict()
        for rank in ranks:
            # If an NaN entry is found
            if row[rank] != row[rank]:
                # Get the previous rank and its classification
                prev_rank = ranks.index(rank) - 1
                prev = ranks[prev_rank]
                # If this is the first time we're adding a suffix include the underscore
                if unknowns == 0:
                    row[rank] = row[prev] + "_X"
                else:
                    row[rank] = row[prev] + "X"
                unknowns += 1
            else:
                # Reset the unknowns counter in case there are intermediate unassigned ranks
                unknowns = 0
        d[bold_bin] = row
    return pd.concat([pd.DataFrame(d).T, others])


@logg
def start(dataf):
    """
    Dummy function to generate a copy of the starting dataframe

    :param dataf: Starting dataframe
    :return: copy of dataframe
    """
    return dataf.copy()


@logg
def extract_bold_bins(dataf):
    """
    Returns the dataframe filtered to only rows with a BOLD BIN id, i.e.
    excluding NaN values

    :param dataf: Input dataframe
    :return: Filtered dataframe
    """
    return dataf.loc[dataf.bold_id == dataf.bold_id]


@logg
def fillna(dataf):
    """
    Fills NaN values in dataframe and returns it

    :param dataf: Input dataframe
    :return: Processed dataframe
    """
    return dataf.fillna("")


@logg
def filter_dataframe(dataf, filter_vals=[], filter_col=None):
    """
    Filter input dataframe based on column->value settings

    :param dataf: Input dataframe
    :param filter_vals: Value to filter with
    :param filter_col: Column in dataframe to filter on
    :return: Filtered dataframe
    """
    if len(filter_vals) > 0:
        sys.stderr.write(
            f"Filtering dataframe to {len(filter_vals)} items in {filter_col}\n"
        )
        return dataf.loc[dataf[filter_col].isin(filter_vals)]
    return dataf


def write_seqs(seq_df, outfile, tmpfile, ranks):
    """
    Writes sequences to file, with taxonomic info in header

    :param seq_df: Sequence dataframe, including BOLD BIN ids and taxonomic info
    :param outfile: Output file to write to
    :param tmpfile: Temporary file
    :param ranks: list of ranks to write in header
    :return: Dataframe of taxonomic information that remains
    """
    old_index = []
    new_index = []
    # Sort sequences by BOLD IDs
    sys.stderr.write("Sorting sequences by BOLD IDs\n")
    seq_df = seq_df.sort_values("bold_id")
    tmpfile = os.path.expandvars(tmpfile)
    outfile = os.path.abspath(outfile)
    with open(tmpfile, "w") as fhout:
        for r in tqdm(
            seq_df.iterrows(),
            desc=f"Writing sequences to temporary directory",
            unit=" seqs",
        ):
            record_id, row = r
            seq = row["seq"]
            desc = ";".join([row[x] for x in ranks + ["bold_id"]])
            fhout.write(f">{record_id} {desc}\n{seq}\n")
    sys.stderr.write(f"Moving {tmpfile} to {outfile}\n")
    shutil.move(tmpfile, outfile)
    return seq_df.drop("seq", axis=1)


def filter(sm):
    """
    Function to filter the BOLD information from GBIF

    :param sm:
    :return:
    """
    genes = sm.params.genes
    taxa = sm.params.filter_taxa
    filter_rank = sm.params.filter_rank
    ranks = sm.params.ranks
    nrows = sm.params.nrows
    logfile = sm.log[0]
    ####################################
    ### Read and process occurrences ###
    ####################################
    sys.stderr.write(f"Reading {sm.input[0]}\n")
    occurrences = pd.read_csv(
        sm.input[0],
        sep="\t",
        usecols=[0, 37],
        names=["record_id", "bold_id"],
        dtype={"bold_id": str},
        nrows=nrows,
    )
    sys.stderr.write(f"{occurrences.shape[0]} records read\n")
    occurrences = occurrences.pipe(start).pipe(extract_bold_bins).pipe(fillna)
    ##################################
    ### Read and process sequences ###
    ##################################
    sys.stderr.write(f"Reading {sm.input[1]}\n")
    seqs = pd.read_csv(
        sm.input[1],
        header=None,
        sep="\t",
        index_col=0,
        names=["record_id", "gene", "seq"],
        usecols=[0, 1, 2],
        nrows=nrows,
    )
    sys.stderr.write(f"{seqs.shape[0]} records read\n")
    seqs = seqs.pipe(start).pipe(filter_dataframe, genes, "gene")
    ##########################################
    ### Read and process backbone taxonomy ###
    ##########################################
    sys.stderr.write(f"Reading {sm.input[2]}\n")
    backbone = pd.read_csv(
        sm.input[2],
        header=0,
        sep="\t",
        nrows=nrows,
        usecols=[0, 2, 5, 7, 11, 17, 18, 19, 20, 21, 22],
    )
    sys.stderr.write(f"{backbone.shape[0]} lines read\n")
    backbone = backbone.pipe(start).pipe(filter_dataframe, taxa, filter_rank)
    #######################################
    ### Merge occurrences and sequences ###
    #######################################
    sys.stderr.write(
        f"Merging occurrences ({occurrences.shape[0]}) and sequence ({seqs.shape[0]}) records\n"
    )
    seq_df = pd.merge(
        occurrences, seqs, left_on="record_id", right_index=True, how="inner"
    )
    sys.stderr.write(f"{seq_df.shape[0]} records remaining\n")
    ################################
    ### Remove duplicate records ###
    ################################
    sys.stderr.write(f"Removing duplicate records\n")
    seq_df_nr = seq_df.groupby("record_id").first().reset_index()
    sys.stderr.write(
        f"{seq_df.shape[0] - seq_df_nr.shape[0]} rows removed, {seq_df_nr.shape[0]} rows remaining\n"
    )
    ##################################
    ### Assign taxonomy to records ###
    ##################################
    # Create new dataframe with scientific name as index
    bin_tax_df = backbone.set_index("scientificName")
    # Extract only BOLD IDs
    bin_tax_df = bin_tax_df.loc[bin_tax_df.index.str.startswith("BOLD:")]
    # Assign default species column
    bin_tax_df = bin_tax_df.assign(
        species=pd.Series([np.nan] * bin_tax_df.shape[0], index=bin_tax_df.index)
    )
    # Extract BOLD ids assigned down to genus level, putative species bins
    species_bins = list(bin_tax_df.loc[bin_tax_df.genus == bin_tax_df.genus].index)
    sys.stderr.write(f"{len(species_bins)} BINs assigned to genus level\n")
    # Extract ids of parents
    parent_ids = bin_tax_df.loc[species_bins, "parentNameUsageID"].values
    # Extract parent dataframe from parent ids
    parent_df = backbone.loc[
        (backbone.taxonID.isin(parent_ids)) & (backbone.taxonRank == "species")
    ]
    sys.stderr.write("Adding species names\n")
    # Attempt to add species to species_bins using parent dataframe
    bin_tax_df = add_species(species_bins, bin_tax_df, parent_df)
    bin_tax_df = bin_tax_df.loc[:, ranks]
    bins_with_species_names = bin_tax_df.loc[
        bin_tax_df.species == bin_tax_df.species
    ].shape[0]
    unique_species_names = len(
        bin_tax_df.loc[bin_tax_df.species == bin_tax_df.species, "species"].unique()
    )
    sys.stderr.write(
        f"Added {unique_species_names} unique species names to {bins_with_species_names} BINS\n"
    )
    # Fill unassigned ranks
    bins_to_fill = list(
        bin_tax_df.loc[bin_tax_df.loc[:, ranks].isna().sum(axis=1) > 0].index
    )
    sys.stderr.write(f"Filling unassigned ranks for {len(bins_to_fill)} BINs\n")
    bin_tax_df = fill_unassigned(bin_tax_df, bins_to_fill, ranks)
    ###############################################
    ### Identify and remove non-unique lineages ###
    ###############################################
    bin_tax_df.index.name = "bold_id"
    dups = find_non_unique_lineages(bin_tax_df, ranks)
    bin_tax_df_cleaned, loglist = clean_up_non_unique_lineages(
        bin_tax_df.reset_index(), dups, ranks
    )
    bin_tax_df_cleaned.set_index("bold_id", inplace=True)
    # Write loglist to file, showing what, if anything, has been done to the taxa
    logdf = pd.DataFrame(loglist).T
    logdf.index.name = "taxa"
    if len(loglist) > 0:
        sys.stderr.write(f"Writing info on taxa name modifications to {logfile}\n")
        logdf.to_csv(logfile, sep="\t")
    ################################################
    ### Merge BIN taxonomy with record dataframe ###
    ################################################
    sys.stderr.write(f"Merging BIN taxonomy with records\n")
    df = pd.merge(
        seq_df_nr,
        bin_tax_df_cleaned.loc[:, ranks],
        left_on="bold_id",
        how="inner",
        right_index=True,
    )
    df.set_index("record_id", inplace=True)
    sys.stderr.write(f"{df.shape[0]} records remaining\n")
    #####################
    ### Write to file ###
    #####################
    # Write seqs to file
    info_df = write_seqs(df, sm.output.fasta, sm.params.tmpf, sm.params.ranks)
    # Write info to file
    info_df.to_csv(sm.output.info, header=True, index=True, sep="\t")


def clean_fasta(sm):
    """
    Reformats fasta headers to strip vsearch specific strings

    :param sm: snakemake object
    :return:
    """
    from Bio import SeqIO
    import re

    with open(sm.input.fasta, "r") as fhin, open(sm.output.fasta, "w") as fhout:
        for record in SeqIO.parse(fhin, "fasta"):
            desc = (record.description).lstrip("centroid=")
            desc = re.split(";seqs=\d+", desc)[0]
            fhout.write(f">{desc}\n{record.seq}\n")


def format_fasta(sm):
    """
    Format a fasta file into two output files, one for use with the
    assignTaxonomy function and one with addSpecies in DADA2.

    The file for assignTaxonomy should have:
    >Level1;Level2;Level3;Level4;Level5;Level6;
    ACCTAGAAAGTCGTAGATCGAAGTTGAAGCATCGCCCGATGATCGTCTGAAGCTGTAGCATGAGTC
    >Level1;Level2;Level3;Level4;Level5;
    CGCTAGAAAGTCGTAGAAGGCTCGGAGGTTTGAAGCATCGCCCGATGGGATCTCGTTGCTGTAGCA

    while assignSpecies should have:
    >ID Genus species
    ACCTAGAAAGTCGTAGATCGAAGTTGAAGCATCGCCCGATGATCGTCTGAAGCTGTAGCATGAGTC
    >ID Genus species
    CGCTAGAAAGTCGTAGAAGGCTCGGAGGTTTGAAGCATCGCCCGATGGGATCTCGTTGCTGTAGCA

    for example:
    >Arthropoda;Insecta;Lepidoptera;Depressariidae;Acraephnes;
    and
    >AANIC216-11 Acraephnes nivea
    :param sm: snakemake object
    :return:
    """
    ranks = sm.params.ranks
    # if "species" in ranks:
    #    ranks.remove("species")
    from Bio import SeqIO

    info = pd.read_csv(sm.input.info, sep="\t", index_col=0, header=0)
    with open(sm.output.assignTaxaFasta, "w") as fh1, open(
        sm.output.addSpeciesFasta, "w"
    ) as fh2:
        for record in SeqIO.parse(sm.input.fasta, "fasta"):
            names = []
            id = record.id
            rec_info = info.loc[id.lstrip("centroid=")]
            # Iterate ranks and add names as long as they are not NaN
            for r in ranks:
                n = rec_info[r]
                if n == n:
                    names.append(n)
                else:
                    break
            id_tax = ";".join(names) + ";"
            fh1.write(f">{id_tax}\n{record.seq}\n")
            species = rec_info["species"]
            bold_id = rec_info["bold_id"]
            id_spec = f"{id} {species}"
            fh2.write(f">{id_spec}\n{record.seq}\n")


def main(sm):
    toolbox = {
        "filter_data": filter,
        "format_dada2": format_fasta,
        "clean": clean_fasta,
    }
    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)  # noqa: F821
