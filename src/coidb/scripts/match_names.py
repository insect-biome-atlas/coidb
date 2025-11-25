#!/usr/bin/env python

from pygbif import species
import polars as pl
from multiprocessing import Pool
from argparse import ArgumentParser
from tqdm import tqdm
from collections import defaultdict


def match_species(
    species_name,
    ranks=["kingdom", "phylum", "class", "order", "family", "genus", "species"],
):
    """
    Matches species names to the GBIF backbone using the species.name_backbone
    function from pygbif.

    Parameters
    ----------
    species_name: str
        Species name to be matched
    ranks: list
        Taxonomic ranks to return for matched species

    Return
    ------
    d: dict
        Dictionary with keys 'name': original species name and <rank>: matched rank
        for the match in GBIF.
    """
    d = defaultdict.fromkeys(ranks, None)
    d["name"] = species_name
    results = species.name_backbone(species_name)
    for rank in ranks:
        try:
            d[rank] = results[rank]
        except KeyError:
            continue
    return d


def get_species_names(f):
    info = pl.scan_csv(f, separator="\t")
    sp = (
        (
            info.filter(~pl.col("species").str.contains(r"_X+$"))
            .select("species")
            .unique()
        )
        .collect()
        .to_series()
        .to_list()
    )
    return sp


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "-i", "--infile", type=str, help="TSV infile with a 'species' column"
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        help="TSV outfile with taxonomic ranks matched to backbone",
    )
    parser.add_argument(
        "-p", dest="cpus", type=int, default=1, help="Number of cpus to use"
    )
    args = parser.parse_args()
    species_names = get_species_names(args.infile)
    with Pool(args.cpus) as p:
        matches = pl.DataFrame(
            list(
                tqdm(
                    p.imap_unordered(match_species, species_names),
                    total=len(species_names),
                    unit=" species names",
                    desc="Matching species names",
                )
            )
        )
    matches.write_csv(args.outfile, separator="\t")
