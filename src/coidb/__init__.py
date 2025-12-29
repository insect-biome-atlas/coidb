# SPDX-FileCopyrightText: 2025-present John Sundh <john.sundh@scilifelab.se>
#
# SPDX-License-Identifier: MIT

import gzip as gz
from tempfile import NamedTemporaryFile


def read_records(f):
    """
    Read records from fasta file
    """
    records = []
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
                records.append(line.lstrip(">").split(" ")[0])
    return records


def series_to_fasta(series, outfile, compress=False):
    """
    Writes a polars Series to a file in fasta format
    """
    if compress:
        with gz.open(outfile, "wt") as fhout:
            fhout.write("\n".join(series.to_list()))
    else:
        with open(outfile, "w") as fhout:
            fhout.write("\n".join(series.to_list()))


def get_header(f):
    """
    Returns header from input file
    """
    with open(f, "r") as fhin:
        return fhin.readline().rstrip().split("\t")


def extract_columns(f, indices):
    """
    Loops over each line in input file and only writes certain columns to a
    temporary file.
    """
    fout = NamedTemporaryFile(mode="w", delete=False)
    records = {}
    with open(f, "r") as fhin:
        for line in fhin:
            line = line.rstrip()
            items = line.split("\t")
            try:
                records[items[0]]
                continue
            except KeyError:
                records[items[0]] = 1
            kept_items = []
            for i in indices:
                try:
                    kept_items.append(items[i])
                except IndexError:
                    kept_items.append("")
            fout.write("\t".join(kept_items) + "\n")
    return fout.name
