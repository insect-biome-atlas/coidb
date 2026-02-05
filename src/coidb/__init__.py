# SPDX-FileCopyrightText: 2025-present John Sundh <john.sundh@scilifelab.se>
#
# SPDX-License-Identifier: MIT

import gzip as gz
from tempfile import NamedTemporaryFile

__version__ = "0.6.0"


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
