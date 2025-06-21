# SPDX-FileCopyrightText: 2025-present John Sundh <john.sundh@scilifelab.se>
#
# SPDX-License-Identifier: MIT

import gzip as gz


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
