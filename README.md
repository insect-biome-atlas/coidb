[![Pixi Badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json)](https://pixi.sh)

# COI DB

## Table of contents

- [Overview](#overview)
- [Installation](#installation)
  - [Install with pixi](#install-with-pixi-from-source)
  - [Install with conda](#install-with-conda)
- [Obtain data](#obtain-data)
- [Running coidb](#running-coidb)
  - [Using a configuration file](#using-a-configuration-file)
  - [Cluster execution](#cluster-execution)
- [Output](#output)
- [How it works](#how-it-works)


## Overview

The coidb package runs a Snakemake workflow under the hood which contains steps
to filter the BOLD public data to only the COI-5P marker gene, remove leading
and trailing gaps and sequences with internal gaps and ambiguous nucleotides. It
also applies a length filtering and only keeps records assigned to a BOLD BIN.
Steps are also taken to ensure that taxonomic lineages are unique by prefixing
duplicated taxonomic labels or removing BOLD BINs with unassigned records. The
filtered sequences are then dereplicated by clustering sequences within each
BOLD BIN using vsearch. A consensus taxonomy is calculated using an 80%
consensus threshold starting from species and moving up in the taxonomy tree. 

Finally, fasta and tab separated files compatible with SINTAX, DADA2 and QIIME2 are generated.

## Installation

### Install with pixi from source

1. First install [pixi](https://pixi.sh/latest/#installation)

```bash
curl -fsSL https://pixi.sh/install.sh | sh
```

2. Then clone the GitHub repository and change into the `coidb` directory

```bash
git clone git@github.com:insect-biome-atlas/coidb.git
cd coidb
```

3. Install the `coidb` package with pixi and start a shell in the installed environment

```bash
pixi shell
```

If the installation worked you should be able to run `coidb run -h`. To test the
installation, run `pixi run test`. Proceed to the [Running
coidb](#running-coidb) section to see usage information.

### Install with Conda

1. Make sure [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) is installed on your system.

2. Create a new environment with `coidb` installed:

```bash
conda create -n coidb -c bioconda coidb
```

3. Activate the `coidb` environment:

```bash
conda activate coidb
```

## Obtain data

The coidb package uses public barcode reference libraries from [BOLD](https://bench.boldsystems.org/index.php) to build reference fasta files compatible with tools like SINTAX, QIIME2 and DADA2. The first thing you need to do is get your hands on a BOLD Public Data Package:

1. Go to the BOLD systems [data package page](https://bench.boldsystems.org/index.php/datapackages/Latest).
2. Login is required to access files on this page, so either login or sign up if you don't already have an account.
3. On the data package page, click the **Data package (tar.gz compressed)** download button, accept the terms and click **Download** to obtain a temporary download link.
4. Use the link to download the data package which will be named `BOLD_Public.<dd>-<Mmm>-<YYYY>.tar.gz`, for example `BOLD_Public.20-Jun-2025.tar.gz`.

> [!TIP]
> To download via the command line you can copy the Download link instead of
> clicking it, then use `wget` or `curl` to download directly to a file of your
> choice. For example, to download the data package to a directory called
> `data/` you could run:
> ```bash
> mkdir data
> wget -O data/BOLD_Public.20-Jun-2025.tar.gz <copied download link>
> ```
> or with `curl`:
> ```bash
> mkdir data 
> curl -o data/BOLD_Public.20-Jun-2025.tar.gz <copied download link>
> ```

The downloaded file can now be used as input to `coidb` by pointing to it with
the `--input-file` argument, or by setting `input_file:
<path-to-downloaded-tar.gz file>` in a [configuration
file](#using-a-configuration-file)

BOLD citation:

Ratnasingham, Sujeevan, and Paul D N Hebert. “bold: The Barcode of Life Data System (http://www.barcodinglife.org).” Molecular ecology notes vol. 7,3 (2007): 355-364. doi:10.1111/j.1471-8286.2007.01678.x

## Running coidb

The general syntax for running `coidb` is:

```bash
coidb run <arguments>
```

A typical run could look like this:

```bash
coidb run -i data/BOLD_Public.04-Jul-2025.tar.gz -o results -c 4
```

In this example, the file `data/BOLD_Public.04-Jul-2025.tar.gz` was downloaded
from [boldsystems.org/](https://boldsystems.org/) (read more about how to obtain
the input data under [Obtain data](#obtain-data)), output is stored in the
`results` directory and 4 threads are used for running `coidb`.

To see a list of all arguments, run `coidb run -h`. The available arguments are listed below:

```bash
--input-file       -i PATH        Input tar.gz archive dowloaded from BOLD.
--output-dir       -o PATH        Folder to store database files in [default: results]
--account          -A TEXT        SLURM compute account [default: None]
--temp-dir            PATH        Folder for temporary files [default: tmp]
--gbif-backbone                   Use GBIF backbone to infer consensus taxonomy for BOLD BINs
--consensus-threshold INTEGER     Threshold (in %) when calculating consensus taxonomy [default: 80]
--consensus-method    [rank|full] Method to use when calculating consensus [default: rank]
--vsearch-identity    FLOAT       Identity at which to cluster sequences per BIN [default: 1.0]
--ranks               TEXT        Ranks to use for calculating consensus and generating fastas [default: kingdom, phylum, class, order, family, genus, species]
--min-len             INTEGER     Minimum length of sequences to include [default: 500]
--batch-size          INTEGER     Number of BOLD BINs per batch for running vsearch [default: 50000]
```

* The `--input-file` `-i` argument must point to a BOLD tar archive that you
  have downloaded from [boldsystems.org](https://boldsystems.org/) (see [Obtain
  data](#obtain-data) below). 
* The `--output-dir` or `-o` argument is a directory in which the output from
  `coidb` will be stored (see details under [Output](#output) below).
* The `--account` or `-A` argument sets a compute account for running on SLURM
  clusters (see [Cluster execution](#cluster-execution) below).
* The `--temp-dir` argument sets a directory to use for storing temporary
  output. This directory can be deleted once `coidb` finishes.
* The `--gbif-backbone` argument instructs `coidb` to use the GBIF backbone
  taxonomy to infer taxonomic information for BOLD BINs. Note that this option
  is currently not reliable because of outdated GBIF data.
* The `--consensus-threshold` specifies a threshold in percent when calculating
  consensus taxonomies for BOLD BINs. 
* The `--consensus-method` argument specifies how the consensus taxonomy is
  calculated. With `rank` (default), a consensus is calculated at each rank
  separately starting from 'species' and moving up in the hierarchy (genus,
  family etc.). If a consensus above the consensus-threshold is found at any
  rank, the taxonomy at that rank and its parent lineage is used as taxonomy for
  the BOLD BIN. With `full`, a consensus is applied by taking into account the
  parent lineages at each rank, so starting with all labels from
  kingdom->species, then kingdom->genus etc.
* The `--vsearch-identity` argument specifies the identity threshold to use when
  clustering sequences with vsearch. The default is `1.0` meaning sequences are
  clustered at 100% identity.
* The `--ranks` argument specifies what taxonomic ranks to use. This applies
  both to what ranks are included in the final output and what ranks are used to
  calculate the consensus taxonomy.
* The `--min-len` argument sets a minimum length for sequences to include in the
  final output.
* The `--batch-size` argument sets the number of BOLD bins to process with
  vsearch in parallell. This is used to reduce the the size of the workflow
  graph by splitting the input sequences into batches with `batch-size` number
  of BOLD bins per file.

In addition to these command line arguments there are some arguments that define how `coidb` runs on your system and which are similar to how you typically interact with Snakemake workflows:

```bash
--config             FILE     Path to snakemake config file. Overrides existing workflow configuration. [default: None] 
--resource        -r PATH     Additional resources to copy from workflow directory at run time.
--profile         -p TEXT     Name of profile to use for configuring Snakemake. [default: None]
--dry             -n          Do not execute anything, and display what would be done.
--lock            -l          Lock the working directory.
--dag             -d PATH     Save directed acyclic graph to file. Must end in .pdf, .png or .svg [default: None]
--cores           -c INTEGER  Set the number of cores to use. If None will use all cores. [default: None]
--no-conda.                   Do not use conda environments.
--keep-resources              Keep resources after pipeline completes.
--keep-snakemake              Keep .snakemake folder after pipeline completes.
--verbose         -v          Run workflow in verbose mode.
--help-snakemake  -hs         Print the snakemake help and exit.
--help            -h                Show this message and exit.
```

* The `--config` argument lets you pass a configuration file in YAML format, as
  an alternative to specifying arguments directly on the command line.
* The `--profile` argument specifies configuration profile to use for running `coidb`.

### Using a configuration file

You can generate a default configuration file by running:

```bash
coidb config > config.yml
```

This creates a new file `config.yml` with the following default parameters:

```yaml
account: ''
batch_size: 50000
consensus_method: rank
consensus_threshold: 80
gbif_backbone: false
input_file: null
min_len: 500
output_dir: results
ranks:
- kingdom
- phylum
- class
- order
- family
- genus
- species
temp_dir: tmp
vsearch_identity: 1.0
```

You can then edit this file and use it with coidb like so:

```bash
coidb run --config config.yml <additional arguments>
```

### Cluster execution

To run `coidb` on a compute cluster you must set the compute account to use with
the `--account` or `-A` argument. In addition, you should use one of the
pre-defined configuration profiles. To see available profiles, run:

```bash
coidb profile list
```

The profiles can be used as-is by adding `--profile <name of profile>` to the
command line call, or you can output the settings for a profile to a file and
modify it you fit your needs. For example, to output the generic SLURM profile
settings, we recommend that you run:

```bash
mkdir my-slurm-profile
coidb profile show slurm > my-slurm-profile/config.yaml
```

Then you can edit the `my-slurm-profile/config.yaml`. Once you're done you can
use this profile by passing `--profile my-slurm-profile` to the `coidb run`
command.

...

## Output

The primary outputs from a run are placed in the directory set by the `--output-dir` command line argument (default: `results/`). These include:

* `coidb.clustered.fasta.gz`: A fasta file with sequences clustered at whatever
   threshold set in the config file (default is 1.0 which means 100% identity).
   Sequence ids in this file correspond to process_ids, _e.g._ `BPALB370-17`,
   and can be looked up in the [BOLD
   portal](portal.boldsystems.org/result?query=BPALB370-17[ids]). The fasta
   header also includes the corresponding BOLD BIN id of the sequence, _e.g._
   `bin_uri:BOLD:AAF7702`.

* `coidb.info.tsv.gz`: This TSV file contains sequence and taxonomic information
  for all records kept after filtering.

* `coidb.BOLD_BIN.consensus_taxonomy.tsv.gz`: This TSV file contains the
  calculated consensus taxonomy of BOLD BINs. If a consensus could not be
  reached at a certain taxonomic rank, the taxonomic label at that rank is
  prefixed with 'unresolved.' followed by the label of the lowest consensus
  rank.

* `sintax/coidb.sintax.fasta.gz`: This fasta file is compatible with the SINTAX classification tool implemented in [vsearch](https://github.com/torognes/vsearch) and has headers with the format:

```
>BPALB370-17;tax=k:Animalia,p:Arthropoda,c:Insecta,o:Lepidoptera,f:Lycaenidae,g:Thersamonia,s:Thersamonia_X,t:BOLD:AAF7702
```

* `dada2/coidb.dada2.toGenus.fasta.gz`, `dada2/coidb.dada2.toSpecies.fasta.gz` and `dada2/coidb.dada2.addSpecies.fasta.gz`: These fasta files are compatible with the `assignTaxonomy` and `addSpecies` functions from [DADA2](https://benjjneb.github.io/dada2/assign.html).

* `qiime2/coidb.qiime2.info.tsv.gz`: This TSV file can be used with QIIME2 to create a taxonomy artifact for use with the [feature-classifier](https://amplicon-docs.qiime2.org/en/latest/references/plugins/feature-classifier.html#q2-plugin-feature-classifier) plugin. Unzip the file then run `qiime tools import --type 'FeatureData[Taxonomy]' --input-format TSVTaxonomyFormat --input-path coidb.qiime2.info.tsv --output-path taxonomy.qza`. The `coidb.clustered.fasta.gz` file can be used to import sequences with `qiime tools import --type 'FeatureData[Sequence]' --input-path coidb.clustered.fasta --output-path seqs.qza`.

### Log files

Log files from a coidb run are stored under `_logs/` in your `--output-dir` directory (default: `results/`).

### Temporary files

Temporary files are stored under the directory defined by `--temp-dir` (default:
`tmp/`). This entire directory can be removed after a successful run of coidb,
but can also be used to get an idea of what has happened to the raw data along
the way.

For example, the `_extract/` directory contains the raw TSV file extracted from
the input file you used, while the `_processed/` directory contains _e.g._ the
`data.filtered.tsv` file which is the output from the filtering step of coidb.

## How it works

### Filtering
Firstly, the input file is extracted and the TSV file with taxonomic information
and sequence data for each record is identified. This TSV file is then filtered
by:

1. Selecting a useful subset of columns
2. Only keeping records with `COI-5P` in the `marker_code` column.
3. Only keeping records assigned to a proper BOLD BIN (with the exception of prokaryotic sequences which are all kept).
4. Removing records with sequences that are too short (as defined by the `--min-len` argument).
5. Stripping any leading and trailing gap (`-`) characters
6. Removing sequences with remaining gaps
7. Removing sequences with non DNA characters.

### Filling missing data

The filtered TSV file is then processed to fill in missing values for taxonomic
ranks. Ranks with missing data is filled with the lowest assigned taxonomic
label, suffixed with `_X`. For example:

| processid    | kingdom   | phylum          | class | order       | family | genus | species | bin_uri      |
|--------------|-----------|-----------------|-------|-------------|--------|-------|---------|--------------|
| DUTCH124-19  | Animalia  | Platyhelminthes | None  | Polycladida | None   |  None |  None   | BOLD:ACC8697 |
| AACTA1367-20 | Animalia  | Arthropoda      | None  |  None       | None   |  None |  None   | BOLD:AED1280 |

becomes:

| processid    | kingdom   | phylum          | class | order       | family | genus | species | bin_uri      |
|--------------|-----------|-----------------|-------|-------------|--------|-------|---------|--------------|
| DUTCH124-19  | Animalia  | Platyhelminthes | Platyhelminthes_X  | Polycladida | Polycladida_X   |  Polycladida_XX |  Polycladida_XXX   | BOLD:ACC8697 |
| AACTA1367-20 | Animalia  | Arthropoda      | Arthropoda_X  |  Arthropoda_XX       | Arthropoda_XXX   |  Arthropoda_XXXX |  Arthropoda_XXXXX   | BOLD:AED1280 |

### Fix non-unique taxa

Some records in the BOLD data may have the same taxonomic labels at a specific
rank, but with different labels for parent ranks. Take for example the
`Hemineura` genus where records may have these conflicting labels for higher
ranks:

| kingdom | phylum | class | order | family | genus |
|---------|--------|-------|-------|--------|-------|
| Animalia | Arthropoda | Insecta | Psocodea | Elipsocidae | Hemineura |
| Protista | Rhodophyta | Florideophyceae | Ceramiales | Delesseriaceae | Hemineura |

This is dealt with by either removing BOLD BINs that lack taxonomic information
for higher ranks, or by prefixing the non-unique rank with the label of the higher taxonomic rank. In the example above, this would generate:

| kingdom | phylum | class | order | family | genus |
|---------|--------|-------|-------|--------|-------|
| Animalia | Arthropoda | Insecta | Psocodea | Elipsocidae | Elipsocidae_Hemineura |
| Protista | Rhodophyta | Florideophyceae | Ceramiales | Delesseriaceae | Delesseriaceae_Hemineura |

### Consensus taxonomy

After these steps, a consensus taxonomy is calculated for BOLD BINs by taking
into account the taxonomic information for all records in each BIN. For example, a BOLD BIN with records labelled:

| kingdom | phylum | class | order | family | genus | species |
|---------|--------|-------|-------|--------|-------|---------|
|  K |  P |  C |  O |  F |  G |  S |
|  K |  P |  C |  O |  F |  G |  S |
|  K |  P |  C |  O |  F |  G |  S |
|  K |  P |  C |  O |  F |  G |  S |
|  K |  P |  C |  O |  F |  G |  S2 |
|  K |  P |  C |  O |  F |  G |  G_X |
|  K |  P |  C |  O |  F |  G |  G_X |
|  K |  P |  C |  O |  F |  G |  G_X |
|  K |  P |  C |  O |  F |  G |  G_X |

This BIN has 4 records labelled species `S`, 1 record labelled species `S2` and
4 records with missing species labels (these were given the genus label suffixed
with `_X`). Starting from species, 40% of records are labelled species `S`, 10%
are labelled `S2` and 40% have ambiguous labels. Using a consensus threshold of
80% (the default) we see that no consensus can be reached for this BIN at
species level. Moving one step up in the hierarchy however gets us 100% of
records labelled genus `G`. Consequently this BIN will receive a consensus
taxonomy:

| kingdom | phylum | class | order | family | genus | species |
|---------|--------|-------|-------|--------|-------|---------|
| K | P | C | O | F | G | unresolved.G |

Using the command line argument `--consensus-exclude-missing-data` all the
records with `G_X` at rank=species would be ignored in which case there are 80%
records with species `S` and 20% with species `S2`. The consensus taxonomy for
the BIN would then become:

| kingdom | phylum | class | order | family | genus | species |
|---------|--------|-------|-------|--------|-------|---------|
| K | P | C | O | F | G | S|