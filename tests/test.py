import subprocess
import pandas as pd
from yaml import safe_load
import pytest
import gzip as gz
import polars as pl
import os


def read_fasta(f):
    seqid = None
    desc = None
    seq = []
    seqs = {}
    if f.endswith(".gz"):
        openf = gz.open
    else:
        openf = open
    with openf(f, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if seq:
                    seqs[seqid] = "".join(seq)
                    seq.clear()
                seqid = line.split()[0][1:]
            else:
                seq.append(line.rstrip())
        if seq:
            seqs[seqid] = "".join(seq)
    return seqs


def run_workflow(config):
    return subprocess.run(["coidb", "run", "--config", config, "--notemp"])


def read_df(f):
    return pd.read_csv(f, sep="\t", index_col=0, header=0)


class Workflow:
    def __init__(self, config):
        self.config = config
        with open(config, "r") as fhin:
            d = safe_load(fhin)
        self.output_dir = d["output_dir"]
        self.returncode = run_workflow(self.config).returncode
        self.info = read_df(f"{self.output_dir}/coidb/coidb.info.tsv.gz")
        self.consensus = read_df(
            f"{self.output_dir}/consensus_taxonomy/coidb.exclNA.tsv.gz"
        )
        self.clustered_fasta = read_fasta(
            f"{self.output_dir}/coidb/coidb.clustered.fasta.gz"
        )
        self.dada2_addspecies = read_fasta(
            f"{self.output_dir}/dada2/coidb.dada2.addSpecies.exclNA.fasta.gz"
        )
        self.dada2_toGenus = read_fasta(
            f"{self.output_dir}/dada2/coidb.dada2.toGenus.exclNA.fasta.gz"
        )
        self.dada2_toSpecies = read_fasta(
            f"{self.output_dir}/dada2/coidb.dada2.toSpecies.exclNA.fasta.gz"
        )
        self.sintax = read_fasta(
            f"{self.output_dir}/sintax/coidb.sintax.exclNA.fasta.gz"
        )


@pytest.fixture
def workflow_runs():
    return [
        Workflow("tests/config1.yml"),
        Workflow("tests/config2.yml"),
        Workflow("tests/config3.yml"),
    ]


@pytest.fixture
def taxdata():
    return pl.DataFrame(
        {
            "kingdom": ["Animalia"] * 3,
            "phylum": ["Arthropoda"] * 3,
            "class": ["Insecta"] * 3,
            "order": ["Lepidoptera"] * 3,
            "family": ["Geometridae", "Geometridae", "Lepidoptera_X"],
            "genus": ["Arhodia-X", "Arhodia", "Lepidoptera_XX"],
            "species": ["Arhodia AH03", "Arhodia AH03", "Lepidoptera_XXX"],
            "n": [2, 6, 2],
            "bin_uri": ["test"] * 3,
        }
    )


@pytest.fixture
def matched_data():
    return pl.DataFrame(
        {
            "name": [
                "Arhodia lasiocamparia",
                "Arhodia AH03",
                "Homo neanderthalensis",
                "Homo sapiens",
            ],
            "kingdom": ["Animalia"] * 4,
            "phylum": ["Mollusca", "Arthropoda", "Chordata", "Chordata"],
            "class": ["Gastropoda", "Insecta", "Mammalia", "Mammalia"],
            "order": ["Nudibranchia", "Lepidoptera", "Primates", "Primates"],
            "family": ["Janolidae", "Geometridae", "Panidae", "Hominidae"],
            "genus": ["Arhodia", "Arhodia", "Palaeoanthropus", "Homo"],
            "species": [
                "Arhodia lasiocamparia",
                "Arhodia AH03",
                "Palaeoanthropus neanderthalensis",
                "Homo sapiens",
            ],
        }
    )


def test_consensus_taxonomy(taxdata):
    ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    from coidb.scripts import consensus_taxonomy

    # with 'full' method all ranks up to the current are used so for taxdata the
    # full method should resolve neither species nor genus
    assert (
        consensus_taxonomy.calculate_consensus(
            taxdata,
            ranks=ranks,
            threshold=80,
            method="full",
            exclude_missing_data=False,
        )
        .select("species")
        .item()
        == "unresolved.Geometridae"
    )
    assert (
        consensus_taxonomy.calculate_consensus(
            taxdata,
            ranks=ranks,
            threshold=80,
            method="full",
            exclude_missing_data=False,
        )
        .select("genus")
        .item()
        == "unresolved.Geometridae"
    )
    assert (
        consensus_taxonomy.calculate_consensus(
            taxdata,
            ranks=ranks,
            threshold=80,
            method="full",
            exclude_missing_data=False,
        )
        .select("family")
        .item()
        == "Geometridae"
    )
    # raising the threshold to 90 should resolve family with method='full' only
    # if missing data is ignored
    assert (
        consensus_taxonomy.calculate_consensus(
            taxdata,
            ranks=ranks,
            threshold=90,
            method="full",
            exclude_missing_data=True,
        )
        .select("family")
        .item()
        == "Geometridae"
    )
    assert (
        consensus_taxonomy.calculate_consensus(
            taxdata,
            ranks=ranks,
            threshold=90,
            method="full",
            exclude_missing_data=False,
        )
        .select("family")
        .item()
        == "unresolved.Lepidoptera"
    )
    # with 'rank' method, there should be a consensus at species because 8/10
    # records have 'Arhodia AH03' already at species
    assert (
        consensus_taxonomy.calculate_consensus(
            taxdata,
            ranks=ranks,
            threshold=80,
            method="rank",
            exclude_missing_data=False,
        )
        .select("species")
        .item()
        == "Arhodia AH03"
    )
    # raising the threshold to 90 should resolve species only if missing data is
    # excluded, otherwise only order should be resolved
    assert (
        consensus_taxonomy.calculate_consensus(
            taxdata,
            ranks=ranks,
            threshold=90,
            method="rank",
            exclude_missing_data=True,
        )
        .select("species")
        .item()
        == "Arhodia AH03"
    )
    assert (
        consensus_taxonomy.calculate_consensus(
            taxdata,
            ranks=ranks,
            threshold=90,
            method="rank",
            exclude_missing_data=False,
        )
        .select("species")
        .item()
        == "unresolved.Lepidoptera"
    )
    assert (
        consensus_taxonomy.calculate_consensus(
            taxdata,
            ranks=ranks,
            threshold=90,
            method="rank",
            exclude_missing_data=False,
        )
        .select("genus")
        .item()
        == "unresolved.Lepidoptera"
    )
    assert (
        consensus_taxonomy.calculate_consensus(
            taxdata,
            ranks=ranks,
            threshold=90,
            method="rank",
            exclude_missing_data=False,
        )
        .select("family")
        .item()
        == "unresolved.Lepidoptera"
    )
    assert (
        consensus_taxonomy.calculate_consensus(
            taxdata,
            ranks=ranks,
            threshold=90,
            method="rank",
            exclude_missing_data=False,
        )
        .select("order")
        .item()
        == "Lepidoptera"
    )
    ranks.pop()

    assert (
        "species"
        not in consensus_taxonomy.calculate_consensus(
            taxdata, ranks=ranks, threshold=80, method="rank"
        ).columns
    )


def test_returncode(workflow_runs):
    assert all(r.returncode == 0 for r in workflow_runs)


def test_files_exist(workflow_runs):
    for r in workflow_runs:
        output_dir = r.output_dir
        assert all(
            [
                os.path.exists(f"{output_dir}/coidb/coidb.info.tsv.gz"),
                os.path.exists(f"{output_dir}/consensus_taxonomy/coidb.exclNA.tsv.gz"),
                os.path.exists(f"{output_dir}/consensus_taxonomy/coidb.inclNA.tsv.gz"),
                os.path.exists(f"{output_dir}/coidb/coidb.clustered.fasta.gz"),
                os.path.exists(
                    f"{output_dir}/dada2/coidb.dada2.addSpecies.exclNA.fasta.gz"
                ),
                os.path.exists(
                    f"{output_dir}/dada2/coidb.dada2.addSpecies.inclNA.fasta.gz"
                ),
                os.path.exists(
                    f"{output_dir}/dada2/coidb.dada2.toGenus.exclNA.fasta.gz"
                ),
                os.path.exists(
                    f"{output_dir}/dada2/coidb.dada2.toGenus.inclNA.fasta.gz"
                ),
                os.path.exists(
                    f"{output_dir}/dada2/coidb.dada2.toSpecies.exclNA.fasta.gz"
                ),
                os.path.exists(
                    f"{output_dir}/dada2/coidb.dada2.toSpecies.inclNA.fasta.gz"
                ),
                os.path.exists(f"{output_dir}/sintax/coidb.sintax.exclNA.fasta.gz"),
                os.path.exists(f"{output_dir}/sintax/coidb.sintax.inclNA.fasta.gz"),
                os.path.exists(f"{output_dir}/qiime2/coidb.qiime2.info.exclNA.tsv.gz"),
                os.path.exists(f"{output_dir}/qiime2/coidb.qiime2.info.inclNA.tsv.gz"),
            ]
        )


def test_short(workflow_runs):
    info_dfs = [r.info for r in workflow_runs]
    assert all(df.loc[df.index.str.endswith("-short")].shape[0] == 0 for df in info_dfs)


def test_nonDNA(workflow_runs):
    info_dfs = [r.info for r in workflow_runs]
    assert all(
        df.loc[df.index.str.endswith("-nonDNA")].shape[0] == 0 for df in info_dfs
    )


def test_noBIN(workflow_runs):
    info_dfs = [r.info for r in workflow_runs]
    assert all(df.loc[df.index.str.endswith("-noBIN")].shape[0] == 0 for df in info_dfs)


def test_deletion(workflow_runs):
    info_dfs = [r.info for r in workflow_runs]
    assert all(
        df.loc[df.index.str.endswith("-deletion")].shape[0] == 0 for df in info_dfs
    )


def test_wrongmarker(workflow_runs):
    info_dfs = [r.info for r in workflow_runs]
    assert all(
        df.loc[df.index.str.endswith("-wrongmarker")].shape[0] == 0 for df in info_dfs
    )


def test_clustering(workflow_runs):
    r1, r2, r3 = workflow_runs
    assert len([x for x in r1.clustered_fasta.keys() if x.startswith("seq1")]) == 2
    assert len([x for x in r2.clustered_fasta.keys() if x.startswith("seq1")]) == 1
    assert len([x for x in r3.clustered_fasta.keys() if x.startswith("seq1")]) == 1

    assert len([x for x in r1.clustered_fasta.keys() if x.startswith("seq2")]) == 4
    assert len([x for x in r2.clustered_fasta.keys() if x.startswith("seq2")]) == 3
    assert len([x for x in r3.clustered_fasta.keys() if x.startswith("seq2")]) == 2


def test_consolidate_taxonomy(matched_data, workflow_runs):
    r = workflow_runs[0]
    output_dir = r.output_dir
    matched = f"{output_dir}/matched.tsv"
    consolidated = f"{output_dir}/consolidated.tsv"
    matched_data.write_csv(matched, separator="\t")
    res = subprocess.run(
        [
            "consolidate-names",
            "-i",
            f"{output_dir}/coidb/coidb.info.tsv",
            "-m",
            matched,
            "-o",
            consolidated,
        ]
    )
    assert res.returncode == 0
    cons_df = pl.read_csv(consolidated, separator="\t")
    df = pl.read_csv(f"{output_dir}/coidb/coidb.info.tsv", separator="\t")
    joined = df.join(cons_df, on="processid", suffix="_cons")
    primates = (
        joined.filter(pl.col("species") == "Homo neanderthalensis")
        .select("genus_cons", "species_cons")
        .unique()
    )
    assert primates.item(0, 1) == "Palaeoanthropus neanderthalensis"
    assert primates.item(0, 0) == "Palaeoanthropus"
    arhodia1 = (
        joined.filter(pl.col("species") == "Arhodia AH03")
        .select("family_cons", "genus_cons", "species_cons")
        .unique()
    )
    # assert that the consolidated taxonomy for Arhodia AH03 is unique at family, genus and species
    assert arhodia1.height == 1
    assert arhodia1.item(0, 2) == "Arhodia AH03"
    assert arhodia1.item(0, 1) == "Arhodia"
    assert arhodia1.item(0, 0) == "Geometridae"
    arhodia2 = (
        joined.filter(pl.col("species") == "Arhodia lasiocamparia")
        .select(
            "phylum_cons",
            "class_cons",
            "order_cons",
            "family_cons",
            "genus_cons",
            "species_cons",
        )
        .unique()
    )
    # assert that the erroneous match for Arhodia lasiocamparia has not been assigned
    assert arhodia2.item(0, 0) == "Arthropoda"
    assert arhodia2.item(0, 1) == "Insecta"
    assert arhodia2.item(0, 2) == "Lepidoptera"
    assert arhodia2.item(0, 3) == "Geometridae"
    assert arhodia2.item(0, 4) == "Arhodia"
