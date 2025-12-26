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
        self.info = read_df(f"{self.output_dir}/coidb.info.tsv.gz")
        self.consensus = read_df(
            f"{self.output_dir}/coidb.BOLD_BIN.consensus_taxonomy.tsv.gz"
        )
        self.clustered_fasta = read_fasta(f"{self.output_dir}/coidb.clustered.fasta.gz")
        self.dada2_addspecies = read_fasta(
            f"{self.output_dir}/dada2/coidb.dada2.addSpecies.fasta.gz"
        )
        self.dada2_toGenus = read_fasta(
            f"{self.output_dir}/dada2/coidb.dada2.toGenus.fasta.gz"
        )
        self.dada2_toSpecies = read_fasta(
            f"{self.output_dir}/dada2/coidb.dada2.toSpecies.fasta.gz"
        )
        self.sintax = read_fasta(f"{self.output_dir}/sintax/coidb.sintax.fasta.gz")


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
            "genus": ["Arhodia", "Arhodia", "Lepidoptera_XX"],
            "species": ["Arhodia lasiocamparia", "Arhodia AH03", "Arhodia AH03"],
            "n": [2, 6, 2],
            "bin_uri": ["test"] * 3,
        }
    )


def test_consensus_taxonomy(taxdata):
    ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    from coidb.scripts import consensus_taxonomy

    assert (
        consensus_taxonomy.calculate_consensus(
            taxdata, ranks=ranks, threshold=80, method="full"
        )
        .select("species")
        .item()
        == "unresolved.Arhodia"
    )

    assert (
        consensus_taxonomy.calculate_consensus(
            taxdata, ranks=ranks, threshold=80, method="rank"
        )
        .select("species")
        .item()
        == "Arhodia AH03"
    )

    assert (
        consensus_taxonomy.calculate_consensus(
            taxdata, ranks=ranks, threshold=80, method="rank"
        )
        .select("genus")
        .item()
        == "Arhodia"
    )

    assert (
        consensus_taxonomy.calculate_consensus(
            taxdata, ranks=ranks, threshold=90, method="rank"
        )
        .select("genus")
        .item()
        == "unresolved.Lepidoptera"
    )
    assert (
        consensus_taxonomy.calculate_consensus(
            taxdata, ranks=ranks, threshold=90, method="full"
        )
        .select("genus")
        .item()
        == "unresolved.Lepidoptera"
    )

    ranks.pop()

    assert (
        consensus_taxonomy.calculate_consensus(
            taxdata, ranks=ranks, threshold=80, method="full"
        )
        .select("genus")
        .item()
        == "Arhodia"
    )

    assert (
        consensus_taxonomy.calculate_consensus(
            taxdata, ranks=ranks, threshold=80, method="rank"
        )
        .select("genus")
        .item()
        == "Arhodia"
    )


def test_returncode(workflow_runs):
    assert all(r.returncode == 0 for r in workflow_runs)


def test_files_exist(workflow_runs):
    for r in workflow_runs:
        output_dir = r.output_dir
        assert all(
            [
                os.path.exists(f"{output_dir}/coidb.info.tsv.gz"),
                os.path.exists(
                    f"{output_dir}/coidb.BOLD_BIN.consensus_taxonomy.tsv.gz"
                ),
                os.path.exists(f"{output_dir}/coidb.clustered.fasta.gz"),
                os.path.exists(f"{output_dir}/dada2/coidb.dada2.addSpecies.fasta.gz"),
                os.path.exists(f"{output_dir}/dada2/coidb.dada2.toGenus.fasta.gz"),
                os.path.exists(f"{output_dir}/dada2/coidb.dada2.toSpecies.fasta.gz"),
                os.path.exists(f"{output_dir}/sintax/coidb.sintax.fasta.gz"),
                os.path.exists(f"{output_dir}/qiime2/coidb.qiime2.info.tsv.gz"),
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
