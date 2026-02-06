# Rules in this Snakefile are meant for matching species names in BOLD to
# accepted species names in GBIF/COL with the aim to clean up erroneous entries
# and synonyms prior to calculating consensus taxonomies for BOLD bins.


localrules:
    download_backbone,


rule download_backbone:
    """
    This rule downloads the latest GBIF backbone release.
    """
    message:
        "Downloading GBIF taxonomy backbone"
    output:
        tsv=temp(expand("{tempdir}/_backbone/Taxon.tsv", tempdir=config["temp_dir"])),
        timestamp=temp(
            expand(
                "{output_dir}/backbone_timestamp.txt", output_dir=config["output_dir"]
            )
        ),
    log:
        os.path.join(config["output_dir"], "_logs/download_backbone.log"),
    params:
        output_dir=lambda wildcards, output: os.path.dirname(output[0]),
        url="https://hosted-datasets.gbif.org/datasets/backbone/current/backbone.zip",
    shell:
        """
        curl -o {params.output_dir}/backbone.zip {params.url} > {log} 2>&1
        echo "GBIF backbone: " > {output.timestamp}
        timestamp=$(unzip -p {params.output_dir}/backbone.zip eml.xml | grep dateStamp | cut -f2 -d '>' | cut -f1 -d '<')
        echo "  - date: $timestamp" >> {output.timestamp}
        echo "  - url: {params.url}" >> {output.timestamp}
        unzip -p {params.output_dir}/backbone.zip Taxon.tsv | cut -f1,3,6,8,12,18- > {output.tsv}
        rm {params.output_dir}/backbone.zip
        """


rule parse_backbone:
    """
    This rule parses the GBIF backbone for downstream use.
    """
    message:
        "Extracting backbone taxonomy for BOLD bins"
    output:
        temp(
            expand(
                "{tempdir}/_backbone/bold_bin.taxonomy.tsv", tempdir=config["temp_dir"]
            )
        ),
    input:
        rules.download_backbone.output[0],
    log:
        os.path.join(config["output_dir"], "_logs/parse_backbone.log"),
    shell:
        """
        parse-backbone {input} {output} 2>{log}
        """


rule match_names:
    """
    Uses the pygbif module to match species names in filtered data to the GBIF
    backbone.
    """
    message:
        "Matching names to GBIF backbone"
    input:
        rules.filter.output.tsv,
    output:
        tsv=temp(
            os.path.join(
                config["temp_dir"], "_processed/data.filtered.names_matched.tsv"
            )
        ),
    log:
        os.path.join(config["output_dir"], "_logs/match_names.log"),
    threads: 4
    shell:
        """
        match-names -i {input} -o {output} -p {threads} > {log} 2>&1
        """


rule consolidate_names:
    """
    Consolidates names from BOLD and GBIF backbone.

    1. If a given BOLD record does not have a GBIF taxonomy match, the BOLD
    assignment is kept.
    
    2. If a higher taxon name is missing in GBIF but not in BOLD, import it from
    BOLD if the next higher rank that is not empty in either GBIF or BOLD has
    the same content. Otherwise leave the higher taxon name empty.
    """
    message:
        "Consolidating names from BOLD and GBIF backbone"
    input:
        filtered=rules.filter.output.tsv,
        matched=rules.match_names.output.tsv,
        backbone=rules.parse_backbone.output[0],
