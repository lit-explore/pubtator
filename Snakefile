"""
lit-explore: PubTator Central Data Preparation

Processes data from PubTator Central into a form that can be easily used in various lit-explore
efforts.

For more information about the source data, see:

https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/README.txt
"""
import os

configfile: "config/config.yml"

# PubTator Central annotation types
ANNOT_TYPES = ["cellline", "chemical", "disease", "gene", "mutation", "species"]

rule all:
    input:
        expand(os.path.join(config["output_dir"], "counts", "{annot}_concept_ids.feather"), annot=ANNOT_TYPES),
        expand(os.path.join(config["output_dir"], "counts", "{annot}_mentions.feather"), annot=ANNOT_TYPES),
        os.path.join(config["output_dir"], "counts", "unique_genes.feather"),
        os.path.join(config["output_dir"], "ids", "human", "symbol2entrez.feather"),
        os.path.join(config["output_dir"], "ids", "human", "ensgene2entrez.feather"),
        os.path.join(config["output_dir"], "json", "human_entrez_pmids.json"),
        os.path.join(config["output_dir"], "json", "disease_pmids.json"),
        os.path.join(config["output_dir"], "json", "chemical_pmids.json"),
        os.path.join(config["output_dir"], "json", "mutation_pmids.json"),
        os.path.join(config["output_dir"], "co-occurrence", "human_gene_gene_comat.feather"),
        os.path.join(config["output_dir"], "co-occurrence", "human_disease_disease_comat.feather"),
        os.path.join(config["output_dir"], "co-occurrence", "human_chemical_chemical_comat.feather"),
        os.path.join(config["output_dir"], "co-occurrence", "human_gene_chemical_comat.feather"),
        os.path.join(config["output_dir"], "co-occurrence", "human_gene_disease_comat.feather"),
        os.path.join(config["output_dir"], "co-occurrence", "human_mutation_disease_comat.feather")

rule create_gene_gene_comat:
    input:
        os.path.join(config["output_dir"], "json", "human_entrez_pmids.json"),
    output:
        os.path.join(config["output_dir"], "co-occurrence", "human_gene_gene_comat.feather")
    script:
        "scripts/create_gene_gene_comat.py"

rule create_chemical_chemical_comat:
    input:
        os.path.join(config["output_dir"], "json", "chemical_pmids.json")
    output:
        os.path.join(config["output_dir"], "co-occurrence", "human_chemical_chemical_comat.feather")
    script:
        "scripts/create_chemical_chemical_comat.py"

rule create_disease_disease_comat:
    input:
        os.path.join(config["output_dir"], "json", "disease_pmids.json")
    output:
        os.path.join(config["output_dir"], "co-occurrence", "human_disease_disease_comat.feather")
    script:
        "scripts/create_disease_disease_comat.py"

rule create_gene_chemical_comat:
    input:
        os.path.join(config["output_dir"], "json", "human_entrez_pmids.json"),
        os.path.join(config["output_dir"], "json", "chemical_pmids.json")
    output:
        os.path.join(config["output_dir"], "co-occurrence", "human_gene_chemical_comat.feather")
    script:
        "scripts/create_gene_chemical_comat.py"

rule create_mutation_disease_comat:
    input:
        os.path.join(config["output_dir"], "json", "mutation_pmids.json"),
        os.path.join(config["output_dir"], "json", "disease_pmids.json")
    output:
        os.path.join(config["output_dir"], "co-occurrence", "human_mutation_disease_comat.feather")
    script:
        "scripts/create_mutation_disease_comat.py"

rule create_gene_disease_comat:
    input:
        os.path.join(config["output_dir"], "json", "human_entrez_pmids.json"),
        os.path.join(config["output_dir"], "json", "disease_pmids.json")
    output:
        os.path.join(config["output_dir"], "co-occurrence", "human_gene_disease_comat.feather")
    script:
        "scripts/create_gene_disease_comat.py"

rule create_chemical_pmid_mapping:
    input:
        os.path.join(config["output_dir"], "filtered", "chemical.feather")
    output:
        os.path.join(config["output_dir"], "json", "chemical_pmids.json")
    script:
        "scripts/create_chemical_pmid_mapping.py"

rule create_disease_pmid_mapping:
    input:
        os.path.join(config["output_dir"], "filtered", "disease.feather")
    output:
        os.path.join(config["output_dir"], "json", "disease_pmids.json")
    script:
        "scripts/create_disease_pmid_mapping.py"

rule create_mutation_pmid_mapping:
    input:
        os.path.join(config["output_dir"], "filtered", "mutation.feather")
    output:
        os.path.join(config["output_dir"], "json", "mutation_pmids.json")
    script:
        "scripts/create_mutation_pmid_mapping.py"

rule create_human_gene_pmid_mapping:
    input:
        os.path.join(config["output_dir"], "filtered", "species.feather"),
        os.path.join(config["output_dir"], "filtered", "gene.feather"),
        os.path.join(config["output_dir"], "ids", "human", "symbol2entrez.feather"),
        os.path.join(config["output_dir"], "ids", "human", "ensgene2entrez.feather"),
    output:
        os.path.join(config["output_dir"], "json", "human_entrez_pmids.json")
    params:
        species="9606"
    script:
        "scripts/create_gene_pmid_mapping.py"

rule count_genes:
    input:
        os.path.join(config["output_dir"], "filtered", "gene.feather")
    output:
        os.path.join(config["output_dir"], "counts", "unique_genes.feather")
    script:
        "scripts/count_genes.py"

rule summarize_counts:
    input:
        os.path.join(config["output_dir"], "filtered", "{annot}.feather")
    output:
        os.path.join(config["output_dir"], "counts", "{annot}_concept_ids.feather"),
        os.path.join(config["output_dir"], "counts", "{annot}_mentions.feather")
    script:
        "scripts/summarize_counts.py"

rule filter_data:
    input:
        os.path.join(config["output_dir"], "raw", "{annot}.gz")
    output:
        os.path.join(config["output_dir"], "filtered", "{annot}.feather")
    script:
        "scripts/filter_data.py"

rule create_symbol2entrez_mapping:
    input:
        "data/symbols.txt"
    output:
        os.path.join(config["output_dir"], "ids", "human", "symbol2entrez.feather")
    script:
        "scripts/create_human_gene_mapping.py"

rule create_ensgene2entrez_mapping:
    input:
        "data/ensgenes.txt"
    output:
        os.path.join(config["output_dir"], "ids", "human", "ensgene2entrez.feather")
    script:
        "scripts/create_human_gene_mapping.py"

rule download_data:
    output:
        os.path.join(config["output_dir"], "raw", "{annot}.gz")
    shell:
        """
        curl --output {output} https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/{wildcards.annot}2pubtator3.gz
        """

# vi:ft=snakemake
