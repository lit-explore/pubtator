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

# Co-occurrence score types
SCORE_TYPES = ["counts", "dice", "npmi"]

# Species configuration
SPECIES = config.get("species", [{"id": "9606", "name": "human", "gene_symbols_file": "data/symbols.txt", "gene_ensembl_file": "data/ensgenes.txt"}])
SPECIES_NAMES = [s["name"] for s in SPECIES]

# Helper to get species config by name
def get_species_config(species_name):
    for s in SPECIES:
        if s["name"] == species_name:
            return s
    raise ValueError(f"Unknown species: {species_name}")

rule all:
    input:
        # Counts and summaries
        expand(os.path.join(config["output_dir"], "counts", "{annot}_concept_ids.feather"), annot=ANNOT_TYPES),
        expand(os.path.join(config["output_dir"], "counts", "{annot}_mentions.feather"), annot=ANNOT_TYPES),
        os.path.join(config["output_dir"], "counts", "unique_genes.feather"),
        # ID mappings (per species)
        expand(os.path.join(config["output_dir"], "ids", "{species}", "symbol2entrez.feather"), species=SPECIES_NAMES),
        expand(os.path.join(config["output_dir"], "ids", "{species}", "ensgene2entrez.feather"), species=SPECIES_NAMES),
        # PMID mappings (gzip-compressed JSON)
        expand(os.path.join(config["output_dir"], "json", "{species}_entrez_pmids.json.gz"), species=SPECIES_NAMES),
        os.path.join(config["output_dir"], "json", "disease_pmids.json.gz"),
        os.path.join(config["output_dir"], "json", "chemical_pmids.json.gz"),
        os.path.join(config["output_dir"], "json", "mutation_pmids.json.gz"),
        # Co-occurrence matrices (sparse .npz format with multiple score types)
        expand(os.path.join(config["output_dir"], "co-occurrence", "{species}_gene_gene_{score}.npz"), species=SPECIES_NAMES, score=SCORE_TYPES),
        expand(os.path.join(config["output_dir"], "co-occurrence", "{species}_gene_gene_labels.json"), species=SPECIES_NAMES),
        expand(os.path.join(config["output_dir"], "co-occurrence", "disease_disease_{score}.npz"), score=SCORE_TYPES),
        os.path.join(config["output_dir"], "co-occurrence", "disease_disease_labels.json"),
        expand(os.path.join(config["output_dir"], "co-occurrence", "chemical_chemical_{score}.npz"), score=SCORE_TYPES),
        os.path.join(config["output_dir"], "co-occurrence", "chemical_chemical_labels.json"),
        expand(os.path.join(config["output_dir"], "co-occurrence", "{species}_gene_chemical_{score}.npz"), species=SPECIES_NAMES, score=SCORE_TYPES),
        expand(os.path.join(config["output_dir"], "co-occurrence", "{species}_gene_chemical_labels.json"), species=SPECIES_NAMES),
        expand(os.path.join(config["output_dir"], "co-occurrence", "{species}_gene_disease_{score}.npz"), species=SPECIES_NAMES, score=SCORE_TYPES),
        expand(os.path.join(config["output_dir"], "co-occurrence", "{species}_gene_disease_labels.json"), species=SPECIES_NAMES),
        expand(os.path.join(config["output_dir"], "co-occurrence", "mutation_disease_{score}.npz"), score=SCORE_TYPES),
        os.path.join(config["output_dir"], "co-occurrence", "mutation_disease_labels.json"),

rule create_gene_gene_comat:
    input:
        os.path.join(config["output_dir"], "json", "{species}_entrez_pmids.json.gz"),
    output:
        expand(os.path.join(config["output_dir"], "co-occurrence", "{{species}}_gene_gene_{score}.npz"), score=SCORE_TYPES),
        os.path.join(config["output_dir"], "co-occurrence", "{species}_gene_gene_labels.json"),
    params:
        prefix=lambda wildcards: f"{wildcards.species}_gene_gene"
    script:
        "scripts/create_gene_gene_comat.py"

rule create_chemical_chemical_comat:
    input:
        os.path.join(config["output_dir"], "json", "chemical_pmids.json.gz")
    output:
        expand(os.path.join(config["output_dir"], "co-occurrence", "chemical_chemical_{score}.npz"), score=SCORE_TYPES),
        os.path.join(config["output_dir"], "co-occurrence", "chemical_chemical_labels.json"),
    script:
        "scripts/create_chemical_chemical_comat.py"

rule create_disease_disease_comat:
    input:
        os.path.join(config["output_dir"], "json", "disease_pmids.json.gz")
    output:
        expand(os.path.join(config["output_dir"], "co-occurrence", "disease_disease_{score}.npz"), score=SCORE_TYPES),
        os.path.join(config["output_dir"], "co-occurrence", "disease_disease_labels.json"),
    script:
        "scripts/create_disease_disease_comat.py"

rule create_gene_chemical_comat:
    input:
        os.path.join(config["output_dir"], "json", "{species}_entrez_pmids.json.gz"),
        os.path.join(config["output_dir"], "json", "chemical_pmids.json.gz")
    output:
        expand(os.path.join(config["output_dir"], "co-occurrence", "{{species}}_gene_chemical_{score}.npz"), score=SCORE_TYPES),
        os.path.join(config["output_dir"], "co-occurrence", "{species}_gene_chemical_labels.json"),
    params:
        prefix=lambda wildcards: f"{wildcards.species}_gene_chemical"
    script:
        "scripts/create_gene_chemical_comat.py"

rule create_mutation_disease_comat:
    input:
        os.path.join(config["output_dir"], "json", "mutation_pmids.json.gz"),
        os.path.join(config["output_dir"], "json", "disease_pmids.json.gz")
    output:
        expand(os.path.join(config["output_dir"], "co-occurrence", "mutation_disease_{score}.npz"), score=SCORE_TYPES),
        os.path.join(config["output_dir"], "co-occurrence", "mutation_disease_labels.json"),
    script:
        "scripts/create_mutation_disease_comat.py"

rule create_gene_disease_comat:
    input:
        os.path.join(config["output_dir"], "json", "{species}_entrez_pmids.json.gz"),
        os.path.join(config["output_dir"], "json", "disease_pmids.json.gz")
    output:
        expand(os.path.join(config["output_dir"], "co-occurrence", "{{species}}_gene_disease_{score}.npz"), score=SCORE_TYPES),
        os.path.join(config["output_dir"], "co-occurrence", "{species}_gene_disease_labels.json"),
    params:
        prefix=lambda wildcards: f"{wildcards.species}_gene_disease"
    script:
        "scripts/create_gene_disease_comat.py"

rule create_chemical_pmid_mapping:
    input:
        os.path.join(config["output_dir"], "filtered", "chemical.feather")
    output:
        os.path.join(config["output_dir"], "json", "chemical_pmids.json.gz")
    script:
        "scripts/create_chemical_pmid_mapping.py"

rule create_disease_pmid_mapping:
    input:
        os.path.join(config["output_dir"], "filtered", "disease.feather")
    output:
        os.path.join(config["output_dir"], "json", "disease_pmids.json.gz")
    script:
        "scripts/create_disease_pmid_mapping.py"

rule create_mutation_pmid_mapping:
    input:
        os.path.join(config["output_dir"], "filtered", "mutation.feather")
    output:
        os.path.join(config["output_dir"], "json", "mutation_pmids.json.gz")
    script:
        "scripts/create_mutation_pmid_mapping.py"

rule create_gene_pmid_mapping:
    input:
        os.path.join(config["output_dir"], "filtered", "species.feather"),
        os.path.join(config["output_dir"], "filtered", "gene.feather"),
        os.path.join(config["output_dir"], "ids", "{species}", "symbol2entrez.feather"),
        os.path.join(config["output_dir"], "ids", "{species}", "ensgene2entrez.feather"),
    output:
        os.path.join(config["output_dir"], "json", "{species}_entrez_pmids.json.gz")
    params:
        species=lambda wildcards: get_species_config(wildcards.species)["id"]
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
        lambda wildcards: get_species_config(wildcards.species)["gene_symbols_file"]
    output:
        os.path.join(config["output_dir"], "ids", "{species}", "symbol2entrez.feather")
    params:
        species_id=lambda wildcards: int(get_species_config(wildcards.species)["id"]),
        scope="symbol",
        colname="symbol"
    script:
        "scripts/create_gene_mapping.py"

rule create_ensgene2entrez_mapping:
    input:
        lambda wildcards: get_species_config(wildcards.species)["gene_ensembl_file"]
    output:
        os.path.join(config["output_dir"], "ids", "{species}", "ensgene2entrez.feather")
    params:
        species_id=lambda wildcards: int(get_species_config(wildcards.species)["id"]),
        scope="ensembl.gene",
        colname="ensgene"
    script:
        "scripts/create_gene_mapping.py"

rule download_data:
    output:
        os.path.join(config["output_dir"], "raw", "{annot}.gz")
    retries: 3
    shell:
        """
        curl --retry 3 --retry-delay 5 --fail --output {output} \
            https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/{wildcards.annot}2pubtator3.gz

        # Verify file is non-empty
        if [ ! -s {output} ]; then
            echo "Error: Downloaded file is empty"
            rm -f {output}
            exit 1
        fi
        """

# vi:ft=snakemake
