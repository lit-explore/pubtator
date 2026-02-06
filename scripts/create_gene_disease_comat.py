"""
Creates a gene-disease co-occurrence matrix
"""

import json

import numpy as np
import pandas as pd

snek = snakemake

# load gene and disease pmid mappings
with open(snek.input[0]) as fp:
    gene_pmid_mapping = json.load(fp)

with open(snek.input[1]) as fp:
    disease_pmid_mapping = json.load(fp)

# pre-convert to sets for faster intersection operations
gene_pmid_sets = {k: set(v) for k, v in gene_pmid_mapping.items()}
disease_pmid_sets = {k: set(v) for k, v in disease_pmid_mapping.items()}

# create empty matrix to store gene-disease co-occurrence counts
entrez_ids = list(gene_pmid_sets.keys())
num_genes = len(entrez_ids)

mesh_ids = list(disease_pmid_sets.keys())
num_diseases = len(mesh_ids)

comat = np.zeros((num_genes, num_diseases), dtype=np.uint32)

print(f"Computing {num_genes} x {num_diseases} co-occurrence matrix...")

# iterate over genes & diseases
for i, gene in enumerate(entrez_ids):
    gene_pmids = gene_pmid_sets[gene]

    for j, disease in enumerate(mesh_ids):
        disease_pmids = disease_pmid_sets[disease]
        comat[i, j] = len(gene_pmids & disease_pmids)

# store gene-disease co-occurrence matrix
comat = pd.DataFrame(comat, index=entrez_ids, columns=mesh_ids)
comat.reset_index().rename(columns={"index": "entrez_id"}).to_feather(snek.output[0])
