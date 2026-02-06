"""
Creates a gene-chemical co-occurrence matrix
"""

import json

import numpy as np
import pandas as pd

snek = snakemake

# load gene and chemical pmid mappings
with open(snek.input[0]) as fp:
    gene_pmid_mapping = json.load(fp)

with open(snek.input[1]) as fp:
    chem_pmid_mapping = json.load(fp)

# pre-convert to sets for faster intersection operations
gene_pmid_sets = {k: set(v) for k, v in gene_pmid_mapping.items()}
chem_pmid_sets = {k: set(v) for k, v in chem_pmid_mapping.items()}

# create empty matrix to store gene-chemical co-occurrence counts
entrez_ids = list(gene_pmid_sets.keys())
num_genes = len(entrez_ids)

mesh_ids = list(chem_pmid_sets.keys())
num_chemicals = len(mesh_ids)

comat = np.zeros((num_genes, num_chemicals), dtype=np.uint32)

print(f"Computing {num_genes} x {num_chemicals} co-occurrence matrix...")

# iterate over genes & chemicals
for i, gene in enumerate(entrez_ids):
    gene_pmids = gene_pmid_sets[gene]

    if i % 100 == 0:
        print(f"gene {i}/{num_genes}...")

    for j, chemical in enumerate(mesh_ids):
        chemical_pmids = chem_pmid_sets[chemical]
        comat[i, j] = len(gene_pmids & chemical_pmids)

# store gene-chemical co-occurrence matrix
comat = pd.DataFrame(comat, index=entrez_ids, columns=mesh_ids)
comat.reset_index().rename(columns={"index": "entrez_id"}).to_feather(snek.output[0])
