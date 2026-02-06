"""
Creates a gene-gene co-occurrence matrix
"""

import json

import numpy as np
import pandas as pd

snek = snakemake

# load gene/pmid mapping
with open(snek.input[0]) as fp:
    gene_pmids = json.load(fp)

# pre-convert to sets for faster intersection operations
gene_pmid_sets = {k: set(v) for k, v in gene_pmids.items()}

# create empty matrix to store gene-gene co-occurrence counts
entrez_ids = list(gene_pmid_sets.keys())
num_genes = len(entrez_ids)

comat = np.zeros((num_genes, num_genes), dtype=np.uint32)

print(f"Computing {num_genes} x {num_genes} co-occurrence matrix...")

# indices of upper triangular matrix
indices = np.triu_indices(num_genes, k=1)

# iterate over pairs of genes and compute overlap
for ind in range(len(indices[0])):
    i = indices[0][ind]
    j = indices[1][ind]

    gene1_pmids = gene_pmid_sets[entrez_ids[i]]
    gene2_pmids = gene_pmid_sets[entrez_ids[j]]

    num_shared = len(gene1_pmids & gene2_pmids)
    comat[i, j] = comat[j, i] = num_shared

# store gene-gene co-occurrence matrix
comat = pd.DataFrame(comat, index=entrez_ids, columns=entrez_ids)
comat.reset_index().rename(columns={"index": "entrez_id"}).to_feather(snek.output[0])
