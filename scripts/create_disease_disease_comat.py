"""
Creates a disease-disease co-occurrence matrix
"""

import json

import numpy as np
import pandas as pd

snek = snakemake

# load disease/pmid mapping
with open(snek.input[0]) as fp:
    disease_pmids = json.load(fp)

# pre-convert to sets for faster intersection operations
disease_pmid_sets = {k: set(v) for k, v in disease_pmids.items()}

# create empty matrix to store disease-disease co-occurrence counts
mesh_ids = list(disease_pmid_sets.keys())
num_diseases = len(mesh_ids)

comat = np.zeros((num_diseases, num_diseases), dtype=np.uint32)

print(f"Computing {num_diseases} x {num_diseases} co-occurrence matrix...")

# get upper triangular matrix indices
ind = np.triu_indices(num_diseases, k=1)

# iterate over pairs of diseases
for cur_ind in range(len(ind[0])):
    i = ind[0][cur_ind]
    j = ind[1][cur_ind]

    disease1_pmids = disease_pmid_sets[mesh_ids[i]]
    disease2_pmids = disease_pmid_sets[mesh_ids[j]]

    num_shared = len(disease1_pmids & disease2_pmids)
    comat[i, j] = comat[j, i] = num_shared

# store disease-disease co-occurrence matrix
comat = pd.DataFrame(comat, index=mesh_ids, columns=mesh_ids)
comat.reset_index().rename(columns={"index": "mesh_id"}).to_feather(snek.output[0])
