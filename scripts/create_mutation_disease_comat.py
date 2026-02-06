"""
Creates a mutation-disease co-occurrence matrix
"""

import json

import numpy as np
import pandas as pd

snek = snakemake

# load mutation and disease pmid mappings
with open(snek.input[0]) as fp:
    mutation_pmid_mapping = json.load(fp)

with open(snek.input[1]) as fp:
    disease_pmid_mapping = json.load(fp)

# pre-convert to sets for faster intersection operations
mutation_pmid_sets = {k: set(v) for k, v in mutation_pmid_mapping.items()}
disease_pmid_sets = {k: set(v) for k, v in disease_pmid_mapping.items()}

# create empty matrix to store mutation-disease co-occurrence counts
mutation_ids = list(mutation_pmid_sets.keys())
num_mutations = len(mutation_ids)

mesh_ids = list(disease_pmid_sets.keys())
num_diseases = len(mesh_ids)

comat = np.zeros((num_mutations, num_diseases), dtype=np.uint32)

print(f"Computing {num_mutations} x {num_diseases} co-occurrence matrix...")

# iterate over mutations & diseases
for i, mutation in enumerate(mutation_ids):
    mutation_pmids = mutation_pmid_sets[mutation]

    for j, disease in enumerate(mesh_ids):
        disease_pmids = disease_pmid_sets[disease]
        comat[i, j] = len(mutation_pmids & disease_pmids)

# store mutation-disease co-occurrence matrix
comat = pd.DataFrame(comat, index=mutation_ids, columns=mesh_ids)
comat.reset_index().rename(columns={"index": "mutation"}).to_feather(snek.output[0])
