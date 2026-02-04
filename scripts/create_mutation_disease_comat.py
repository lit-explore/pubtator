"""
Creates a mutation-disease co-occurrence matrix
"""
import json
import pandas as pd
import numpy as np

snek = snakemake

# load mutation and drug pmid mappings
with open(snek.input[0]) as fp:
    mutation_pmid_mapping = json.load(fp)

with open(snek.input[1]) as fp:
    disease_pmid_mapping = json.load(fp)

# create empty matrix to store mutation-disease co-occurrence counts
mutation_ids = mutation_pmid_mapping.keys()
num_mutations = len(mutation_ids)

mesh_ids = disease_pmid_mapping.keys()
num_diseases = len(mesh_ids)

comat = np.zeros((num_mutations, num_diseases), dtype=np.uint32)

# iterate over mutations & diseases
for i, mutation in enumerate(mutation_ids):
    # get pubmed ids associated with mutation
    mutation_pmids = mutation_pmid_mapping[mutation]

    for j, disease in enumerate(mesh_ids):
        # get pubmed ids associated with disease
        disease_pmids = disease_pmid_mapping[disease]

        # compute mutation-disease co-occurrence count
        num_shared = len(set(mutation_pmids).intersection(disease_pmids))

        comat[i, j] = num_shared

# store mutation-disease co-occurrence matrix
comat = pd.DataFrame(comat, index=mutation_ids, columns=mesh_ids)
comat.reset_index().rename(columns={'index': 'mutation'}).to_feather(snek.output[0])
