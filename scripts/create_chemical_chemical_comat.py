"""
Creates a chemical-chemical co-occurrence matrix
"""

import json

import numpy as np
import pandas as pd

snek = snakemake

# load chemical/pmid mapping
with open(snek.input[0]) as fp:
    unfiltered_chemical_pmids = json.load(fp)

# filter any chemicals with less than N citations; helps keep the co-occurrence
# matrix size more manageable while only losing relatively low-information
# low-citation counts..
concept_id_min_freq = snek.config["filtering"]["chem_comat_concept_id_min_freq"]

# pre-convert to sets and filter in one pass
chemical_pmid_sets = {
    mesh_id: set(pmids)
    for mesh_id, pmids in unfiltered_chemical_pmids.items()
    if len(pmids) > concept_id_min_freq
}

# create empty matrix to store chemical-chemical co-occurrence counts
mesh_ids = list(chemical_pmid_sets.keys())
num_chemicals = len(mesh_ids)

comat = np.zeros((num_chemicals, num_chemicals), dtype=np.uint32)

print(f"Computing {num_chemicals} x {num_chemicals} co-occurrence matrix...")

# get upper triangular matrix indices
ind = np.triu_indices(num_chemicals, k=1)

# iterate over pairs of chemicals
for cur_ind in range(len(ind[0])):
    i = ind[0][cur_ind]
    j = ind[1][cur_ind]

    chemical1_pmids = chemical_pmid_sets[mesh_ids[i]]
    chemical2_pmids = chemical_pmid_sets[mesh_ids[j]]

    num_shared = len(chemical1_pmids & chemical2_pmids)
    comat[i, j] = comat[j, i] = num_shared

# store chemical-chemical co-occurrence matrix
comat = pd.DataFrame(comat, index=mesh_ids, columns=mesh_ids)
comat.reset_index().rename(columns={"index": "mesh_id"}).to_feather(snek.output[0])
