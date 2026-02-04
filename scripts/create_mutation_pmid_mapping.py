"""
Create a mapping from mutations to pubmed article ids
"""
import json
import numpy as np
import pandas as pd

snek = snakemake

# load mutation data
mutation_dat = pd.read_feather(snek.input[0])

# convert concept id to pyarrow/string type (~7x faster in testing..)
mutation_dat.concept_id = mutation_dat.concept_id.astype('string[pyarrow]')

mutations = list(mutation_dat.concept_id.unique())
num_mutations = len(mutations)

# iterate over mutations and retrieve associated pubmed ids for each
mutation_pmids = {}

for mutation in mutations:
    mask = mutation_dat.concept_id == mutation

    pmids = list(set(mutation_dat[mask].pmid.values))

    if len(pmids) > 0:
        mutation_pmids[mutation] = pmids

def encoder(obj) -> int|list:
    """
    encoder to convert int64 elements to generic ints and sets to lists during
    json serialization
    """
    if isinstance(obj, np.generic):
        return obj.item()
    return obj

# store mutation -> pmid mapping as json
with open(snek.output[0], "wt", encoding="utf-8") as fp:
    fp.write(json.dumps(mutation_pmids, default=encoder))
