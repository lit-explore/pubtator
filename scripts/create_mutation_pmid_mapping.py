"""
Create a mapping from mutations to pubmed article ids
"""

import json

import pandas as pd
from utils import numpy_encoder

snek = snakemake

# load mutation data
mutation_dat = pd.read_feather(snek.input[0])

# build mutation -> pmid mapping using groupby (much faster than iteration)
mutation_pmids = (
    mutation_dat.groupby("concept_id")["pmid"]
    .apply(lambda x: x.unique().tolist())
    .to_dict()
)

print(f"Mapped {len(mutation_pmids)} mutations to PMIDs")

# store mutation -> pmid mapping as json
with open(snek.output[0], "wt", encoding="utf-8") as fp:
    json.dump(mutation_pmids, fp, default=numpy_encoder)
