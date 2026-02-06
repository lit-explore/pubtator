"""
Create a mapping from chemical mesh ids to pubmed article ids
"""

import json

import pandas as pd
from utils import numpy_encoder

snek = snakemake

# load chemical data
chemical_dat = pd.read_feather(snek.input[0])

# build chemical -> pmid mapping using groupby (much faster than iteration)
chemical_pmids = (
    chemical_dat.groupby("concept_id")["pmid"]
    .apply(lambda x: x.unique().tolist())
    .to_dict()
)

print(f"Mapped {len(chemical_pmids)} chemicals to PMIDs")

# store chemical -> pmid mapping as json
with open(snek.output[0], "wt", encoding="utf-8") as fp:
    json.dump(chemical_pmids, fp, default=numpy_encoder)
