"""
Create a mapping from disease mesh ids to pubmed article ids
"""

import json

import pandas as pd
from utils import numpy_encoder

snek = snakemake

# load disease data
disease_dat = pd.read_feather(snek.input[0])

# build disease -> pmid mapping using groupby (much faster than iteration)
disease_pmids = (
    disease_dat.groupby("concept_id")["pmid"]
    .apply(lambda x: x.unique().tolist())
    .to_dict()
)

print(f"Mapped {len(disease_pmids)} diseases to PMIDs")

# store disease -> pmid mapping as json
with open(snek.output[0], "wt", encoding="utf-8") as fp:
    json.dump(disease_pmids, fp, default=numpy_encoder)
