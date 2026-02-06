"""
Create a mapping from chemical mesh ids to pubmed article ids
"""

import pandas as pd
from sparse_utils import save_pmid_mapping

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

# store chemical -> pmid mapping as gzip-compressed json
save_pmid_mapping(chemical_pmids, snek.output[0])
