"""
Create a mapping from disease mesh ids to pubmed article ids
"""

import pandas as pd
from sparse_utils import save_pmid_mapping

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

# store disease -> pmid mapping as gzip-compressed json
save_pmid_mapping(disease_pmids, snek.output[0])
