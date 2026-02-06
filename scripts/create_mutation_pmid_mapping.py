"""
Create a mapping from mutations to pubmed article ids
"""

import pandas as pd
from sparse_utils import save_pmid_mapping

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

# store mutation -> pmid mapping as gzip-compressed json
save_pmid_mapping(mutation_pmids, snek.output[0])
