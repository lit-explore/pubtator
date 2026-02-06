"""
Create a mapping from entrez gene ids to pubmed article ids
"""

import json

import pandas as pd
from utils import numpy_encoder

snek = snakemake

SPECIES: str = snek.params["species"]

# load species data
species_df = pd.read_feather(snek.input[0])

# load gene data
gene_df = pd.read_feather(snek.input[1])

# load entrez id mappings
symbol2entrez = pd.read_feather(snek.input[2])
ensgene2entrez = pd.read_feather(snek.input[3])

# create a list of all unique entrez ids associated with entries from either mapping
entrez_ids = set(symbol2entrez.entrezgene.values)
entrez_ids = entrez_ids.union(set(ensgene2entrez.entrezgene.values))

# limit to articles mentioning target species
species_pmids = set(species_df.pmid[species_df.concept_id == SPECIES].values)
gene_df = gene_df[gene_df.pmid.isin(species_pmids)]

# limit to species-specific entrez ids
gene_df = gene_df[gene_df.concept_id.isin(entrez_ids)]

# build gene -> pmid mapping using groupby (much faster than iteration)
gene_pmids = (
    gene_df.groupby("concept_id")["pmid"].apply(lambda x: x.unique().tolist()).to_dict()
)

print(f"Mapped {len(gene_pmids)} genes to PMIDs")

# store gene -> pmid mapping as json
with open(snek.output[0], "wt", encoding="utf-8") as fp:
    json.dump(gene_pmids, fp, default=numpy_encoder)
