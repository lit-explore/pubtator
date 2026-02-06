"""
Use Biothings API to construct a mapping between gene ID types.

Supports mapping from gene symbols or Ensembl gene IDs to Entrez gene IDs
for any species supported by the Biothings API.
"""

from typing import Any

import pandas as pd
from biothings_client import get_client

snek = snakemake

mg = get_client("gene")

# Get parameters from snakemake
species_id: int = snek.params["species_id"]
scope: str = snek.params["scope"]  # "symbol" or "ensembl.gene"
colname: str = snek.params["colname"]  # "symbol" or "ensgene"

# load input id list
with open(snek.input[0], "rt", encoding="utf-8") as fp:
    input_ids = [x.strip() for x in fp.readlines()]


def chunks(lst: list[Any], chunk_size: int):
    """Split list into batches."""
    for j in range(0, len(lst), chunk_size):
        yield lst[j : j + chunk_size]


batches = list(chunks(input_ids, 100))
res = []

for i, batch in enumerate(batches):
    print(f"Querying Biothings ({i + 1}/{len(batches)}) for species {species_id}...")

    mapping = mg.querymany(
        batch, scopes=scope, species=species_id, fields="entrezgene", as_dataframe=True
    )

    try:
        mapping = mapping[~mapping.entrezgene.isna()]
        res.append(mapping)
    except Exception as e:
        print(f"mapping failed for batch {i}: {e}")
        print(batch)

if not res:
    print("Warning: No mappings found!")
    df = pd.DataFrame(columns=[colname, "entrezgene"])
else:
    # drop unneeded columns
    df = pd.concat([x["entrezgene"] for x in res]).to_frame()
    df = df.reset_index().rename(columns={"query": colname})

print(f"Successfully mapped {df.shape[0]}/{len(input_ids)}...")

df.to_feather(snek.output[0])
