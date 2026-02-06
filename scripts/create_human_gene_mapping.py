"""
Use Biothings API to construct a between human gene id types
"""

from typing import Any, List

import pandas as pd
from biothings_client import get_client

snek = snakemake

mg = get_client("gene")

# load input id l ist
with open(snek.input[0], "rt", encoding="utf-8") as fp:
    input_ids = [x.strip() for x in fp.readlines()]


def chunks(lst: List[Any], num_chunks: int):
    """
    split ids into batches
    https://stackoverflow.com/a/312464/554531
    """
    for j in range(0, len(lst), num_chunks):
        yield lst[j : j + num_chunks]


batches = list(chunks(input_ids, 100))

res = []

if snek.input[0] == "data/symbols.txt":
    scopes = "symbol"
    colname = "symbol"
else:
    scopes = "ensembl.gene"
    colname = "ensgene"

for i, batch in enumerate(batches):
    print(f"Querying Biothings ({i + 1}/{len(batches)})")

    mapping = mg.querymany(
        batch, scopes=scopes, species=9606, fields="entrezgene", as_dataframe=True
    )

    try:
        mapping = mapping[~mapping.entrezgene.isna()]
        res.append(mapping)
    except Exception as e:
        print(f"mapping failed for batch {i}: {e}")
        print(batch)

# drop unneeded columns
df = pd.concat([x["entrezgene"] for x in res]).to_frame()

df = df.reset_index().rename(columns={"query": colname})

print(f"Successfully mapped {df.shape[0]}/{len(input_ids)}...")

df.to_feather(snek.output[0])
