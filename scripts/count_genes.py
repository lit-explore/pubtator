"""
Counts the number of unique genes associated with each article
"""

import pandas as pd

snek = snakemake

df = pd.read_feather(snek.input[0])

counts = df.groupby("pmid").concept_id.nunique()
counts.rename("n").reset_index().to_feather(snek.output[0])
