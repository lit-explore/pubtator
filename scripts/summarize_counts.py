"""
Summarizes concept_id and mentions counts
"""

import pandas as pd

snek = snakemake

df = pd.read_feather(snek.input[0])

# count concept ids (pandas 2.0+ returns 'count' column)
concepts = df.concept_id.value_counts().reset_index()
concepts.columns = ["concept_id", "n"]

# count mentions
mentions = df.mentions.value_counts().reset_index()
mentions.columns = ["mentions", "n"]

# store results
concepts.to_feather(snek.output[0])
mentions.to_feather(snek.output[1])
