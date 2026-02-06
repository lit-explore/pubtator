"""
Creates a mutation-disease co-occurrence matrix using sparse matrix multiplication.

Outputs:
    - mutation_disease_counts.npz: Raw co-occurrence counts
    - mutation_disease_dice.npz: Dice coefficient scores
    - mutation_disease_npmi.npz: NPMI scores
    - mutation_disease_labels.json: Row (mutation) and column (disease) labels
"""

from pathlib import Path

from sparse_utils import (
    build_cross_cooccurrence,
    compute_cross_scores,
    load_pmid_mapping,
    save_cross_cooccurrence_matrices,
)

snek = snakemake

# load mutation and disease pmid mappings
mutation_pmids = load_pmid_mapping(snek.input[0])
disease_pmids = load_pmid_mapping(snek.input[1])

print(f"Building co-occurrence matrix for {len(mutation_pmids)} mutations x {len(disease_pmids)} diseases...")

# build co-occurrence matrix via sparse multiplication
comat, mutation_ids, disease_ids, total_docs = build_cross_cooccurrence(mutation_pmids, disease_pmids)

print(f"Computing normalized scores (Dice, NPMI)...")

# compute normalized scores
scores = compute_cross_scores(comat, mutation_pmids, disease_pmids, total_docs)

# save all matrices
output_dir = Path(snek.output[0]).parent
save_cross_cooccurrence_matrices(
    scores=scores,
    row_labels=mutation_ids,
    col_labels=disease_ids,
    output_dir=output_dir,
    prefix="mutation_disease",
)

print(f"Saved mutation-disease co-occurrence matrices to {output_dir}")
