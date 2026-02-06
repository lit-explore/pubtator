"""
Creates a disease-disease co-occurrence matrix using sparse matrix multiplication.

Outputs:
    - disease_disease_counts.npz: Raw co-occurrence counts (upper triangle)
    - disease_disease_dice.npz: Dice coefficient scores (upper triangle)
    - disease_disease_npmi.npz: NPMI scores (upper triangle)
    - disease_disease_labels.json: Entity labels (shared across all matrices)
"""

from pathlib import Path

from sparse_utils import (
    build_incidence_matrix,
    compute_cooccurrence,
    compute_scores,
    load_pmid_mapping,
    save_cooccurrence_matrices,
)

snek = snakemake

# load disease/pmid mapping
disease_pmids = load_pmid_mapping(snek.input[0])

print(f"Building incidence matrix for {len(disease_pmids)} diseases...")

# build sparse incidence matrix (PMIDs x diseases)
incidence, entities, all_pmids = build_incidence_matrix(disease_pmids)

print(f"Computing co-occurrence matrix via sparse multiplication...")

# compute co-occurrence via matrix multiplication
comat = compute_cooccurrence(incidence)

print(f"Computing normalized scores (Dice, NPMI)...")

# compute normalized scores
scores = compute_scores(comat, total_docs=len(all_pmids))

# save all matrices (upper triangle only for symmetric)
output_dir = Path(snek.output[0]).parent
save_cooccurrence_matrices(
    scores=scores,
    labels=entities,
    output_dir=output_dir,
    prefix="disease_disease",
    symmetric=True,
)

print(f"Saved disease-disease co-occurrence matrices to {output_dir}")
