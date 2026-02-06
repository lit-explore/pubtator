"""
Creates a gene-disease co-occurrence matrix using sparse matrix multiplication.

Outputs:
    - gene_disease_counts.npz: Raw co-occurrence counts
    - gene_disease_dice.npz: Dice coefficient scores
    - gene_disease_npmi.npz: NPMI scores
    - gene_disease_labels.json: Row (gene) and column (disease) labels
"""

from pathlib import Path

from sparse_utils import (
    build_cross_cooccurrence,
    compute_cross_scores,
    load_pmid_mapping,
    save_cross_cooccurrence_matrices,
)

snek = snakemake

# load gene and disease pmid mappings
gene_pmids = load_pmid_mapping(snek.input[0])
disease_pmids = load_pmid_mapping(snek.input[1])

print(f"Building co-occurrence matrix for {len(gene_pmids)} genes x {len(disease_pmids)} diseases...")

# build co-occurrence matrix via sparse multiplication
comat, gene_ids, disease_ids, total_docs = build_cross_cooccurrence(gene_pmids, disease_pmids)

print(f"Computing normalized scores (Dice, NPMI)...")

# compute normalized scores
scores = compute_cross_scores(comat, gene_pmids, disease_pmids, total_docs)

# save all matrices
output_dir = Path(snek.output[0]).parent
prefix = getattr(snek.params, "prefix", "gene_disease")
save_cross_cooccurrence_matrices(
    scores=scores,
    row_labels=gene_ids,
    col_labels=disease_ids,
    output_dir=output_dir,
    prefix=prefix,
)

print(f"Saved {prefix} co-occurrence matrices to {output_dir}")
