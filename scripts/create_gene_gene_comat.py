"""
Creates a gene-gene co-occurrence matrix using sparse matrix multiplication.

Outputs:
    - gene_gene_counts.npz: Raw co-occurrence counts (upper triangle)
    - gene_gene_dice.npz: Dice coefficient scores (upper triangle)
    - gene_gene_npmi.npz: NPMI scores (upper triangle)
    - gene_gene_labels.json: Entity labels (shared across all matrices)
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

# load gene/pmid mapping
gene_pmids = load_pmid_mapping(snek.input[0])

print(f"Building incidence matrix for {len(gene_pmids)} genes...")

# build sparse incidence matrix (PMIDs x genes)
incidence, entities, all_pmids = build_incidence_matrix(gene_pmids)

print(f"Computing co-occurrence matrix via sparse multiplication...")

# compute co-occurrence via matrix multiplication
comat = compute_cooccurrence(incidence)

print(f"Computing normalized scores (Dice, NPMI)...")

# compute normalized scores
scores = compute_scores(comat, total_docs=len(all_pmids))

# save all matrices (upper triangle only for symmetric)
output_dir = Path(snek.output[0]).parent
prefix = getattr(snek.params, "prefix", "gene_gene")
save_cooccurrence_matrices(
    scores=scores,
    labels=entities,
    output_dir=output_dir,
    prefix=prefix,
    symmetric=True,
)

print(f"Saved {prefix} co-occurrence matrices to {output_dir}")
