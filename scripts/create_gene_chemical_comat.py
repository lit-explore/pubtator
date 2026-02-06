"""
Creates a gene-chemical co-occurrence matrix using sparse matrix multiplication.

Outputs:
    - gene_chemical_counts.npz: Raw co-occurrence counts
    - gene_chemical_dice.npz: Dice coefficient scores
    - gene_chemical_npmi.npz: NPMI scores
    - gene_chemical_labels.json: Row (gene) and column (chemical) labels
"""

from pathlib import Path

from sparse_utils import (
    build_cross_cooccurrence,
    compute_cross_scores,
    load_pmid_mapping,
    save_cross_cooccurrence_matrices,
)

snek = snakemake

# load gene and chemical pmid mappings
gene_pmids = load_pmid_mapping(snek.input[0])
chemical_pmids = load_pmid_mapping(snek.input[1])

print(f"Building co-occurrence matrix for {len(gene_pmids)} genes x {len(chemical_pmids)} chemicals...")

# build co-occurrence matrix via sparse multiplication
comat, gene_ids, chemical_ids, total_docs = build_cross_cooccurrence(gene_pmids, chemical_pmids)

print(f"Computing normalized scores (Dice, NPMI)...")

# compute normalized scores
scores = compute_cross_scores(comat, gene_pmids, chemical_pmids, total_docs)

# save all matrices
output_dir = Path(snek.output[0]).parent
prefix = getattr(snek.params, "prefix", "gene_chemical")
save_cross_cooccurrence_matrices(
    scores=scores,
    row_labels=gene_ids,
    col_labels=chemical_ids,
    output_dir=output_dir,
    prefix=prefix,
)

print(f"Saved {prefix} co-occurrence matrices to {output_dir}")
