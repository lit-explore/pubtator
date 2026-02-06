"""
Creates a chemical-chemical co-occurrence matrix using sparse matrix multiplication.

Outputs:
    - chemical_chemical_counts.npz: Raw co-occurrence counts (upper triangle)
    - chemical_chemical_dice.npz: Dice coefficient scores (upper triangle)
    - chemical_chemical_npmi.npz: NPMI scores (upper triangle)
    - chemical_chemical_labels.json: Entity labels (shared across all matrices)
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

# load chemical/pmid mapping
unfiltered_chemical_pmids = load_pmid_mapping(snek.input[0])

# filter chemicals with less than N citations to keep matrix size manageable
concept_id_min_freq = snek.config["filtering"]["chem_comat_concept_id_min_freq"]

chemical_pmids = {
    mesh_id: pmids
    for mesh_id, pmids in unfiltered_chemical_pmids.items()
    if len(pmids) > concept_id_min_freq
}

print(
    f"Filtered to {len(chemical_pmids)}/{len(unfiltered_chemical_pmids)} chemicals "
    f"with > {concept_id_min_freq} citations"
)

print(f"Building incidence matrix for {len(chemical_pmids)} chemicals...")

# build sparse incidence matrix (PMIDs x chemicals)
incidence, entities, all_pmids = build_incidence_matrix(chemical_pmids)

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
    prefix="chemical_chemical",
    symmetric=True,
)

print(f"Saved chemical-chemical co-occurrence matrices to {output_dir}")
