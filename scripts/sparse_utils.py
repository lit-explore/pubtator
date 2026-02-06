"""
Sparse matrix utilities for PubTator co-occurrence computation.

This module provides efficient sparse matrix operations for building and
analyzing entity co-occurrence matrices from PMID mappings.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix, diags, load_npz, save_npz, triu


def numpy_encoder(obj: object) -> int | list:
    """
    JSON encoder to convert numpy int64 elements to generic Python ints.

    Usage:
        json.dumps(data, default=numpy_encoder)
    """
    if isinstance(obj, np.generic):
        return obj.item()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def load_pmid_mapping(path: str | Path) -> dict[str, list[int]]:
    """Load PMID mapping from JSON file (supports .json and .json.gz)."""
    path = Path(path)
    if path.suffix == ".gz" or str(path).endswith(".json.gz"):
        with gzip.open(path, "rt", encoding="utf-8") as fp:
            return json.load(fp)
    else:
        with open(path, "rt", encoding="utf-8") as fp:
            return json.load(fp)


def save_pmid_mapping(data: dict[str, list[int]], path: str | Path) -> None:
    """Save PMID mapping to gzip-compressed JSON file."""
    path = Path(path)
    with gzip.open(path, "wt", encoding="utf-8") as fp:
        json.dump(data, fp, default=numpy_encoder)


def build_incidence_matrix(
    entity_pmids: dict[str, list[int]],
) -> tuple[csr_matrix, list[str], list[int]]:
    """
    Build sparse PMID x Entity incidence matrix.

    Args:
        entity_pmids: Mapping from entity IDs to lists of PMIDs

    Returns:
        Tuple of (incidence_matrix, entity_ids, pmid_list)
        - incidence_matrix: Sparse binary matrix (n_pmids x n_entities)
        - entity_ids: List of entity IDs (column labels)
        - pmid_list: List of PMIDs (row labels)
    """
    entities = list(entity_pmids.keys())
    entity_to_idx = {e: i for i, e in enumerate(entities)}

    # Collect all unique PMIDs
    all_pmids = sorted(set(pmid for pmids in entity_pmids.values() for pmid in pmids))
    pmid_to_idx = {p: i for i, p in enumerate(all_pmids)}

    # Build sparse matrix using COO format data
    rows: list[int] = []
    cols: list[int] = []

    for entity, pmids in entity_pmids.items():
        col = entity_to_idx[entity]
        for pmid in pmids:
            rows.append(pmid_to_idx[pmid])
            cols.append(col)

    data = np.ones(len(rows), dtype=np.int8)
    incidence = csr_matrix(
        (data, (rows, cols)), shape=(len(all_pmids), len(entities)), dtype=np.int32
    )

    return incidence, entities, all_pmids


def compute_cooccurrence(incidence: csr_matrix) -> csr_matrix:
    """
    Compute co-occurrence matrix via sparse matrix multiplication.

    Args:
        incidence: Sparse PMID x Entity incidence matrix

    Returns:
        Sparse Entity x Entity co-occurrence matrix.
        Diagonal contains entity frequencies (number of PMIDs per entity).
        Off-diagonal contains co-occurrence counts.
    """
    return (incidence.T @ incidence).tocsr()


def compute_scores(
    comat: csr_matrix, total_docs: int
) -> dict[str, csr_matrix]:
    """
    Compute normalized co-occurrence scores.

    Args:
        comat: Co-occurrence matrix (diagonal = entity frequencies)
        total_docs: Total number of documents in corpus

    Returns:
        Dict with 'counts', 'dice', and 'npmi' sparse matrices
    """
    # Entity frequencies from diagonal
    freqs = np.array(comat.diagonal(), dtype=np.float64).flatten()

    # Get non-zero indices
    rows, cols = comat.nonzero()
    counts = np.array(comat[rows, cols], dtype=np.float64).flatten()

    # Dice coefficient: 2 * |A ∩ B| / (|A| + |B|)
    freq_sums = freqs[rows] + freqs[cols]
    dice_data = np.divide(
        2.0 * counts, freq_sums, out=np.zeros_like(counts), where=freq_sums > 0
    )
    dice = csr_matrix((dice_data, (rows, cols)), shape=comat.shape)

    # NPMI: PMI / -log(P(A,B))
    # PMI = log(P(A,B) / (P(A) * P(B)))
    # P(A,B) = cooccur / N, P(A) = freq_a / N, P(B) = freq_b / N
    p_joint = counts / total_docs
    p_marginal_prod = (freqs[rows] / total_docs) * (freqs[cols] / total_docs)

    # Avoid log(0) and division by zero
    with np.errstate(divide="ignore", invalid="ignore"):
        pmi = np.log(p_joint / p_marginal_prod)
        npmi_data = pmi / -np.log(p_joint)

    # Replace NaN/inf with 0
    npmi_data = np.nan_to_num(npmi_data, nan=0.0, posinf=0.0, neginf=0.0)
    npmi = csr_matrix((npmi_data, (rows, cols)), shape=comat.shape)

    return {"counts": comat, "dice": dice, "npmi": npmi}


def extract_upper_triangle(matrix: csr_matrix) -> csr_matrix:  # type: ignore[return-value]
    """
    Extract upper triangular portion of a symmetric matrix.

    Includes diagonal (entity frequencies).
    """
    return triu(matrix, k=0, format="csr")


def reconstruct_symmetric(upper: csr_matrix) -> csr_matrix:
    """
    Reconstruct full symmetric matrix from upper triangular portion.
    """
    # upper + upper.T - diagonal (to avoid doubling diagonal)
    diagonal = diags(upper.diagonal(), format="csr", dtype=upper.dtype)
    return (upper + upper.T - diagonal).tocsr()


def save_sparse_matrix(
    matrix: csr_matrix,
    labels: list[str],
    output_prefix: str | Path,
) -> None:
    """
    Save sparse matrix with labels.

    Creates two files:
        - {output_prefix}.npz: Sparse matrix in scipy format
        - {output_prefix}_labels.json: Row/column labels
    """
    output_prefix = Path(output_prefix)
    save_npz(output_prefix.with_suffix(".npz"), matrix)

    with open(f"{output_prefix}_labels.json", "w", encoding="utf-8") as fp:
        json.dump(labels, fp)


def load_sparse_matrix(input_prefix: str | Path) -> tuple[csr_matrix, list[str]]:
    """
    Load sparse matrix with labels.

    Args:
        input_prefix: Path prefix (without .npz extension)

    Returns:
        Tuple of (matrix, labels)
    """
    input_prefix = Path(input_prefix)
    matrix = load_npz(input_prefix.with_suffix(".npz"))

    with open(f"{input_prefix}_labels.json", "r", encoding="utf-8") as fp:
        labels = json.load(fp)

    return matrix, labels


def save_cooccurrence_matrices(
    scores: dict[str, csr_matrix],
    labels: list[str],
    output_dir: str | Path,
    prefix: str,
    symmetric: bool = False,
) -> None:
    """
    Save all co-occurrence score matrices.

    Args:
        scores: Dict with 'counts', 'dice', 'npmi' matrices
        labels: Entity labels (row/column names)
        output_dir: Output directory
        prefix: Filename prefix (e.g., 'gene_gene')
        symmetric: If True, store only upper triangle
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save labels once (shared across all score types)
    labels_path = output_dir / f"{prefix}_labels.json"
    with open(labels_path, "w", encoding="utf-8") as fp:
        json.dump(labels, fp)

    # Save each score matrix
    for score_type, matrix in scores.items():
        if symmetric:
            matrix = extract_upper_triangle(matrix)

        matrix_path = output_dir / f"{prefix}_{score_type}.npz"
        save_npz(matrix_path, matrix)


def build_cross_cooccurrence(
    entity1_pmids: dict[str, list[int]],
    entity2_pmids: dict[str, list[int]],
) -> tuple[csr_matrix, list[str], list[str], int]:
    """
    Build co-occurrence matrix between two different entity types.

    Args:
        entity1_pmids: First entity type PMID mapping (rows)
        entity2_pmids: Second entity type PMID mapping (columns)

    Returns:
        Tuple of (comat, entity1_ids, entity2_ids, total_docs)
    """
    # Build incidence matrices for each entity type
    entities1 = list(entity1_pmids.keys())
    entities2 = list(entity2_pmids.keys())

    # Find common PMIDs
    pmids1 = set(pmid for pmids in entity1_pmids.values() for pmid in pmids)
    pmids2 = set(pmid for pmids in entity2_pmids.values() for pmid in pmids)
    all_pmids = sorted(pmids1 | pmids2)
    pmid_to_idx = {p: i for i, p in enumerate(all_pmids)}

    # Build incidence matrix for entity type 1
    rows1, cols1 = [], []
    for i, pmids in enumerate(entity1_pmids.values()):
        for pmid in pmids:
            rows1.append(pmid_to_idx[pmid])
            cols1.append(i)

    inc1 = csr_matrix(
        (np.ones(len(rows1), dtype=np.int32), (rows1, cols1)),
        shape=(len(all_pmids), len(entities1)),
    )

    # Build incidence matrix for entity type 2
    rows2, cols2 = [], []
    for i, pmids in enumerate(entity2_pmids.values()):
        for pmid in pmids:
            rows2.append(pmid_to_idx[pmid])
            cols2.append(i)

    inc2 = csr_matrix(
        (np.ones(len(rows2), dtype=np.int32), (rows2, cols2)),
        shape=(len(all_pmids), len(entities2)),
    )

    # Co-occurrence: entity1 x entity2
    comat = (inc1.T @ inc2).tocsr()

    return comat, entities1, entities2, len(all_pmids)


def compute_cross_scores(
    comat: csr_matrix,
    entity1_pmids: dict[str, list[int]],
    entity2_pmids: dict[str, list[int]],
    total_docs: int,
) -> dict[str, csr_matrix]:
    """
    Compute normalized scores for cross-entity co-occurrence matrix.

    Args:
        comat: Entity1 x Entity2 co-occurrence matrix
        entity1_pmids: Row entity PMID mapping (for frequencies)
        entity2_pmids: Column entity PMID mapping (for frequencies)
        total_docs: Total number of documents

    Returns:
        Dict with 'counts', 'dice', 'npmi' sparse matrices
    """
    # Entity frequencies
    freqs1 = np.array([len(pmids) for pmids in entity1_pmids.values()], dtype=np.float64)
    freqs2 = np.array([len(pmids) for pmids in entity2_pmids.values()], dtype=np.float64)

    # Get non-zero indices
    rows, cols = comat.nonzero()
    counts = np.array(comat[rows, cols], dtype=np.float64).flatten()

    # Dice coefficient
    freq_sums = freqs1[rows] + freqs2[cols]
    dice_data = np.divide(
        2.0 * counts, freq_sums, out=np.zeros_like(counts), where=freq_sums > 0
    )
    dice = csr_matrix((dice_data, (rows, cols)), shape=comat.shape)

    # NPMI
    p_joint = counts / total_docs
    p_marginal_prod = (freqs1[rows] / total_docs) * (freqs2[cols] / total_docs)

    with np.errstate(divide="ignore", invalid="ignore"):
        pmi = np.log(p_joint / p_marginal_prod)
        npmi_data = pmi / -np.log(p_joint)

    npmi_data = np.nan_to_num(npmi_data, nan=0.0, posinf=0.0, neginf=0.0)
    npmi = csr_matrix((npmi_data, (rows, cols)), shape=comat.shape)

    return {"counts": comat, "dice": dice, "npmi": npmi}


def save_cross_cooccurrence_matrices(
    scores: dict[str, csr_matrix],
    row_labels: list[str],
    col_labels: list[str],
    output_dir: str | Path,
    prefix: str,
) -> None:
    """
    Save cross-entity co-occurrence matrices (non-symmetric).

    Args:
        scores: Dict with 'counts', 'dice', 'npmi' matrices
        row_labels: Row entity labels
        col_labels: Column entity labels
        output_dir: Output directory
        prefix: Filename prefix (e.g., 'gene_disease')
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save labels
    labels_path = output_dir / f"{prefix}_labels.json"
    with open(labels_path, "w", encoding="utf-8") as fp:
        json.dump({"rows": row_labels, "cols": col_labels}, fp)

    # Save each score matrix
    for score_type, matrix in scores.items():
        matrix_path = output_dir / f"{prefix}_{score_type}.npz"
        save_npz(matrix_path, matrix)
