"""
Unit tests for sparse_utils module.
"""

from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from sparse_utils import (
    build_cross_cooccurrence,
    build_incidence_matrix,
    compute_cooccurrence,
    compute_cross_scores,
    compute_scores,
    extract_upper_triangle,
    load_pmid_mapping,
    numpy_encoder,
    reconstruct_symmetric,
    save_pmid_mapping,
)


class TestNumpyEncoder:
    """Tests for numpy_encoder function."""

    def test_encodes_numpy_int64(self) -> None:
        result = numpy_encoder(np.int64(42))
        assert result == 42
        assert isinstance(result, int)

    def test_encodes_numpy_float64(self) -> None:
        result = numpy_encoder(np.float64(3.14))
        assert result == pytest.approx(3.14)
        assert isinstance(result, float)

    def test_raises_for_non_numpy(self) -> None:
        with pytest.raises(TypeError, match="not JSON serializable"):
            numpy_encoder({"key": "value"})


class TestPmidMappingIO:
    """Tests for PMID mapping save/load functions."""

    def test_save_and_load_roundtrip(
        self, tmp_path: Path, sample_entity_pmids: dict[str, list[int]]
    ) -> None:
        filepath = tmp_path / "test.json.gz"
        save_pmid_mapping(sample_entity_pmids, filepath)

        loaded = load_pmid_mapping(filepath)
        assert loaded == sample_entity_pmids

    def test_load_gzip_file(self, tmp_json_gz: Path) -> None:
        loaded = load_pmid_mapping(tmp_json_gz)
        assert "entity_a" in loaded
        assert loaded["entity_a"] == [1, 2, 3, 4, 5]


class TestBuildIncidenceMatrix:
    """Tests for build_incidence_matrix function."""

    def test_builds_correct_shape(
        self, sample_entity_pmids: dict[str, list[int]]
    ) -> None:
        incidence, entities, pmids = build_incidence_matrix(sample_entity_pmids)

        # 4 entities
        assert len(entities) == 4
        # 12 unique PMIDs (1-12)
        assert len(pmids) == 12
        # Shape should be (n_pmids, n_entities)
        assert incidence.shape == (12, 4)

    def test_incidence_values_are_binary(
        self, sample_entity_pmids: dict[str, list[int]]
    ) -> None:
        incidence, _, _ = build_incidence_matrix(sample_entity_pmids)
        unique_values = set(incidence.data)
        assert unique_values == {1}

    def test_entity_order_preserved(
        self, sample_entity_pmids: dict[str, list[int]]
    ) -> None:
        incidence, entities, _ = build_incidence_matrix(sample_entity_pmids)
        expected_entities = list(sample_entity_pmids.keys())
        assert entities == expected_entities


class TestComputeCooccurrence:
    """Tests for compute_cooccurrence function."""

    def test_diagonal_contains_frequencies(
        self, sample_entity_pmids: dict[str, list[int]]
    ) -> None:
        incidence, entities, _ = build_incidence_matrix(sample_entity_pmids)
        comat = compute_cooccurrence(incidence)

        # Diagonal should contain entity frequencies
        for i, entity in enumerate(entities):
            expected_freq = len(sample_entity_pmids[entity])
            assert comat[i, i] == expected_freq

    def test_cooccurrence_counts_correct(
        self, sample_entity_pmids: dict[str, list[int]]
    ) -> None:
        incidence, entities, _ = build_incidence_matrix(sample_entity_pmids)
        comat = compute_cooccurrence(incidence)

        # entity_a and entity_b share PMIDs 3, 4, 5 -> count = 3
        idx_a = entities.index("entity_a")
        idx_b = entities.index("entity_b")
        assert comat[idx_a, idx_b] == 3

        # entity_a and entity_d share no PMIDs -> count = 0
        idx_d = entities.index("entity_d")
        assert comat[idx_a, idx_d] == 0

    def test_matrix_is_symmetric(
        self, sample_entity_pmids: dict[str, list[int]]
    ) -> None:
        incidence, _, _ = build_incidence_matrix(sample_entity_pmids)
        comat = compute_cooccurrence(incidence)

        diff = comat - comat.T
        assert diff.nnz == 0  # No non-zero differences


class TestComputeScores:
    """Tests for compute_scores function."""

    def test_returns_all_score_types(
        self, sample_entity_pmids: dict[str, list[int]]
    ) -> None:
        incidence, _, pmids = build_incidence_matrix(sample_entity_pmids)
        comat = compute_cooccurrence(incidence)
        scores = compute_scores(comat, total_docs=len(pmids))

        assert "counts" in scores
        assert "dice" in scores
        assert "npmi" in scores

    def test_dice_bounded_zero_one(
        self, sample_entity_pmids: dict[str, list[int]]
    ) -> None:
        incidence, _, pmids = build_incidence_matrix(sample_entity_pmids)
        comat = compute_cooccurrence(incidence)
        scores = compute_scores(comat, total_docs=len(pmids))

        dice_values = scores["dice"].data
        assert all(0 <= v <= 1 for v in dice_values)

    def test_npmi_bounded(self, sample_entity_pmids: dict[str, list[int]]) -> None:
        incidence, _, pmids = build_incidence_matrix(sample_entity_pmids)
        comat = compute_cooccurrence(incidence)
        scores = compute_scores(comat, total_docs=len(pmids))

        npmi_values = scores["npmi"].data
        # NPMI should be bounded [-1, 1] for most cases
        assert all(-1.01 <= v <= 1.01 for v in npmi_values)


class TestUpperTriangle:
    """Tests for upper triangle extraction and reconstruction."""

    def test_extract_upper_preserves_diagonal(self) -> None:
        data = np.array([[5, 2, 0], [2, 3, 1], [0, 1, 4]])
        matrix = csr_matrix(data)

        upper = extract_upper_triangle(matrix)

        # Check diagonal preserved
        assert upper[0, 0] == 5
        assert upper[1, 1] == 3
        assert upper[2, 2] == 4

    def test_extract_upper_zeros_lower(self) -> None:
        data = np.array([[5, 2, 0], [2, 3, 1], [0, 1, 4]])
        matrix = csr_matrix(data)

        upper = extract_upper_triangle(matrix)

        # Lower triangle should be zero
        assert upper[1, 0] == 0
        assert upper[2, 0] == 0
        assert upper[2, 1] == 0

    def test_reconstruct_symmetric_roundtrip(self) -> None:
        data = np.array([[5, 2, 0], [2, 3, 1], [0, 1, 4]])
        matrix = csr_matrix(data)

        upper = extract_upper_triangle(matrix)
        reconstructed = reconstruct_symmetric(upper)

        diff = matrix - reconstructed
        assert diff.nnz == 0


class TestCrossCooccurrence:
    """Tests for cross-entity co-occurrence computation."""

    def test_cross_cooccurrence_shape(
        self,
        sample_gene_pmids: dict[str, list[int]],
        sample_disease_pmids: dict[str, list[int]],
    ) -> None:
        comat, genes, diseases, _ = build_cross_cooccurrence(
            sample_gene_pmids, sample_disease_pmids
        )

        assert comat.shape == (len(genes), len(diseases))
        assert len(genes) == 3
        assert len(diseases) == 3

    def test_cross_cooccurrence_values(
        self,
        sample_gene_pmids: dict[str, list[int]],
        sample_disease_pmids: dict[str, list[int]],
    ) -> None:
        comat, genes, diseases, _ = build_cross_cooccurrence(
            sample_gene_pmids, sample_disease_pmids
        )

        # Gene 1234 (PMIDs 100-104) and Disease 1 (PMIDs 100,101,105,106,107)
        # Overlap: 100, 101 -> count = 2
        gene_idx = genes.index("1234")
        disease_idx = diseases.index("MESH:D001234")
        assert comat[gene_idx, disease_idx] == 2

        # Gene 9999 (PMIDs 200-202) and Disease 3 (PMIDs 300-302)
        # No overlap -> count = 0
        gene_idx = genes.index("9999")
        disease_idx = diseases.index("MESH:D009999")
        assert comat[gene_idx, disease_idx] == 0

    def test_cross_scores_computed(
        self,
        sample_gene_pmids: dict[str, list[int]],
        sample_disease_pmids: dict[str, list[int]],
    ) -> None:
        comat, _, _, total_docs = build_cross_cooccurrence(
            sample_gene_pmids, sample_disease_pmids
        )
        scores = compute_cross_scores(
            comat, sample_gene_pmids, sample_disease_pmids, total_docs
        )

        assert "counts" in scores
        assert "dice" in scores
        assert "npmi" in scores
