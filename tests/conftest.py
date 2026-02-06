"""
Shared pytest fixtures for PubTator pipeline tests.
"""

import sys
from pathlib import Path

import pytest

# Add scripts directory to path for imports
SCRIPTS_DIR = Path(__file__).parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))


@pytest.fixture
def sample_entity_pmids() -> dict[str, list[int]]:
    """Sample entity -> PMID mapping for testing."""
    return {
        "entity_a": [1, 2, 3, 4, 5],
        "entity_b": [3, 4, 5, 6, 7],
        "entity_c": [5, 6, 7, 8, 9],
        "entity_d": [10, 11, 12],  # No overlap with entity_a
    }


@pytest.fixture
def sample_gene_pmids() -> dict[str, list[int]]:
    """Sample gene -> PMID mapping for testing."""
    return {
        "1234": [100, 101, 102, 103, 104],  # Gene 1234
        "5678": [102, 103, 104, 105, 106],  # Gene 5678 (overlaps with 1234)
        "9999": [200, 201, 202],  # Gene 9999 (no overlap)
    }


@pytest.fixture
def sample_disease_pmids() -> dict[str, list[int]]:
    """Sample disease -> PMID mapping for testing."""
    return {
        "MESH:D001234": [100, 101, 105, 106, 107],  # Disease 1
        "MESH:D005678": [102, 103, 108, 109],  # Disease 2
        "MESH:D009999": [300, 301, 302],  # Disease 3 (no overlap with genes)
    }


@pytest.fixture
def tmp_json_gz(tmp_path: Path, sample_entity_pmids: dict[str, list[int]]) -> Path:
    """Create a temporary gzip-compressed JSON file."""
    import gzip
    import json

    filepath = tmp_path / "test_mapping.json.gz"
    with gzip.open(filepath, "wt", encoding="utf-8") as fp:
        json.dump(sample_entity_pmids, fp)
    return filepath
