# AGENTS.md

## Project Overview

This is a Snakemake pipeline for processing [PubTator 3](https://www.ncbi.nlm.nih.gov/research/pubtator/) biomedical annotation data. It downloads PubTator annotation files, filters them, and generates co-occurrence matrices for genes, diseases, chemicals, and mutations.

## Commands

```bash
# Run the full pipeline (4 parallel jobs)
uv run snakemake -j4 --configfile config/config.yml

# Dry-run to see what will be executed
uv run snakemake -n --configfile config/config.yml

# Run a specific rule
uv run snakemake -j1 --configfile config/config.yml <rule_name>

# Force re-run a rule (even if outputs exist)
uv run snakemake -j1 --configfile config/config.yml -f <rule_name>
```

## Architecture

**Pipeline flow:**
1. `download_data` - Downloads raw annotation files from PubTator 3 FTP
2. `filter_data` - Filters each annotation type (removes problematic entries, low-frequency mentions/concepts)
3. `create_*_pmid_mapping` - Creates JSON mappings of entity IDs to PubMed article IDs
4. `create_*_comat` - Generates co-occurrence matrices from PMID mappings
5. `summarize_counts` / `count_genes` - Generates count statistics

**Annotation types:** cellline, chemical, disease, gene, mutation, species

**Data formats:**
- Raw input: gzip-compressed TSV (no header)
- Intermediate/output: Feather files (for DataFrames), JSON (for ID→PMID mappings)

**Scripts pattern:** All scripts in `scripts/` use `snakemake` global variable (aliased to `snek`) to access:
- `snek.input` - Input file paths
- `snek.output` - Output file paths
- `snek.config` - Config dictionary
- `snek.params` - Rule parameters

## Configuration

Copy `config/config.example.yml` to `config/config.yml` and set:
- `output_dir`: Where to store all pipeline outputs
- `filtering.*`: Thresholds for data filtering (max annotations per article, min frequency for mentions/concepts)

## Data Dependencies

The `data/` directory contains static input files:
- `symbols.txt` - Human gene symbols from genenames.org
- `ensgenes.txt` - Human Ensembl gene IDs from genenames.org
