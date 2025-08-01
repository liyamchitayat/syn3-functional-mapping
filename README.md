# Syn3A Functional Mapping Pipeline

BLAST-based pipeline for detecting bacterial homologs in the Mycoplasma laboratorium syn3A minimal proteome.

## Overview

This project analyzes all 438 proteins from syn3A to identify which have non-mycoplasmal bacterial homologs, helping understand the functional composition of this minimal synthetic organism.

## Project Structure

```
syn3a-blast-pipeline/
├── src/
│   ├── pipeline/       # BLAST pipeline implementations
│   ├── analysis/       # Analysis and post-processing scripts
│   └── utils/          # Utility functions and data download
├── data/               # Input data (FASTA files)
├── results/            # Pipeline outputs
├── cluster/            # HPC/cluster scripts
├── tests/              # Test scripts
├── docs/               # Documentation
└── run_pipeline.py     # Main entry point
```

## Quick Start

### 1. Installation

```bash
git clone <repository-url>
cd syn3a-blast-pipeline
pip install -r requirements.txt
```

### 2. Download Data

```bash
# Download syn3A proteins from NCBI
python src/utils/download_syn3a.py
```

### 3. Run Analysis

```bash
# Test run (5 proteins)
python run_pipeline.py data/syn3A_proteins.fasta --test

# Full analysis with 8 workers
python run_pipeline.py data/syn3A_proteins.fasta --workers 8

# Resume from specific protein
python run_pipeline.py data/syn3A_proteins.fasta --start 100
```

## Cluster Deployment

For HPC/cluster environments:

```bash
# Example SLURM script
sbatch run_blast_cluster.sh
```

The pipeline uses web-based BLAST, so it only needs:
- Python 3.7+
- Internet connection
- No GPUs required (BLAST is CPU-based)

## Parameters

- `--workers`: Number of parallel BLAST jobs (default: 6, max: 50 for clusters)
- `--output-dir`: Output directory for results
- `--start`: Start from specific protein index
- `--max`: Maximum proteins to process

## Results

Results are saved in CSV format with:
- Protein ID and description
- Number of bacterial hits
- Best hit details (E-value, identity, coverage)
- Processing statistics

## Note on Performance

- ~3-4 proteins/hour per worker
- Full 438 proteins: ~15 hours with 8 workers
- Cluster with 32 workers: ~4 hours

🤖 Generated with [Claude Code](https://claude.ai/code)

Co-Authored-By: Claude <noreply@anthropic.com>