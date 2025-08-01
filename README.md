# Syn3A Functional Mapping Pipeline

BLAST-based pipeline for detecting bacterial homologs in the Mycoplasma laboratorium syn3A minimal proteome.

## Overview

This project analyzes all 438 proteins from syn3A to identify which have non-mycoplasmal bacterial homologs, helping understand the functional composition of this minimal synthetic organism.

## Pipeline Features

- Web-based BLAST (no local installation required)
- Parallel processing with configurable workers
- Automatic retry and error handling
- Progress tracking and ETA estimation
- Comprehensive result filtering

## Quick Start

### Installation

```bash
git clone <repository-url>
cd syn3_funcitonal_mapping
pip install -r requirements.txt
```

### Download syn3A proteins

```bash
# The syn3A proteins file is required but not included in git
# Download from NCBI or use provided script
python scripts/download_syn3a.py
```

### Run Analysis

```bash
# Test run (5 proteins)
python blast_pipeline/scripts/fast_blast_pipeline.py syn3A_proteins.fasta --test

# Full analysis with 6 workers
python blast_pipeline/scripts/fast_blast_pipeline.py syn3A_proteins.fasta --workers 6

# Run on cluster with more workers
python blast_pipeline/scripts/fast_blast_pipeline.py syn3A_proteins.fasta --workers 32
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