# BLAST Pipeline for Syn3A Bacterial Homology Analysis

## Overview

This is a comprehensive, production-ready BLAST pipeline for analyzing sequence homology between all 438 syn3A proteins and bacterial proteins (excluding mycoplasma). The pipeline performs true sequence-based homology searches using NCBI's BLAST service.

## Features

- ✅ **True sequence homology** using BLASTP against NCBI nr database
- ✅ **Automated mycoplasma exclusion** using taxonomic filters
- ✅ **Parallel processing** with rate limiting for NCBI compliance
- ✅ **Comprehensive result parsing** with statistical analysis
- ✅ **Robust error handling** and progress logging
- ✅ **Configurable parameters** for E-value, identity, coverage thresholds
- ✅ **Resume capability** for interrupted runs
- ✅ **Production-ready** with proper logging and monitoring

## Directory Structure

```
blast_pipeline/
├── README.md                          # This file
├── blast_config.json                  # Configuration file
├── scripts/
│   ├── setup_blast_environment.sh     # Environment setup
│   ├── run_blast_pipeline.py          # Main pipeline script
│   └── prepare_queries.py             # Query preparation utility
├── databases/
│   ├── bacterial_taxids.json          # Bacterial taxonomy information
│   └── logs/                          # Database setup logs
├── queries/                           # Input protein sequences (FASTA files)
├── results/                           # BLAST results and analysis
│   ├── blast_pipeline_results.csv     # Comprehensive results
│   ├── blast_pipeline_summary.json    # Summary statistics
│   ├── protein_*_blast.xml           # Individual BLAST XML results
│   └── protein_*_hits.csv            # Parsed hits for each protein
└── logs/                              # Pipeline execution logs
```

## Installation & Setup

### 1. Run Environment Setup

```bash
cd blast_pipeline
./scripts/setup_blast_environment.sh
```

This will:
- Install BLAST+ if not present
- Check Python dependencies
- Download bacterial taxonomy information
- Create configuration files
- Test BLAST installation

### 2. Prepare Query Sequences

```bash
# Create individual protein files (recommended for large runs)
./scripts/prepare_queries.py ../syn3A_proteins.fasta --output-dir queries --individual

# OR create chunked files (for batch processing)
./scripts/prepare_queries.py ../syn3A_proteins.fasta --output-dir queries --chunk-size 10
```

### 3. Configure Pipeline (Optional)

Edit `blast_config.json` to adjust parameters:

```json
{
  "blast_parameters": {
    "evalue": 0.01,           # E-value threshold
    "max_target_seqs": 100,   # Maximum hits per query
    "word_size": 6            # BLAST word size
  },
  "filtering": {
    "min_identity": 30.0,     # Minimum identity percentage
    "min_coverage": 50.0,     # Minimum query coverage
    "max_evalue": 0.01        # Maximum E-value for filtering
  },
  "resources": {
    "max_parallel_jobs": 4    # Parallel jobs (respect NCBI limits)
  }
}
```

## Usage

### Test Run (Recommended First)

```bash
# Test with first 3 proteins
./scripts/run_blast_pipeline.py --test

# Check test results
ls results/
```

### Full Pipeline Run

```bash
# Run complete analysis
./scripts/run_blast_pipeline.py

# With custom parameters
./scripts/run_blast_pipeline.py --workers 2 --config blast_config.json
```

### Resume Interrupted Run

The pipeline automatically resumes from where it left off. Simply rerun:

```bash
./scripts/run_blast_pipeline.py
```

## BLAST Parameters

### Search Parameters
- **Program**: BLASTP (protein-protein BLAST)
- **Database**: NCBI nr (non-redundant protein database)
- **E-value threshold**: 0.01 (configurable)
- **Maximum targets**: 100 hits per query
- **Word size**: 6 (for sensitivity)
- **Matrix**: BLOSUM62 (standard for proteins)

### Filtering Criteria
- **Taxonomic exclusion**: `NOT mycoplasma[Organism]`
- **Minimum identity**: 30% (configurable)
- **Minimum coverage**: 50% of query length (configurable)
- **Maximum E-value**: 0.01 (configurable)

### Rate Limiting
- **Maximum parallel jobs**: 4 (NCBI recommended)
- **Delay between queries**: 10-15 seconds
- **Timeout per query**: 15 minutes maximum
- **Retry logic**: Automatic retry for failed queries

## Output Files

### Primary Results

1. **`blast_pipeline_results.csv`** - Comprehensive results for all proteins
   ```csv
   query_name,query_def,query_len,status,hit_count,has_bacterial_homologs,best_hit_evalue,best_hit_identity,best_hit_description
   ```

2. **`blast_pipeline_summary.json`** - Statistical summary
   ```json
   {
     "total_queries": 438,
     "successful_searches": 420,
     "queries_with_bacterial_hits": 285,
     "percentage_with_hits": 67.9
   }
   ```

### Individual Results

3. **`protein_XXX_blast.xml`** - Raw BLAST XML output for each protein
4. **`protein_XXX_hits.csv`** - Parsed hits with E-values, identity, coverage

### Analysis Files

5. **Pipeline logs** - Detailed execution logs with timestamps
6. **Configuration backup** - Parameters used for the analysis

## Expected Runtime

### Time Estimates
- **Per protein**: 2-5 minutes average (including NCBI queue time)
- **Full 438 proteins**: 15-35 hours total
- **Test run (3 proteins)**: 10-15 minutes

### Resource Usage
- **CPU**: Low (waiting for NCBI responses)
- **Memory**: <1 GB
- **Disk**: ~500 MB for all results
- **Network**: Moderate (BLAST submissions and results)

## Quality Control

### Validation Checks
- ✅ Sequence integrity validation
- ✅ BLAST submission success verification
- ✅ Result parsing validation
- ✅ Taxonomic filter verification
- ✅ Statistical threshold enforcement

### Error Handling
- 🔄 Automatic retry for network failures
- 📝 Comprehensive error logging
- ⏸️ Resume capability for interrupted runs
- 🚨 Alert system for persistent failures

## Interpreting Results

### Homology Determination
A protein is considered to have bacterial homologs if:
1. **E-value** ≤ 0.01
2. **Identity** ≥ 30%
3. **Coverage** ≥ 50% of query length
4. **At least one hit** from non-mycoplasma bacteria

### Confidence Levels
- **High confidence**: E-value < 1e-10, Identity > 50%
- **Medium confidence**: E-value < 1e-5, Identity > 40%
- **Low confidence**: E-value < 0.01, Identity > 30%

### Expected Results
Based on functional predictions:
- **~280 proteins** expected to have bacterial homologs (63%)
- **~160 proteins** expected to be mycoplasma-specific (37%)
- **Universal proteins** (DnaA, ribosomal, etc.) should show strong homology
- **Hypothetical proteins** results will be most informative

## Troubleshooting

### Common Issues

1. **BLAST+ not found**
   ```bash
   # Install BLAST+
   brew install blast  # macOS
   sudo apt-get install ncbi-blast+  # Ubuntu
   ```

2. **Network timeouts**
   - Check internet connection
   - Reduce parallel workers: `--workers 1`
   - Increase timeout in config

3. **No query files found**
   ```bash
   # Prepare queries first
   ./scripts/prepare_queries.py ../syn3A_proteins.fasta --output-dir queries --individual
   ```

4. **Rate limiting errors**
   - Reduce parallel workers
   - Increase delays between queries
   - Use `--test` mode first

### Support Files

- **Setup log**: `logs/setup_YYYYMMDD_HHMMSS.log`
- **Pipeline log**: `logs/blast_pipeline_TIMESTAMP.log`
- **Error reports**: Individual protein error messages in logs

## Citation

If you use this pipeline in your research, please cite:

```
BLAST Pipeline for Syn3A Bacterial Homology Analysis
Developed for minimal genome functional analysis
https://github.com/your-repo/syn3a-blast-pipeline
```

## License

This pipeline is provided as-is for academic research purposes. BLAST+ and NCBI services have their own usage policies and restrictions.

---

**Pipeline Version**: 1.0  
**Last Updated**: July 30, 2025  
**Compatible with**: BLAST+ 2.10+, Python 3.7+  
**NCBI Compliance**: Follows NCBI usage guidelines and rate limits