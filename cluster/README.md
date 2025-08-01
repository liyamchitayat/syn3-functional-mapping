# Cluster Deployment

This directory contains scripts for running the pipeline on HPC clusters.

## Usage

1. Edit `run_blast_cluster.sh` for your cluster's scheduler (SLURM, PBS, etc.)
2. Adjust resources (CPUs, memory, time limit)
3. Submit job:
   ```bash
   sbatch run_blast_cluster.sh
   ```

## Performance Tips

- Use 32-64 workers for optimal performance
- Each worker needs ~50MB RAM
- No GPUs required (CPU-only)
- Ensure stable internet connection for NCBI access

## Example Times

- 8 workers: ~55 hours
- 32 workers: ~14 hours  
- 64 workers: ~7 hours