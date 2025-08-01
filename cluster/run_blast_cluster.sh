#!/bin/bash
#SBATCH --job-name=syn3a_blast
#SBATCH --output=blast_%j.out
#SBATCH --error=blast_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G

# Load Python module (adjust for your cluster)
module load python/3.9

# Create virtual environment
python -m venv blast_env
source blast_env/bin/activate

# Install dependencies
pip install -r requirements.txt

# Run BLAST pipeline with many workers
python run_pipeline.py data/syn3A_proteins.fasta \
    --workers 32 \
    --output-dir cluster_results \
    2>&1 | tee blast_cluster.log

echo "Job completed at $(date)"