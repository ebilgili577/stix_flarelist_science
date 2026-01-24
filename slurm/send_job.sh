#!/bin/bash
#SBATCH --partition=performance
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --nodelist=calc-g-006


echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Running: stix-train $@"
echo "Start time: $(date)"
echo "=========================================="

# Change to project directory 
cd ~/stix_flarelist_science

# Initialize conda
eval "$(conda shell.bash hook)"
conda activate stix

# Run training with all passed arguments
stix-train "$@"

EXIT_CODE=$?

echo "=========================================="
echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"
echo "=========================================="

exit $EXIT_CODE
