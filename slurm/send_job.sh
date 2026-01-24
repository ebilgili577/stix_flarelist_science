#!/bin/bash
#SBATCH --partition=performance
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=08:00:00


echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Running: stix-train $@"
echo "Start time: $(date)"
echo "=========================================="

# Change to project directory (use absolute path from script location)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_DIR"
echo "Working directory: $PROJECT_DIR"

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

