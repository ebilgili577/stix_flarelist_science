#!/usr/bin/env bash
#SBATCH --partition=performance,p4500,p3080
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00

echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Running: python -m stix_train.cli $*"
echo "Start time: $(date)"
echo "=========================================="

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_DIR"
echo "Working directory: $PROJECT_DIR"
echo "Node: $(hostname)"

eval "$(conda shell.bash hook)"
conda activate stix

echo "[ENV] CONDA_PREFIX=${CONDA_PREFIX:-<unset>}"
echo "[ENV] python=$(command -v python || true)"
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi || true
fi

python -m stix_train.cli "$@"
EXIT_CODE=$?

echo "=========================================="
echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"
echo "=========================================="
exit $EXIT_CODE
