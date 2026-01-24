#!/bin/bash
# Interactive worker node for checking stuff
# Spawns an interactive bash session on a performance CPU node

# Get project root (directory where this script is located)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"

echo "=========================================="
echo "Requesting interactive node..."
echo "Partition: performance (CPU)"
echo "Project directory: $PROJECT_DIR"
echo "=========================================="
echo ""

# --partition=performance: CPU partition
# --cpus-per-task=8: 8 CPUs
# --mem=16G: 16GB memory
# --time=4:00:00: 4 hour time limit
# --pty: pseudo-terminal for interactive session
srun \
    --partition=performance \
    --cpus-per-task=8 \
    --mem=16G \
    --time=4:00:00 \
    --pty \
    bash -c "cd $PROJECT_DIR && eval \"\$(conda shell.bash hook)\" && conda activate stix && echo 'Conda environment activated: stix' && echo 'Working directory: $(pwd)' && echo '' && exec bash"

