#!/usr/bin/env python3
"""
Interactive job submission for STIX training experiments.

Entry point for the stix-submit command.
Prompts for experiment parameters and submits a SLURM job.
"""

import subprocess
from pathlib import Path

from .config import TrainConfig, DEFAULT_HIDDEN_DIMS, EXPERIMENTS_DIR


def ask(prompt: str, default=None) -> str:
    """Prompt user for input with optional default."""
    suffix = f" [{default}]" if default is not None else ""
    response = input(f"{prompt}{suffix}: ").strip()
    return response if response else str(default) if default is not None else ""


def main():
    """Main entry point for stix-submit command."""
    # Gather experiment parameters
    exp = ask("experiment name", "exp1")
    mode = ask("mode (real/syn/mixed/finetune/all)", "all")
    col_count = ask("how many subcols?", TrainConfig.col_count)
    sidelobes_threshold = ask("sidelobes threshold (1 for no filtering)", TrainConfig.sidelobes_threshold)

    
    # Synthetic data params (only ask if needed)
    if mode != "real":
        n_samples = ask("synthetic n_samples", TrainConfig.n_samples)
        fov_big = ask("synthetic fov_big", TrainConfig.fov_big)
    else:
        n_samples = str(TrainConfig.n_samples)
        fov_big = str(TrainConfig.fov_big)
    
    # Training params
    x_fov = ask("x_fov", TrainConfig.x_fov)
    y_fov = ask("y_fov", TrainConfig.y_fov)
    lr = ask("learning rate", TrainConfig.learning_rate)
    epochs = ask("epochs", TrainConfig.epochs)
    batch = ask("batch size", TrainConfig.batch_size)
    patience = ask("early stop patience", TrainConfig.patience)
    n_runs = ask("n_runs", TrainConfig.n_runs)
    default_hidden_dims_str = ",".join(map(str, DEFAULT_HIDDEN_DIMS))
    hidden_dims = ask("hidden layer sizes (comma-separated)", default_hidden_dims_str)
    description = ask("description of experiment", "")
    
    # Create log directory for SLURM output
    log_dir = EXPERIMENTS_DIR / exp / mode
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # Save description if provided
    if description:
        desc_path = EXPERIMENTS_DIR / exp / "description.txt"
        desc_path.parent.mkdir(parents=True, exist_ok=True)
        with open(desc_path, "w") as f:
            f.write(description)
        print(f"Description saved to {desc_path}")
    
    # Build SLURM command with absolute paths
    package_root = Path(__file__).resolve().parent.parent
    slurm_script = package_root / "slurm" / "send_job.sh"
    log_path = EXPERIMENTS_DIR / exp / mode / "train_%j.log"
    cmd = [
        "sbatch",
        "-J", exp,
        "-o", str(log_path),
        str(slurm_script),
        "--experiment", exp,
        "--n-samples", str(n_samples),
        "--fov-big", str(fov_big),
        "--x-fov", str(x_fov),
        "--y-fov", str(y_fov),
        "--lr", str(lr),
        "--epochs", str(epochs),
        "--batch", str(batch),
        "--patience", str(patience),
        "--mode", mode,
        "--n-runs", str(n_runs),
        "--col-count", str(col_count),
        "--sidelobes-threshold", str(sidelobes_threshold),
        "--hidden-dims", hidden_dims,
    ]
    
    print("\nSubmitting:", " ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
