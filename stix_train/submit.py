#!/usr/bin/env python3
"""
Interactive job submission for STIX training experiments.

Entry point for the stix-submit command.
Prompts for experiment parameters and submits a SLURM job.
"""

import subprocess
import os
import re
from pathlib import Path

from .config import TrainConfig, DEFAULT_HIDDEN_DIMS


def get_project_root() -> Path:
    """Get the project root directory (where setup.py is located)."""
    # Get the directory of this file
    current_file = Path(__file__).resolve()
    # Go up from stix_train/submit.py to stix_flarelist_science/
    project_root = current_file.parent.parent
    return project_root


def ask(prompt: str, default=None) -> str:
    """Prompt user for input with optional default."""
    suffix = f" [{default}]" if default is not None else ""
    response = input(f"{prompt}{suffix}: ").strip()
    return response if response else str(default) if default is not None else ""


def main():
    """Main entry point for stix-submit command."""
    # Change to project root directory
    project_root = get_project_root()
    os.chdir(project_root)
    print(f"Working directory: {project_root}")
    
    # Gather experiment parameters
    exp = ask("experiment name", "exp1")
    mode = ask("mode (real/syn/mixed/finetune/all)", "all")
    
    #if user enters "-finetune", train all except finetune
    no_finetune = (mode == "-finetune")
    if no_finetune:
        mode = "all"  # Will be handled by --no-finetune flag
    
    col_count = ask("how many subcols?", TrainConfig.col_count)
    sidelobes_threshold = ask("sidelobes threshold (1 for no filtering)", TrainConfig.sidelobes_threshold)

    
    # Synthetic data params (only ask if needed)
    needs_synthetic = no_finetune or mode in ["syn", "mixed", "finetune", "all"]
    if needs_synthetic:
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
    extended_metrics = ask("Extended metrics (q95, error distribution, spatial error map)? (y/n)", "n").lower().startswith("y")
    create_plots = ask("Create plots after training? (y/n)", "n").lower().startswith("y")
    
    # Create log directory for SLURM output
    log_dir = Path("experiments") / exp / mode
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # Save description if provided
    if description:
        desc_path = Path("experiments") / exp / "description.txt"
        desc_path.parent.mkdir(parents=True, exist_ok=True)
        with open(desc_path, "w") as f:
            f.write(description)
        print(f"Description saved to {desc_path}")
    
    # Build training arguments
    train_args = [
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
    if no_finetune:
        train_args.append("--no-finetune")
    if extended_metrics:
        train_args.append("--extended-metrics")
    
    # Build SLURM command (always use sbatch for automatic logging)
    cmd = [
        "sbatch",
        "-J", exp,
        "-o", f"experiments/{exp}/{mode}/train_%j.log",
        "slurm/send_job.sh",
    ] + train_args
    
    print("\nSubmitting:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    
    # Extract job ID from sbatch output
    job_id = None
    for line in result.stdout.split('\n'):
        # Look for "Submitted batch job 12345"
        match = re.search(r'Submitted batch job (\d+)', line)
        if match:
            job_id = match.group(1)
            print(f"Job submitted: {job_id}")
            break
    
    if job_id:
        # Save plot configuration if requested
        if create_plots:
            plot_config = {
                "experiment": exp,
                "n_runs": int(n_runs),
                "extended_metrics": extended_metrics,
            }
            plot_config_path = project_root / "experiments" / exp / "plot_config.json"
            plot_config_path.parent.mkdir(parents=True, exist_ok=True)
            import json
            with open(plot_config_path, 'w') as f:
                json.dump(plot_config, f, indent=2)
            print(f"\nPlot configuration saved. Plots will be generated automatically after training completes.")
            print(f"  (Will plot all available modes, saved to experiments/{exp}/plots/comparison_performance.png)")


if __name__ == "__main__":
    main()
