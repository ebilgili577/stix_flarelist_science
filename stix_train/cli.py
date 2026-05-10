#!/usr/bin/env python3
"""
STIX Localization Model Training CLI.

Entry point for the stix-train command.
"""

# Suppress TensorFlow/XLA verbose logging before any imports
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'  
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

import argparse
import numpy as np

from .config import TrainConfig, EXPERIMENTS_DIR
from .data import load_data, ensure_synthetic_data
from .training import train_single_mode, save_experiment_config


def parse_args() -> TrainConfig:
    """Parse command-line arguments and return TrainConfig."""
    parser = argparse.ArgumentParser(description="Train STIX localization models")
    
    parser.add_argument("--experiment", type=str, required=True, help="Experiment name")
    parser.add_argument("--n-samples", type=int, default=1_000_000, help="Number of synthetic samples")
    parser.add_argument("--x-min", type=float, default=-2295.0, help="X min")
    parser.add_argument("--x-max", type=float, default=2295.0, help="X max")
    parser.add_argument("--y-min", type=float, default=-1878.0, help="y min ")
    parser.add_argument("--y-max", type=float, default=1878.0, help="Y max")
    parser.add_argument("--lr", type=float, default=1e-3, help="Learning rate")
    parser.add_argument("--epochs", type=int, default=1000, help="Number of epochs")
    parser.add_argument("--batch", type=int, default=512, help="Batch size")
    parser.add_argument("--patience", type=int, default=50, help="Early stopping patience")
    parser.add_argument(
        "--mode",
        type=str,
        default="all",
        choices=["real", "syn", "mixed", "finetune", "all"],
        help="Training mode"
    )
    parser.add_argument("--n-runs", type=int, default=1, help="Number of runs for statistics")
    parser.add_argument("--col-count", type=int, default=12, help="How many subcols to use")
    parser.add_argument("--sidelobes-threshold", type=float, default=0.84, help="Sidelobes threshold")
    parser.add_argument(
        "--hidden-dims",
        type=str,
        default="1024,512,256,128,64",
        help="Comma-separated list of hidden layer sizes (e.g., '1024,512,256,128,64')"
    )
    
    args = parser.parse_args()
    
    # Parse hidden_dims string to list
    hidden_dims_list = [int(x.strip()) for x in args.hidden_dims.split(",")]
    
    return TrainConfig(
        experiment=args.experiment,
        mode=args.mode,
        n_samples=args.n_samples,
        x_min=args.x_min,
        x_max=args.x_max,
        y_min=args.y_min,
        y_max=args.y_max,
        learning_rate=args.lr,
        epochs=args.epochs,
        batch_size=args.batch,
        patience=args.patience,
        n_runs=args.n_runs,
        col_count=args.col_count,
        sidelobes_threshold=args.sidelobes_threshold,
        hidden_dims=hidden_dims_list,
    )


def main():
    """Main entry point for stix-train command."""
    print("[DEBUG] Starting stix-train", flush=True)
    
    # Parse arguments
    config = parse_args()
    print(f"[DEBUG] Parsed config: {config}", flush=True)
    
    # Determine which modes to train
    if config.mode == "all":
        modes_to_train = ["real", "syn", "mixed", "finetune"]
    else:
        modes_to_train = [config.mode]
    
    # Load synthetic data if needed
    needs_synthetic = any(m in ["syn", "mixed", "finetune"] for m in modes_to_train)
    syn_data_path = None
    if needs_synthetic:
        syn_data_path = ensure_synthetic_data(n_samples=config.n_samples, min_x=config.x_min, max_x=config.x_max, min_y=config.y_min, max_y=config.y_max)

    # Load data
    print("[DEBUG] Loading data...", flush=True)
    data = load_data(syn_data_path=syn_data_path, x_min=config.x_min,x_max=config.x_max, y_min=config.y_min, y_max=config.y_max, sidelobes_threshold=config.sidelobes_threshold, col_count=config.col_count)
    print("[DEBUG] Data loaded", flush=True)
    
    # Train models
    results = {}
    syn_model_path = None
    
    for mode in modes_to_train:
        pretrained_path = syn_model_path if mode == "finetune" else None
        
        mae = train_single_mode(mode, data, config, pretrained_path)
        results[mode] = mae
        
        # Track syn model path for finetune reuse
        if mode == "syn":
            syn_dir = EXPERIMENTS_DIR / config.experiment / "syn"
            if config.n_runs > 1:
                syn_model_path = str(syn_dir / "model_run0.keras")
            else:
                syn_model_path = str(syn_dir / "model.keras")
    
    # Print comparison if multiple modes
    if len(modes_to_train) > 1:
        print("\n" + "=" * 60)
        print("FINAL COMPARISON (MAE on real test set)")
        print("=" * 60)
        for name, mae in results.items():
            if isinstance(mae, np.ndarray) and len(mae.shape) == 1:
                print(f"{name:10s}: X={mae[0]:7.2f}, Y={mae[1]:7.2f}, Mean={np.mean(mae):7.2f}")
            else:
                print(f"{name:10s}: {mae}")
        
        if config.n_runs == 1:
            best_name = min(results.keys(), key=lambda k: np.mean(results[k]))
            print(f"\nBest method: {best_name} (mean MAE: {np.mean(results[best_name]):.2f})")
    
    # Save configuration
    save_experiment_config(config, syn_data_path)


if __name__ == "__main__":
    main()
