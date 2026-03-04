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
from typing import List, Optional

from .config import TrainConfig, EXPERIMENTS_DIR
from .data import load_data, ensure_synthetic_data
from .training import train_single_mode, save_experiment_config


def parse_args() -> TrainConfig:
    """Parse command-line arguments and return TrainConfig."""
    parser = argparse.ArgumentParser(description="Train STIX localization models")
    
    parser.add_argument("--experiment", type=str, required=True, help="Experiment name")
    parser.add_argument("--n-samples", type=int, default=1_000_000, help="Number of synthetic samples")
    parser.add_argument("--fov-big", type=int, default=7200, help="FOV_big for synthetic data")
    parser.add_argument("--x-fov", type=float, default=99999.0, help="X FOV bound (99999 = no filter)")
    parser.add_argument("--y-fov", type=float, default=99999.0, help="Y FOV bound (99999 = no filter)")
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
        fov_big=args.fov_big,
        x_fov=args.x_fov,
        y_fov=args.y_fov,
        learning_rate=args.lr,
        epochs=args.epochs,
        batch_size=args.batch,
        patience=args.patience,
        n_runs=args.n_runs,
        col_count=args.col_count,
        sidelobes_threshold=args.sidelobes_threshold,
        hidden_dims=hidden_dims_list
    )


def validate_paths(config: TrainConfig, modes_to_train: list) -> Optional[str]:
    """
    Validate that all required data files exist before starting heavy work.
    Returns the synthetic data path if needed, None otherwise.
    Raises FileNotFoundError early if anything is missing.
    """
    from .config import REAL_DATA_PATH, SYNTHETIC_DATA_DIR

    if not REAL_DATA_PATH.exists():
        raise FileNotFoundError(
            f"Real data file not found: {REAL_DATA_PATH}\n"
            f"Expected at: {REAL_DATA_PATH}"
        )

    needs_synthetic = any(m in ["syn", "mixed", "finetune"] for m in modes_to_train)
    if needs_synthetic:
        return ensure_synthetic_data(config.n_samples, config.fov_big)
    return None


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
    
    # Validate all paths early (before importing TensorFlow)
    syn_data_path = validate_paths(config, modes_to_train)
    
    # Save config immediately so it persists even if training crashes
    save_experiment_config(config, syn_data_path)
    
    # Load data
    print("[DEBUG] Loading data...", flush=True)
    data = load_data(syn_data_path, config.x_fov, config.y_fov, config.sidelobes_threshold, config.col_count)
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
        print("FINAL COMPARISON (on real test set)")
        print("=" * 60)
        print(f"{'Mode':10s}  {'MAE_X':>8s}  {'MAE_Y':>8s}  {'MAE':>8s}  {'Eucl':>8s}  {'Q95':>8s}")
        print("-" * 60)
        for name, m in results.items():
            print(f"{name:10s}  {m['mae_x']:8.2f}  {m['mae_y']:8.2f}  "
                  f"{m['mae_mean']:8.2f}  {m.get('euclidean_mean', 0):8.2f}  "
                  f"{m.get('euclidean_q95', 0):8.2f}")
        
        best_name = min(results.keys(), key=lambda k: results[k].get("euclidean_mean", float("inf")))
        print(f"\nBest method: {best_name} "
              f"(Eucl mean: {results[best_name].get('euclidean_mean', 0):.2f}, "
              f"MAE mean: {results[best_name]['mae_mean']:.2f})")
    


if __name__ == "__main__":
    main()
