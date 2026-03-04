#!/usr/bin/env python3
"""
STIX Model Evaluation CLI.

Entry point for the stix-test command.
Lists experiments, lets the user pick one, and runs evaluation
to generate extended metrics.
"""

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

import json
import sys
from pathlib import Path
from typing import List, Optional, Tuple

from .config import EXPERIMENTS_DIR


MODES_ORDER = ["real", "syn", "mixed", "finetune"]


def find_models(exp_dir: Path) -> List[Tuple[str, Optional[int], Path]]:
    """
    Discover trained models in an experiment directory.
    
    Returns:
        List of (mode, run_id_or_None, model_path) tuples.
    """
    models = []
    for mode in MODES_ORDER:
        mode_dir = exp_dir / mode
        if not mode_dir.exists():
            continue

        single = mode_dir / "model.keras"
        if single.exists():
            models.append((mode, None, single))
            continue

        run_id = 0
        while True:
            p = mode_dir / f"model_run{run_id}.keras"
            if not p.exists():
                break
            models.append((mode, run_id, p))
            run_id += 1

    return models


def metrics_name_for(mode: str, run_id: Optional[int]) -> str:
    """Build the metrics filename prefix matching plot_results.py conventions."""
    base = "finetune" if mode == "finetune" else f"train_on_{mode}"
    if run_id is not None:
        return f"{base}_run{run_id}"
    return base


def _label(mode: str, run_id: Optional[int]) -> str:
    if run_id is not None:
        return f"{mode} (run {run_id})"
    return mode


def print_detailed_report(all_results: List[Tuple[str, Optional[int], dict]]):
    """Print a comprehensive per-model report with all metrics."""
    W = 72

    for mode, run_id, m in all_results:
        label = _label(mode, run_id)
        print(f"\n{'=' * W}")
        print(f"  {label.upper()}")
        print(f"{'=' * W}")
        print(f"  Test samples: {m.get('n_test_samples', '?')}")

        # --- Absolute Error ---
        print(f"\n  {'Absolute Error':30s} {'X':>12s} {'Y':>12s} {'Mean':>12s}")
        print(f"  {'-'*66}")
        print(f"  {'Mean  (MAE)':30s} {m['mae_x']:12.2f} {m['mae_y']:12.2f} {m['mae_mean']:12.2f}")
        print(f"  {'Median (MedAE)':30s} {m.get('median_ae_x', 0):12.2f} {m.get('median_ae_y', 0):12.2f} {m.get('median_ae_mean', 0):12.2f}")
        print(f"  {'MSE':30s} {m.get('mse_x', 0):12.1f} {m.get('mse_y', 0):12.1f} {(m.get('mse_x', 0) + m.get('mse_y', 0)) / 2:12.1f}")
        print(f"  {'RMSE':30s} {m.get('rmse_x', 0):12.2f} {m.get('rmse_y', 0):12.2f} {(m.get('rmse_x', 0) + m.get('rmse_y', 0)) / 2:12.2f}")

        # --- Euclidean Distance ---
        print(f"\n  {'Euclidean Distance':30s} {'Value':>12s}")
        print(f"  {'-'*42}")
        print(f"  {'Mean':30s} {m.get('euclidean_mean', 0):12.2f}")
        print(f"  {'Median':30s} {m.get('euclidean_median', 0):12.2f}")
        print(f"  {'Std':30s} {m.get('euclidean_std', 0):12.2f}")

        # --- Percentiles ---
        print(f"\n  {'Percentiles':30s} {'X':>12s} {'Y':>12s} {'Eucl':>12s}")
        print(f"  {'-'*66}")
        print(f"  {'Q90  (90th)':30s} {m.get('q90_x', 0):12.2f} {m.get('q90_y', 0):12.2f} {m.get('euclidean_q90', 0):12.2f}")
        print(f"  {'Q95  (95th)':30s} {m.get('q95_x', 0):12.2f} {m.get('q95_y', 0):12.2f} {m.get('euclidean_q95', 0):12.2f}")
        print(f"  {'Q99  (99th)':30s} {m.get('q99_x', 0):12.2f} {m.get('q99_y', 0):12.2f} {m.get('euclidean_q99', 0):12.2f}")
        print(f"  {'Max':30s} {m.get('max_x', 0):12.2f} {m.get('max_y', 0):12.2f} {m.get('euclidean_max', 0):12.2f}")

        # --- Quick interpretation ---
        median_eucl = m.get('euclidean_median', 0)
        q95_eucl = m.get('euclidean_q95', 0)
        ratio = q95_eucl / median_eucl if median_eucl > 0 else 0
        print(f"\n  Q95/Median ratio: {ratio:.1f}x", end="")
        if ratio > 5:
            print("  (heavy tail -- many large outliers)")
        elif ratio > 3:
            print("  (moderate outliers)")
        else:
            print("  (tight distribution)")

    # --- Comparison table if multiple models ---
    if len(all_results) > 1:
        print(f"\n\n{'=' * W}")
        print("  COMPARISON SUMMARY")
        print(f"{'=' * W}")
        header = (f"  {'Mode':16s} {'MAE_X':>7s} {'MAE_Y':>7s} {'MAE':>7s} "
                  f"{'MedAE':>7s} {'Eucl':>7s} {'Q95':>7s} {'Max':>8s}")
        print(header)
        print(f"  {'-' * (W - 2)}")
        for mode, run_id, m in all_results:
            label = _label(mode, run_id)
            print(f"  {label:16s} {m['mae_x']:7.1f} {m['mae_y']:7.1f} {m['mae_mean']:7.1f} "
                  f"{m.get('median_ae_mean', 0):7.1f} {m.get('euclidean_mean', 0):7.1f} "
                  f"{m.get('euclidean_q95', 0):7.1f} {m.get('euclidean_max', 0):8.1f}")
        print(f"  {'-' * (W - 2)}")

        best = min(all_results, key=lambda x: x[2].get('euclidean_mean', float('inf')))
        print(f"\n  Best (Eucl mean): {_label(best[0], best[1])} = {best[2].get('euclidean_mean', 0):.2f}")
        print(f"{'=' * W}")


def evaluate_experiment(exp_dir: Path, force: bool = False):
    """Load models for an experiment and compute extended metrics."""
    from tensorflow import keras
    from .data import load_data, denormalize_locations
    from .training import compute_extended_metrics

    config_path = exp_dir / "config.json"
    if not config_path.exists():
        print(f"Error: config.json not found in {exp_dir}")
        sys.exit(1)

    with open(config_path, "r") as f:
        config = json.load(f)

    models = find_models(exp_dir)
    if not models:
        print(f"No trained models found in {exp_dir}")
        sys.exit(1)

    print(f"\nExperiment: {exp_dir.name}")
    print(f"Config: col_count={config.get('col_count')}, "
          f"sidelobes={config.get('sidelobes_threshold')}, "
          f"x_fov={config.get('x_fov')}, y_fov={config.get('y_fov')}")
    print(f"Found {len(models)} model(s): "
          + ", ".join(f"{m}{'_run'+str(r) if r is not None else ''}" for m, r, _ in models))

    # Check which models already have metrics
    to_evaluate = []
    already_done = []
    for mode, run_id, model_path in models:
        name = metrics_name_for(mode, run_id)
        metrics_file = model_path.parent / f"{name}_extended_metrics.json"
        if metrics_file.exists() and not force:
            with open(metrics_file, "r") as f:
                already_done.append((mode, run_id, json.load(f)))
        else:
            to_evaluate.append((mode, run_id, model_path))

    if already_done and not to_evaluate:
        print(f"\nAll {len(already_done)} model(s) already evaluated.")
        print("Use --force to re-evaluate.")
        print_detailed_report(already_done)
        return

    if already_done:
        print(f"\n{len(already_done)} model(s) already evaluated, "
              f"{len(to_evaluate)} remaining.")

    # Load test data once (only need real test data for evaluation)
    print("\nLoading test data...")
    data = load_data(
        syn_data_path=None,
        x_fov=config.get("x_fov", 99999.0),
        y_fov=config.get("y_fov", 99999.0),
        sidelobes_threshold=config.get("sidelobes_threshold", 0.84),
        col_count=config.get("col_count", 12),
    )
    X_test, y_test = data["test_real"]
    print(f"Test set: {X_test.shape[0]} samples, {X_test.shape[1]} features")

    new_results = []
    for mode, run_id, model_path in to_evaluate:
        name = metrics_name_for(mode, run_id)
        run_label = f" (run {run_id})" if run_id is not None else ""
        print(f"\nEvaluating {mode}{run_label}...")
        print(f"  Loading model: {model_path}")

        model = keras.models.load_model(model_path)
        y_pred = denormalize_locations(model.predict(X_test, verbose=0))
        y_true = denormalize_locations(y_test)

        metrics = compute_extended_metrics(y_true, y_pred, model_path.parent, name)
        new_results.append((mode, run_id, metrics))

    all_results = already_done + new_results
    all_results.sort(key=lambda x: (MODES_ORDER.index(x[0]) if x[0] in MODES_ORDER else 99, x[1] or 0))
    print_detailed_report(all_results)


def main():
    """Main entry point for stix-test command."""
    import argparse

    parser = argparse.ArgumentParser(description="Evaluate trained STIX models")
    parser.add_argument("--experiment", type=str, default=None,
                        help="Experiment name (skip interactive picker)")
    parser.add_argument("--force", action="store_true",
                        help="Re-evaluate even if metrics already exist")
    args = parser.parse_args()

    if args.experiment:
        exp_dir = EXPERIMENTS_DIR / args.experiment
        if not exp_dir.exists():
            print(f"Error: Experiment not found: {exp_dir}")
            sys.exit(1)
        evaluate_experiment(exp_dir, force=args.force)
        return

    # Interactive experiment picker
    if not EXPERIMENTS_DIR.exists():
        print(f"No experiments directory found at {EXPERIMENTS_DIR}")
        sys.exit(1)

    experiments = sorted(
        [d for d in EXPERIMENTS_DIR.iterdir() if d.is_dir()],
        key=lambda d: d.name,
    )
    if not experiments:
        print("No experiments found.")
        sys.exit(1)

    print("Available experiments:\n")
    for i, exp in enumerate(experiments, 1):
        models = find_models(exp)
        modes = sorted(set(m for m, _, _ in models))
        n_models = len(models)
        config_path = exp / "config.json"
        config_info = ""
        if config_path.exists():
            with open(config_path, "r") as f:
                c = json.load(f)
            config_info = f"  cols={c.get('col_count', '?')}, sidelobes={c.get('sidelobes_threshold', '?')}"
        status = f"{n_models} model(s) [{', '.join(modes)}]" if models else "no models"
        print(f"  [{i:2d}] {exp.name:40s}  {status}{config_info}")

    try:
        choice = input("\nSelect experiment (number): ").strip()
        idx = int(choice) - 1
        if idx < 0 or idx >= len(experiments):
            print("Invalid selection.")
            sys.exit(1)
    except (ValueError, EOFError, KeyboardInterrupt):
        print("\nCancelled.")
        sys.exit(1)

    evaluate_experiment(experiments[idx], force=args.force)


if __name__ == "__main__":
    main()
