#!/usr/bin/env python3
"""
Plot results from training experiments.

Usage: python -m stix_train.plot_results <experiment_name>
"""

import json
import os
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import matplotlib.pyplot as plt

from .config import EXPERIMENTS_DIR


def get_project_root() -> Path:
    """Get the project root directory (where setup.py is located)."""
    current_file = Path(__file__).resolve()
    project_root = current_file.parent.parent
    return project_root


def load_metrics(exp_dir: Path, mode: str, run_id: Optional[int] = None) -> Optional[Dict]:
    """Load extended metrics for a specific run."""
    if run_id is not None:
        metrics_file = exp_dir / mode / f"train_on_{mode}_run{run_id}_extended_metrics.json"
    else:
        metrics_file = exp_dir / mode / f"train_on_{mode}_extended_metrics.json"

    if not metrics_file.exists():
        return None

    with open(metrics_file, "r") as f:
        return json.load(f)


def load_error_data(exp_dir: Path, mode: str, run_id: Optional[int] = None) -> Optional[Dict]:
    """Load error distribution data for a specific run."""
    if run_id is not None:
        error_file = exp_dir / mode / f"train_on_{mode}_run{run_id}_errors.npz"
    else:
        error_file = exp_dir / mode / f"train_on_{mode}_errors.npz"

    if not error_file.exists():
        return None

    return dict(np.load(error_file))


def load_spatial_data(exp_dir: Path, mode: str, run_id: Optional[int] = None) -> Optional[Dict]:
    """Load spatial error map data for a specific run."""
    if run_id is not None:
        spatial_file = exp_dir / mode / f"train_on_{mode}_run{run_id}_spatial_errors.npz"
    else:
        spatial_file = exp_dir / mode / f"train_on_{mode}_spatial_errors.npz"

    if not spatial_file.exists():
        return None

    return dict(np.load(spatial_file))


def aggregate_metrics(all_metrics: List[Dict]) -> Dict:
    """Aggregate metrics from multiple runs."""
    aggregated: Dict = {}
    for key in all_metrics[0].keys():
        values = [m[key] for m in all_metrics if key in m]
        if values:
            aggregated[key] = {
                "mean": float(np.mean(values)),
                "std": float(np.std(values)),
                "min": float(np.min(values)),
                "max": float(np.max(values)),
            }
    return aggregated


def _label_inside_bar(ax, bar, text):
    h = bar.get_height()
    x = bar.get_x() + bar.get_width() / 2
    if h > 5:
        ax.text(x, h * 0.82, text, ha="center", va="top",
                fontsize=8, fontweight="bold", color="white")
    else:
        ax.text(x, h + 2, text, ha="center", va="bottom",
                fontsize=8, fontweight="bold")


def _pick_summary_mode(available_modes: List[str]) -> str:
    """
    Pick the mode to summarize (matches the reference intent: focus on Syn+Real Mixed if available).
    Preference: mixed > finetune > real > syn > first
    """
    for m in ["mixed", "finetune", "real", "syn"]:
        if m in available_modes:
            return m
    return available_modes[0]


def _get_mean(metrics: Dict, key: str, default: float = 0.0) -> float:
    v = metrics.get(key, default)
    if isinstance(v, dict) and "mean" in v:
        return float(v["mean"])
    try:
        return float(v)
    except Exception:
        return float(default)


def generate_metrics_from_model(exp_dir: Path, mode: str, run_id: Optional[int] = None) -> bool:
    """
    Generate extended metrics by loading model and predicting on test set.
    This is a fallback when metrics files don't exist.
    """
    try:
        from .training import compute_extended_metrics
        from .data import load_data, denormalize_locations
        from tensorflow import keras

        config_file = exp_dir / "config.json"
        if not config_file.exists():
            print(f"Error: config.json not found in {exp_dir}")
            return False

        with open(config_file, "r") as f:
            config_dict = json.load(f)

        model_dir = exp_dir / mode
        if run_id is not None:
            model_path = model_dir / f"model_run{run_id}.keras"
            metrics_name = f"train_on_{mode}_run{run_id}"
        else:
            model_path = model_dir / "model.keras"
            metrics_name = f"train_on_{mode}"

        if not model_path.exists():
            print(f"Error: Model not found at {model_path}")
            return False

        print(f"Loading model from {model_path}...")
        model = keras.models.load_model(model_path)

        print("Loading test data...")
        syn_data_path = config_dict.get("synthetic_data_path")
        if mode == "real":
            syn_data_path = None

        data = load_data(
            syn_data_path,
            config_dict.get("x_fov", 2200.0),
            config_dict.get("y_fov", 1800.0),
            config_dict.get("sidelobes_threshold", 0.84),
            config_dict.get("col_count", 12),
        )

        X_test, y_test = data["test_real"]

        print("Running predictions on test set...")
        y_pred_norm = model.predict(X_test, verbose=0)
        y_pred = denormalize_locations(y_pred_norm)
        y_true = denormalize_locations(y_test)

        print("Computing extended metrics...")
        compute_extended_metrics(y_true, y_pred, model_dir, metrics_name)

        print(f"✓ Successfully generated metrics for {mode} mode" +
              (f" (run {run_id})" if run_id is not None else ""))
        return True

    except Exception as e:
        print(f"Error generating metrics: {e}")
        import traceback
        traceback.print_exc()
        return False


def create_comparison_plot(
    exp_dir: Path,
    modes: List[str],
    n_runs: int,
    output_dir: Path,
    generate_if_missing: bool = False,
):
    """
    Create comparison plot with all available modes side-by-side.

    Styling is matched to the provided reference script as closely as possible.
    Only intentional difference: panel (d) uses P95 instead of P90.
    """
    # Load metrics for all modes
    all_modes_metrics: Dict[str, Dict] = {}
    all_modes_error_data: Dict[str, Optional[Dict]] = {}
    # Store raw metrics from all runs for boxplots
    all_modes_raw_metrics: Dict[str, List[Dict]] = {}

    for mode in modes:
        if n_runs > 1:
            all_metrics: List[Dict] = []
            missing_runs: List[int] = []
            for run_id in range(n_runs):
                metrics = load_metrics(exp_dir, mode, run_id)
                if metrics:
                    all_metrics.append(metrics)
                else:
                    missing_runs.append(run_id)

            if missing_runs and generate_if_missing:
                print(f"\nMissing metrics for {mode} mode, runs: {missing_runs}")
                print("Generating metrics from models...")
                for run_id in missing_runs:
                    if generate_metrics_from_model(exp_dir, mode, run_id):
                        metrics = load_metrics(exp_dir, mode, run_id)
                        if metrics:
                            all_metrics.append(metrics)

            if all_metrics:
                all_modes_metrics[mode] = aggregate_metrics(all_metrics)
                all_modes_raw_metrics[mode] = all_metrics  # Store raw metrics for boxplots
                all_modes_error_data[mode] = load_error_data(exp_dir, mode, 0)
                if not all_modes_error_data[mode] and generate_if_missing:
                    if generate_metrics_from_model(exp_dir, mode, 0):
                        all_modes_error_data[mode] = load_error_data(exp_dir, mode, 0)
        else:
            metrics = load_metrics(exp_dir, mode, None)
            if not metrics and generate_if_missing:
                print(f"\nMetrics not found for {mode} mode. Generating from model...")
                if generate_metrics_from_model(exp_dir, mode, None):
                    metrics = load_metrics(exp_dir, mode, None)

            if metrics:
                all_modes_metrics[mode] = metrics
                all_modes_raw_metrics[mode] = [metrics]  # Single run as list for consistency
                all_modes_error_data[mode] = load_error_data(exp_dir, mode, None)
                if not all_modes_error_data[mode] and generate_if_missing:
                    if generate_metrics_from_model(exp_dir, mode, None):
                        all_modes_error_data[mode] = load_error_data(exp_dir, mode, None)

    if not all_modes_metrics:
        print("Error: No metrics found for any mode")
        if not generate_if_missing:
            print("\nTip: Run with --generate-metrics to generate metrics from saved models")
        return

    mode_labels = {
        "syn": "Synthetic\nOnly",
        "mixed": "Syn+Real\nMixed",
        "real": "Real Only\n(Baseline)",
        "finetune": "Syn+Real\nFinetune",
    }

    model_colors = {
        "syn": ["#1e8449", "#58d68d"],
        "mixed": ["#1a5276", "#5dade2"],
        "real": ["#922b21", "#ec7063"],
        "finetune": ["#7d3c98", "#bb8fce"],
    }

    main_colors = {
        "syn": "#27ae60",
        "mixed": "#3498db",
        "real": "#e74c3c",
        "finetune": "#9b59b6",
    }

    # Only modes that actually have metrics
    available_modes = [m for m in modes if m in all_modes_metrics]
    if not available_modes:
        print("Error: No modes with metrics available")
        return

    fig, axes = plt.subplots(2, 3, figsize=(16, 11))
    fig.suptitle(
        f"Model Comparison: {n_runs} Training Runs (Mean ± Std)" if n_runs > 1 else "Model Comparison",
        fontsize=16,
        fontweight="bold",
    )

    x = np.arange(len(available_modes))
    w = 0.35

    summary_mode = _pick_summary_mode(available_modes)

    ax = axes[0, 0]
    for i, m in enumerate(available_modes):
        metrics = all_modes_metrics[m]
        mae_x = _get_mean(metrics, "mae_x")
        mae_y = _get_mean(metrics, "mae_y")

        bx = ax.bar(i - w/2, mae_x, w, color=model_colors[m][0], edgecolor='white', linewidth=1.5)
        by = ax.bar(i + w/2, mae_y, w, color=model_colors[m][1], edgecolor='white', linewidth=1.5)
        _label_inside_bar(ax, bx[0], f"{mae_x:.1f}")
        _label_inside_bar(ax, by[0], f"{mae_y:.1f}")
        
        # Add boxplots if we have multiple runs
        if n_runs > 1 and m in all_modes_raw_metrics:
            raw_metrics = all_modes_raw_metrics[m]
            mae_x_values = [rm.get("mae_x", 0) for rm in raw_metrics if "mae_x" in rm]
            mae_y_values = [rm.get("mae_y", 0) for rm in raw_metrics if "mae_y" in rm]
            
            if len(mae_x_values) > 1 and len(mae_y_values) > 1:
                bp_x = ax.boxplot(
                    [mae_x_values],
                    positions=[i - w/2],
                    widths=w * 0.6,
                    patch_artist=False,
                    showfliers=True,
                    flierprops=dict(marker='o', markersize=3, alpha=0.6, markeredgecolor='black', markerfacecolor='black'),
                    boxprops=dict(color='black', linewidth=1.5),
                    medianprops=dict(color='black', linewidth=2),
                    whiskerprops=dict(color='black', linewidth=1.5),
                    capprops=dict(color='black', linewidth=1.5),
                )
                
                bp_y = ax.boxplot(
                    [mae_y_values],
                    positions=[i + w/2],
                    widths=w * 0.6,
                    patch_artist=False,
                    showfliers=True,
                    flierprops=dict(marker='o', markersize=3, alpha=0.6, markeredgecolor='black', markerfacecolor='black'),
                    boxprops=dict(color='black', linewidth=1.5),
                    medianprops=dict(color='black', linewidth=2),
                    whiskerprops=dict(color='black', linewidth=1.5),
                    capprops=dict(color='black', linewidth=1.5),
                )

    ax.set_title("(a) Mean Absolute Error")
    ax.set_ylabel("MAE (arcsec)")
    ax.set_xticks(x)
    ax.set_xticklabels([mode_labels.get(m, m) for m in available_modes])
    ax.grid(axis="y", alpha=0.3)

    ax = axes[0, 1]
    for i, m in enumerate(available_modes):
        metrics = all_modes_metrics[m]
        med_x = _get_mean(metrics, "median_ae_x")
        med_y = _get_mean(metrics, "median_ae_y")

        bx = ax.bar(i - w/2, med_x, w, color=model_colors[m][0], edgecolor='white', linewidth=1.5)
        by = ax.bar(i + w/2, med_y, w, color=model_colors[m][1], edgecolor='white', linewidth=1.5)
        _label_inside_bar(ax, bx[0], f"{med_x:.1f}")
        _label_inside_bar(ax, by[0], f"{med_y:.1f}")
        
        # Add boxplots if we have multiple runs
        if n_runs > 1 and m in all_modes_raw_metrics:
            raw_metrics = all_modes_raw_metrics[m]
            med_x_values = [rm.get("median_ae_x", 0) for rm in raw_metrics if "median_ae_x" in rm]
            med_y_values = [rm.get("median_ae_y", 0) for rm in raw_metrics if "median_ae_y" in rm]
            
            if len(med_x_values) > 1 and len(med_y_values) > 1:
                bp_x = ax.boxplot(
                    [med_x_values],
                    positions=[i - w/2],
                    widths=w * 0.6,
                    patch_artist=False,
                    showfliers=True,
                    flierprops=dict(marker='o', markersize=3, alpha=0.6, markeredgecolor='black', markerfacecolor='black'),
                    boxprops=dict(color='black', linewidth=1.5),
                    medianprops=dict(color='black', linewidth=2),
                    whiskerprops=dict(color='black', linewidth=1.5),
                    capprops=dict(color='black', linewidth=1.5),
                )
                
                bp_y = ax.boxplot(
                    [med_y_values],
                    positions=[i + w/2],
                    widths=w * 0.6,
                    patch_artist=False,
                    showfliers=True,
                    flierprops=dict(marker='o', markersize=3, alpha=0.6, markeredgecolor='black', markerfacecolor='black'),
                    boxprops=dict(color='black', linewidth=1.5),
                    medianprops=dict(color='black', linewidth=2),
                    whiskerprops=dict(color='black', linewidth=1.5),
                    capprops=dict(color='black', linewidth=1.5),
                )

    ax.set_title("(b) Median Absolute Error")
    ax.set_ylabel("Median AE (arcsec)")
    ax.set_xticks(x)
    ax.set_xticklabels([mode_labels.get(m, m) for m in available_modes])
    ax.grid(axis="y", alpha=0.3)

    ax = axes[0, 2]
    for i, m in enumerate(available_modes):
        metrics = all_modes_metrics[m]
        euc = _get_mean(metrics, "euclidean_mae")

        b = ax.bar(i, euc, 0.6, color=main_colors[m], edgecolor='white', linewidth=1.5)
        _label_inside_bar(ax, b[0], f"{euc:.1f}")
        
        # Add boxplots if we have multiple runs
        if n_runs > 1 and m in all_modes_raw_metrics:
            raw_metrics = all_modes_raw_metrics[m]
            euc_values = [rm.get("euclidean_mae", 0) for rm in raw_metrics if "euclidean_mae" in rm]
            
            if len(euc_values) > 1:
                bp = ax.boxplot(
                    [euc_values],
                    positions=[i],
                    widths=0.6 * 0.6,
                    patch_artist=False,
                    showfliers=True,
                    flierprops=dict(marker='o', markersize=3, alpha=0.6, markeredgecolor='black', markerfacecolor='black'),
                    boxprops=dict(color='black', linewidth=1.5),
                    medianprops=dict(color='black', linewidth=2),
                    whiskerprops=dict(color='black', linewidth=1.5),
                    capprops=dict(color='black', linewidth=1.5),
                )

    ax.set_title("(c) Euclidean Distance Error")
    ax.set_ylabel("Euclidean Error (arcsec)")
    ax.set_xticks(x)
    ax.set_xticklabels([mode_labels.get(m, m) for m in available_modes])
    ax.grid(axis="y", alpha=0.3)

    ax = axes[1, 0]
    for i, m in enumerate(available_modes):
        metrics = all_modes_metrics[m]
        p95_x = _get_mean(metrics, "p95_x")
        p95_y = _get_mean(metrics, "p95_y")

        bx = ax.bar(i - w/2, p95_x, w, color=model_colors[m][0], edgecolor='white', linewidth=1.5)
        by = ax.bar(i + w/2, p95_y, w, color=model_colors[m][1], edgecolor='white', linewidth=1.5)
        _label_inside_bar(ax, bx[0], f"{p95_x:.1f}")
        _label_inside_bar(ax, by[0], f"{p95_y:.1f}")
        
        # Add boxplots if we have multiple runs
        if n_runs > 1 and m in all_modes_raw_metrics:
            raw_metrics = all_modes_raw_metrics[m]
            p95_x_values = [rm.get("p95_x", 0) for rm in raw_metrics if "p95_x" in rm]
            p95_y_values = [rm.get("p95_y", 0) for rm in raw_metrics if "p95_y" in rm]
            
            if len(p95_x_values) > 1 and len(p95_y_values) > 1:
                bp_x = ax.boxplot(
                    [p95_x_values],
                    positions=[i - w/2],
                    widths=w * 0.6,
                    patch_artist=False,
                    showfliers=True,
                    flierprops=dict(marker='o', markersize=3, alpha=0.6, markeredgecolor='black', markerfacecolor='black'),
                    boxprops=dict(color='black', linewidth=1.5),
                    medianprops=dict(color='black', linewidth=2),
                    whiskerprops=dict(color='black', linewidth=1.5),
                    capprops=dict(color='black', linewidth=1.5),
                )
                
                bp_y = ax.boxplot(
                    [p95_y_values],
                    positions=[i + w/2],
                    widths=w * 0.6,
                    patch_artist=False,
                    showfliers=True,
                    flierprops=dict(marker='o', markersize=3, alpha=0.6, markeredgecolor='black', markerfacecolor='black'),
                    boxprops=dict(color='black', linewidth=1.5),
                    medianprops=dict(color='black', linewidth=2),
                    whiskerprops=dict(color='black', linewidth=1.5),
                    capprops=dict(color='black', linewidth=1.5),
                )

    ax.set_title("(d) 95th Percentile Error")
    ax.set_ylabel("P95 Error (arcsec)")
    ax.set_xticks(x)
    ax.set_xticklabels([mode_labels.get(m, m) for m in available_modes])
    ax.grid(axis="y", alpha=0.3)

    ax = axes[1, 1]
    bins = np.linspace(0, 500, 60)
    for m in available_modes:
        ed = all_modes_error_data.get(m) or {}
        errors = ed.get("euclidean_errors", None)
        if errors is None:
            continue
        # Filled histogram with transparency
        ax.hist(
            errors,
            bins=bins,
            alpha=0.4,
            histtype="stepfilled",
            color=main_colors[m],
            edgecolor=main_colors[m],
            linewidth=1.5,
            label=mode_labels.get(m, m),
        )
        # Add outline for better definition
        ax.hist(
            errors,
            bins=bins,
            alpha=1.0,
            histtype="step",
            color=main_colors[m],
            edgecolor=main_colors[m],
            linewidth=1.5,
            label=None,
        )

    ax.set_yscale("log")
    ax.set_title("(e) Error Distribution")
    ax.set_xlabel("Euclidean Error (arcsec)")
    ax.set_ylabel("Count (log)")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.95, edgecolor="gray")
    ax.grid(alpha=0.3)

    ax = axes[1, 2]
    spatial_mode = "mixed" if "mixed" in available_modes else summary_mode
    spatial_data = load_spatial_data(exp_dir, spatial_mode, 0 if n_runs > 1 else None)

    if not spatial_data and generate_if_missing:
        if generate_metrics_from_model(exp_dir, spatial_mode, 0 if n_runs > 1 else None):
            spatial_data = load_spatial_data(exp_dir, spatial_mode, 0 if n_runs > 1 else None)

    if spatial_data:
        true_x = spatial_data.get("true_x")
        true_y = spatial_data.get("true_y")
        euclidean_error = spatial_data.get("euclidean_error")

        if true_x is None or true_y is None or euclidean_error is None:
            ax.text(0.5, 0.5, "Spatial error map data incomplete", ha="center", va="center",
                    transform=ax.transAxes, fontsize=12)
            ax.set_title("(f) Spatial Error Map")
        else:
            ax.scatter(true_x, true_y, s=4, alpha=0.05, color="black", label="True")

            hb = ax.hexbin(
                true_x, true_y,
                C=euclidean_error,
                reduce_C_function=np.mean,
                gridsize=55,
                mincnt=10,
                vmin=0,
                vmax=100,  # Set max to 100 to better show errors around 20
            )
            cb = fig.colorbar(hb, ax=ax)
            cb.set_label("Mean Euclidean Error (arcsec)")

            if spatial_mode == "mixed":
                ax.set_title("(f) Spatial Error Map (Syn+Real Mixed)")
            else:
                ax.set_title(f"(f) Spatial Error Map ({mode_labels.get(spatial_mode, spatial_mode)})")

            ax.set_xlabel("X (arcsec)")
            ax.set_ylabel("Y (arcsec)")
            ax.set_aspect("equal")
            ax.grid(alpha=0.3)
    else:
        ax.text(0.5, 0.5, "Spatial error map data not available", ha="center", va="center",
                transform=ax.transAxes, fontsize=12)
        ax.set_title("(f) Spatial Error Map")

    plt.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "comparison_performance.png"
    plt.savefig(output_path, dpi=300)
    print(f"Saved: {output_path}")
    plt.close()


def find_available_modes(exp_dir: Path) -> List[str]:
    """Find available training modes in experiment directory."""
    modes: List[str] = []
    for mode in ["real", "syn", "mixed", "finetune"]:
        mode_dir = exp_dir / mode
        if mode_dir.exists():
            if any(mode_dir.glob("*_extended_metrics.json")):
                modes.append(mode)
    return modes


def main():
    parser = argparse.ArgumentParser(description="Plot training results")
    parser.add_argument("experiment", type=str, nargs="?", help="Experiment name")
    parser.add_argument("--mode", type=str, default=None, help="Training mode to plot (ignored - always plots all)")
    parser.add_argument(
        "--generate-metrics",
        action="store_true",
        help="Generate metrics from saved models if they don't exist (rare fallback option)",
    )

    args = parser.parse_args()

    # Change to project root directory
    project_root = get_project_root()
    os.chdir(project_root)

    if not args.experiment:
        experiments = [d.name for d in EXPERIMENTS_DIR.iterdir() if d.is_dir()]
        if not experiments:
            print("No experiments found")
            sys.exit(1)

        print("Available experiments:")
        for i, exp in enumerate(experiments, 1):
            print(f"  [{i}] {exp}")

        try:
            choice = input("\nSelect experiment (number): ").strip()
            idx = int(choice) - 1
            if idx < 0 or idx >= len(experiments):
                print("Invalid selection")
                sys.exit(1)
            experiment_name = experiments[idx]
        except (ValueError, EOFError, KeyboardInterrupt):
            print("Cancelled")
            sys.exit(1)
    else:
        experiment_name = args.experiment

    exp_dir = EXPERIMENTS_DIR / experiment_name
    if not exp_dir.exists():
        print(f"Error: Experiment directory not found: {exp_dir}")
        sys.exit(1)

    available_modes = find_available_modes(exp_dir)
    if not available_modes:
        print(f"Error: No training results found in {exp_dir}")
        sys.exit(1)

    # Detect n_runs
    n_runs = 1
    for run_id in range(10):
        p = exp_dir / available_modes[0] / f"train_on_{available_modes[0]}_run{run_id}_extended_metrics.json"
        if p.exists():
            n_runs = max(n_runs, run_id + 1)

    print(f"\nPlotting experiment: {experiment_name}")
    print(f"Modes: {', '.join(available_modes)}")
    print(f"Runs: {n_runs}")

    output_dir = exp_dir / "plots"
    create_comparison_plot(
        exp_dir,
        available_modes,
        n_runs,
        output_dir,
        generate_if_missing=args.generate_metrics,
    )

    print(f"\nDone! Plot saved to: {output_dir}/comparison_performance.png")


if __name__ == "__main__":
    main()
