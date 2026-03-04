#!/usr/bin/env python3
"""
Plot results from training experiments - FIXED VERSION with proper boxplot scaling.

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
    if mode == "finetune":
        # Finetune uses different naming: finetune_run{id}_extended_metrics.json
        if run_id is not None:
            metrics_file = exp_dir / mode / f"finetune_run{run_id}_extended_metrics.json"
        else:
            metrics_file = exp_dir / mode / f"finetune_extended_metrics.json"
    else:
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
    if mode == "finetune":
        # Finetune uses different naming: finetune_run{id}_errors.npz
        if run_id is not None:
            error_file = exp_dir / mode / f"finetune_run{run_id}_errors.npz"
        else:
            error_file = exp_dir / mode / f"finetune_errors.npz"
    else:
        if run_id is not None:
            error_file = exp_dir / mode / f"train_on_{mode}_run{run_id}_errors.npz"
        else:
            error_file = exp_dir / mode / f"train_on_{mode}_errors.npz"

    if not error_file.exists():
        return None

    return dict(np.load(error_file))


def load_spatial_data(exp_dir: Path, mode: str, run_id: Optional[int] = None) -> Optional[Dict]:
    """Load spatial error map data for a specific run."""
    if mode == "finetune":
        # Finetune uses different naming: finetune_run{id}_spatial_errors.npz
        if run_id is not None:
            spatial_file = exp_dir / mode / f"finetune_run{run_id}_spatial_errors.npz"
        else:
            spatial_file = exp_dir / mode / f"finetune_spatial_errors.npz"
    else:
        if run_id is not None:
            spatial_file = exp_dir / mode / f"train_on_{mode}_run{run_id}_spatial_errors.npz"
        else:
            spatial_file = exp_dir / mode / f"train_on_{mode}_spatial_errors.npz"

    if not spatial_file.exists():
        return None

    return dict(np.load(spatial_file))


def aggregate_error_data(exp_dir: Path, mode: str, n_runs: int) -> Optional[Dict]:
    """Aggregate error distribution data from multiple runs."""
    all_errors = []
    
    for run_id in range(n_runs):
        data = load_error_data(exp_dir, mode, run_id)
        if data and "euclidean_errors" in data:
            all_errors.append(data["euclidean_errors"])
            
    if not all_errors:
        return None
        
    return {"euclidean_errors": np.concatenate(all_errors)}


def aggregate_spatial_data(exp_dir: Path, mode: str, n_runs: int) -> Optional[Dict]:
    """Aggregate spatial error map data from multiple runs."""
    all_true_x = []
    all_true_y = []
    all_pred_x = []
    all_pred_y = []
    all_errors = []
    
    for run_id in range(n_runs):
        data = load_spatial_data(exp_dir, mode, run_id)
        if data:
            if "true_x" in data and "true_y" in data and "euclidean_error" in data:
                all_true_x.append(data["true_x"])
                all_true_y.append(data["true_y"])
                all_errors.append(data["euclidean_error"])
                
                # Collect predictions if available (for MAE calculation)
                if "pred_x" in data and "pred_y" in data:
                    all_pred_x.append(data["pred_x"])
                    all_pred_y.append(data["pred_y"])
    
    if not all_true_x:
        return None
        
    result = {
        "true_x": np.concatenate(all_true_x),
        "true_y": np.concatenate(all_true_y),
        "euclidean_error": np.concatenate(all_errors)
    }
    
    if all_pred_x and len(all_pred_x) == len(all_true_x):
        result["pred_x"] = np.concatenate(all_pred_x)
        result["pred_y"] = np.concatenate(all_pred_y)
        
    return result


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


def _add_panel_label(ax, label: str, x: float = -0.02, y: float = 1.02) -> None:
    """Add panel label (a), (b), etc. outside panel, aligned with title."""
    ax.text(x, y, f"{label})", transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom", ha="right")


def _label_inside_bar(ax, bar, text):
    h = bar.get_height()
    x = bar.get_x() + bar.get_width() / 2
    ax.text(x, h * 0.5, text, ha="center", va="center",
            fontsize=10, fontweight="bold", color="white")


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
        if mode == "finetune":
            if run_id is not None:
                model_path = model_dir / f"model_run{run_id}.keras"
                metrics_name = f"finetune_run{run_id}"
            else:
                model_path = model_dir / "model.keras"
                metrics_name = "finetune"
        else:
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

        x_fov = config_dict.get("x_fov", 2200.0)
        y_fov = config_dict.get("y_fov", 1800.0)
        location_norm = config_dict.get("location_norm", "global")
        counts_norm = config_dict.get("counts_norm", "max")

        data = load_data(
            syn_data_path,
            x_fov,
            y_fov,
            config_dict.get("sidelobes_threshold", 0.84),
            config_dict.get("col_count", 12),
            counts_norm=counts_norm,
            location_norm=location_norm,
        )

        X_test, y_test = data["test_real"]

        print("Running predictions on test set...")
        y_pred_norm = model.predict(X_test, verbose=0)
        y_pred = denormalize_locations(y_pred_norm, method=location_norm, x_fov=x_fov, y_fov=y_fov)
        y_true = denormalize_locations(y_test, method=location_norm, x_fov=x_fov, y_fov=y_fov)

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
    specific_run: Optional[int] = None,
    figsize: Optional[tuple] = None,
):
    """
    Create comparison plots with all available modes side-by-side.
    Split into two figures for better readability:
    1. comparison_metrics.png: Panels (a)-(d) (Metrics)
    2. comparison_analysis.png: Panels (e)-(f) (Distributions & Maps)
    """
    # Load metrics for all modes
    all_modes_metrics: Dict[str, Dict] = {}
    all_modes_error_data: Dict[str, Optional[Dict]] = {}
    # Store raw metrics from all runs for boxplots
    all_modes_raw_metrics: Dict[str, List[Dict]] = {}

    for mode in modes:
        if specific_run is not None:
            # Plot only the specific run
            metrics = load_metrics(exp_dir, mode, specific_run)
            if not metrics and generate_if_missing:
                print(f"\nMissing metrics for {mode} mode, run {specific_run}")
                print("Generating metrics from model...")
                if generate_metrics_from_model(exp_dir, mode, specific_run):
                    metrics = load_metrics(exp_dir, mode, specific_run)
            if metrics:
                all_modes_metrics[mode] = metrics
                all_modes_raw_metrics[mode] = [metrics]
                all_modes_error_data[mode] = load_error_data(exp_dir, mode, specific_run)
                if not all_modes_error_data[mode] and generate_if_missing:
                    if generate_metrics_from_model(exp_dir, mode, specific_run):
                        all_modes_error_data[mode] = load_error_data(exp_dir, mode, specific_run)
        elif n_runs > 1:
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
                
                # Load aggregated error data for distribution plots
                agg_errors = aggregate_error_data(exp_dir, mode, n_runs)
                if agg_errors:
                    all_modes_error_data[mode] = agg_errors
                else:
                    # Fallback to run 0 if aggregation fails (e.g. missing files)
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
        "mixed": ["#d35400", "#e59866"],  # Changed from Blue to Orange for better contrast
        "real": ["#922b21", "#ec7063"],
        "finetune": ["#7d3c98", "#bb8fce"],
    }

    main_colors = {
        "syn": "#27ae60",
        "mixed": "#e67e22",  # Changed from Blue to Orange
        "real": "#e74c3c",
        "finetune": "#9b59b6",
    }

    # Only modes that actually have metrics
    available_modes = [m for m in modes if m in all_modes_metrics]
    if not available_modes:
        print("Error: No modes with metrics available")
        return

    x = np.arange(len(available_modes))
    w = 0.38  # Increased from 0.35 to make bars wider for better label fit
    summary_mode = _pick_summary_mode(available_modes)
    output_dir.mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # FIGURE 1: Metrics (Panels a, b, c, d) - 2x2 Grid
    # =========================================================================
    
    # Use 2x2 grid, scale figsize up for readability
    fig1_size = figsize if figsize else (14, 12)
    fig1, axes1 = plt.subplots(2, 2, figsize=fig1_size)
    
    if specific_run is not None:
        title = f"Model Comparison: Run {specific_run}"
    elif n_runs > 1:
        title = f"Model Comparison: {n_runs} Training Runs (Mean ± Std)"
    else:
        title = "Model Comparison"
    fig1.suptitle(title, fontsize=18, fontweight="bold")

    # Panel (a): Mean Absolute Error
    ax = axes1[0, 0]
    ax.set_ylim(0, 220)  # Set limits first for correct label positioning
    for i, m in enumerate(available_modes):
        metrics = all_modes_metrics[m]
        mae_x = _get_mean(metrics, "mae_x")
        mae_y = _get_mean(metrics, "mae_y")

        bx = ax.bar(i - w/2, mae_x, w, color=model_colors[m][0], edgecolor='white', linewidth=1.5)
        by = ax.bar(i + w/2, mae_y, w, color=model_colors[m][1], edgecolor='white', linewidth=1.5)
        _label_inside_bar(ax, bx[0], f"{mae_x:.1f}")
        _label_inside_bar(ax, by[0], f"{mae_y:.1f}")
        
        if n_runs > 1 and m in all_modes_raw_metrics:
            raw_metrics = all_modes_raw_metrics[m]
            mae_x_values = [rm.get("mae_x", 0) for rm in raw_metrics if "mae_x" in rm]
            mae_y_values = [rm.get("mae_y", 0) for rm in raw_metrics if "mae_y" in rm]
            if len(mae_x_values) > 1 and len(mae_y_values) > 1:
                std_x = np.std(mae_x_values)
                std_y = np.std(mae_y_values)
                ax.errorbar(i - w/2, mae_x, yerr=std_x, fmt='none', ecolor='black', elinewidth=2, capsize=4, capthick=2)
                ax.errorbar(i + w/2, mae_y, yerr=std_y, fmt='none', ecolor='black', elinewidth=2, capsize=4, capthick=2)

    _add_panel_label(ax, "a")
    ax.set_title("Mean Absolute Error", fontsize=14)
    ax.set_ylabel("MAE (arcsec)", fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels([mode_labels.get(m, m) for m in available_modes], fontsize=11)
    ax.grid(axis="y", alpha=0.3)

    # Panel (b): Median Absolute Error
    ax = axes1[0, 1]
    ax.set_ylim(0, 70)  # Set limits first
    for i, m in enumerate(available_modes):
        metrics = all_modes_metrics[m]
        med_x = _get_mean(metrics, "median_ae_x")
        med_y = _get_mean(metrics, "median_ae_y")

        bx = ax.bar(i - w/2, med_x, w, color=model_colors[m][0], edgecolor='white', linewidth=1.5)
        by = ax.bar(i + w/2, med_y, w, color=model_colors[m][1], edgecolor='white', linewidth=1.5)
        _label_inside_bar(ax, bx[0], f"{med_x:.1f}")
        _label_inside_bar(ax, by[0], f"{med_y:.1f}")
        
        if n_runs > 1 and m in all_modes_raw_metrics:
            raw_metrics = all_modes_raw_metrics[m]
            med_x_values = [rm.get("median_ae_x", 0) for rm in raw_metrics if "median_ae_x" in rm]
            med_y_values = [rm.get("median_ae_y", 0) for rm in raw_metrics if "median_ae_y" in rm]
            if len(med_x_values) > 1 and len(med_y_values) > 1:
                std_x = np.std(med_x_values)
                std_y = np.std(med_y_values)
                ax.errorbar(i - w/2, med_x, yerr=std_x, fmt='none', ecolor='black', elinewidth=2, capsize=4, capthick=2)
                ax.errorbar(i + w/2, med_y, yerr=std_y, fmt='none', ecolor='black', elinewidth=2, capsize=4, capthick=2)

    _add_panel_label(ax, "b")
    ax.set_title("Median Absolute Error", fontsize=14)
    ax.set_ylabel("Median AE (arcsec)", fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels([mode_labels.get(m, m) for m in available_modes], fontsize=11)
    ax.grid(axis="y", alpha=0.3)

    # Panel (c): Euclidean Distance Error (Moved to row 2, col 1)
    ax = axes1[1, 0]
    for i, m in enumerate(available_modes):
        metrics = all_modes_metrics[m]
        euc_x = _get_mean(metrics, "mae_x")
        euc_y = _get_mean(metrics, "mae_y")
        euc_eucl = _get_mean(metrics, "euclidean_mean")

        bx = ax.bar(i - w, euc_x, w, color=model_colors[m][0], edgecolor='white', linewidth=1.5)
        by = ax.bar(i, euc_y, w, color=model_colors[m][1], edgecolor='white', linewidth=1.5)
        be = ax.bar(i + w, euc_eucl, w, color=main_colors[m], edgecolor='white', linewidth=1.5, alpha=0.7)
        _label_inside_bar(ax, bx[0], f"{euc_x:.1f}")
        _label_inside_bar(ax, by[0], f"{euc_y:.1f}")
        _label_inside_bar(ax, be[0], f"{euc_eucl:.1f}")

        if n_runs > 1 and m in all_modes_raw_metrics:
            raw_metrics = all_modes_raw_metrics[m]
            ex = [rm.get("mae_x", 0) for rm in raw_metrics if "mae_x" in rm]
            ey = [rm.get("mae_y", 0) for rm in raw_metrics if "mae_y" in rm]
            ee = [rm.get("euclidean_mean", 0) for rm in raw_metrics if "euclidean_mean" in rm]
            if len(ex) > 1:
                ax.errorbar(i - w, np.mean(ex), yerr=np.std(ex), fmt='none', ecolor='black', elinewidth=2, capsize=4, capthick=2)
            if len(ey) > 1:
                ax.errorbar(i, np.mean(ey), yerr=np.std(ey), fmt='none', ecolor='black', elinewidth=2, capsize=4, capthick=2)
            if len(ee) > 1:
                ax.errorbar(i + w, np.mean(ee), yerr=np.std(ee), fmt='none', ecolor='black', elinewidth=2, capsize=4, capthick=2)

    _add_panel_label(ax, "c")
    ax.set_title("Euclidean Distance Error (X / Y / Eucl)", fontsize=14)
    ax.set_ylabel("Error (arcsec)", fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels([mode_labels.get(m, m) for m in available_modes], fontsize=11)
    ax.grid(axis="y", alpha=0.3)

    # Panel (d): 95th Percentile Error (Moved to row 2, col 2)
    ax = axes1[1, 1]
    for i, m in enumerate(available_modes):
        metrics = all_modes_metrics[m]
        p95_x = _get_mean(metrics, "q95_x")
        p95_y = _get_mean(metrics, "q95_y")
        p95_eucl = _get_mean(metrics, "euclidean_q95")

        bx = ax.bar(i - w, p95_x, w, color=model_colors[m][0], edgecolor='white', linewidth=1.5)
        by = ax.bar(i, p95_y, w, color=model_colors[m][1], edgecolor='white', linewidth=1.5)
        be = ax.bar(i + w, p95_eucl, w, color=main_colors[m], edgecolor='white', linewidth=1.5, alpha=0.7)
        _label_inside_bar(ax, bx[0], f"{p95_x:.1f}")
        _label_inside_bar(ax, by[0], f"{p95_y:.1f}")
        _label_inside_bar(ax, be[0], f"{p95_eucl:.1f}")

        if n_runs > 1 and m in all_modes_raw_metrics:
            raw_metrics = all_modes_raw_metrics[m]
            p95_x_values = [rm.get("q95_x", 0) for rm in raw_metrics if "q95_x" in rm]
            p95_y_values = [rm.get("q95_y", 0) for rm in raw_metrics if "q95_y" in rm]
            p95_e_values = [rm.get("euclidean_q95", 0) for rm in raw_metrics if "euclidean_q95" in rm]
            if len(p95_x_values) > 1:
                ax.errorbar(i - w, np.mean(p95_x_values), yerr=np.std(p95_x_values), fmt='none', ecolor='black', elinewidth=2, capsize=4, capthick=2)
            if len(p95_y_values) > 1:
                ax.errorbar(i, np.mean(p95_y_values), yerr=np.std(p95_y_values), fmt='none', ecolor='black', elinewidth=2, capsize=4, capthick=2)
            if len(p95_e_values) > 1:
                ax.errorbar(i + w, np.mean(p95_e_values), yerr=np.std(p95_e_values), fmt='none', ecolor='black', elinewidth=2, capsize=4, capthick=2)

    _add_panel_label(ax, "d")
    ax.set_title("95th Percentile Error (X / Y / Eucl)", fontsize=14)
    ax.set_ylabel("P95 Error (arcsec)", fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels([mode_labels.get(m, m) for m in available_modes], fontsize=11)
    ax.grid(axis="y", alpha=0.3)

    plt.figure(fig1.number)
    plt.tight_layout()
    # Reserve top space for title and labels
    top_frac = 0.12
    for ax in axes1.flat:
        pos = ax.get_position()
        delta = pos.height * top_frac
        ax.set_position([pos.x0, pos.y0, pos.width, pos.height - delta])

    if specific_run is not None:
        path1 = output_dir / f"comparison_metrics_run{specific_run}.png"
        path1_pdf = output_dir / f"comparison_metrics_run{specific_run}.pdf"
    else:
        path1 = output_dir / "comparison_metrics.png"
        path1_pdf = output_dir / "comparison_metrics.pdf"
    plt.savefig(path1, dpi=300, bbox_inches="tight")
    plt.savefig(path1_pdf, bbox_inches="tight")
    print(f"Saved: {path1} and {path1_pdf}")
    plt.close(fig1)

    # =========================================================================
    # FIGURE 2: Analysis (Panels e, f) - 1x2 Grid
    # =========================================================================
    
    # 1x2 grid
    fig2_size = (16, 6) # Wide format
    fig2, axes2 = plt.subplots(1, 2, figsize=fig2_size)
    # fig2.suptitle("Error Distribution and Spatial Analysis", fontsize=16, fontweight="bold")

    # Panel (e): Error Distribution
    ax = axes2[0]
    bins = np.linspace(0, 500, 60)
    for m in available_modes:
        ed = all_modes_error_data.get(m) or {}
        errors = ed.get("euclidean_errors", None)
        if errors is None:
            continue
            
        # Make "mixed" appear on top of others
        zorder = 10 if m == "mixed" else 3
            
        ax.hist(errors, bins=bins, alpha=0.4, histtype="stepfilled", color=main_colors[m], edgecolor=main_colors[m], linewidth=1.5, label=mode_labels.get(m, m), zorder=zorder)
        ax.hist(errors, bins=bins, alpha=1.0, histtype="step", color=main_colors[m], edgecolor=main_colors[m], linewidth=1.5, label=None, zorder=zorder)

    ax.set_yscale("log")
    _add_panel_label(ax, "e")
    ax.set_title("Error Distribution", fontsize=14)
    ax.set_xlabel("Euclidean Error (arcsec)", fontsize=12)
    ax.set_ylabel("Count (log)", fontsize=12)
    ax.legend(loc="upper right", fontsize=10, framealpha=0.95, edgecolor="gray")
    ax.grid(alpha=0.3)

    # Panel (f): Spatial Error Map
    ax = axes2[1]
    spatial_mode = "mixed" if "mixed" in available_modes else summary_mode
    
    if specific_run is not None:
        spatial_data = load_spatial_data(exp_dir, spatial_mode, specific_run)
        if not spatial_data and generate_if_missing:
            if generate_metrics_from_model(exp_dir, spatial_mode, specific_run):
                spatial_data = load_spatial_data(exp_dir, spatial_mode, specific_run)
    elif n_runs > 1:
        # Aggregate spatial data from all runs
        spatial_data = aggregate_spatial_data(exp_dir, spatial_mode, n_runs)
        if not spatial_data and generate_if_missing:
            # Try generating for all runs
            for r in range(n_runs):
                generate_metrics_from_model(exp_dir, spatial_mode, r)
            spatial_data = aggregate_spatial_data(exp_dir, spatial_mode, n_runs)
    else:
        spatial_data = load_spatial_data(exp_dir, spatial_mode, None)
        if not spatial_data and generate_if_missing:
            if generate_metrics_from_model(exp_dir, spatial_mode, None):
                spatial_data = load_spatial_data(exp_dir, spatial_mode, None)

    if spatial_data:
        true_x = spatial_data.get("true_x")
        true_y = spatial_data.get("true_y")
        euclidean_error = spatial_data.get("euclidean_error")

        if true_x is None or true_y is None or euclidean_error is None:
            ax.text(0.5, 0.5, "Spatial error map data incomplete", ha="center", va="center", transform=ax.transAxes, fontsize=12)
            _add_panel_label(ax, "f")
            ax.set_title("Spatial Error Map", fontsize=14)
        else:
            # ax.scatter(true_x, true_y, s=4, alpha=0.05, color="black", label="True")
            hb = ax.hexbin(
                true_x, true_y,
                C=euclidean_error,
                reduce_C_function=np.mean,
                gridsize=55,
                mincnt=1,
                vmin=0,
                vmax=100,
                cmap="plasma",
            )
            cb = fig2.colorbar(hb, ax=ax)
            cb.set_label("Mean Euclidean Error (arcsec)", fontsize=11)

            _add_panel_label(ax, "f")
            if spatial_mode == "mixed":
                ax.set_title("Spatial Error Map (Syn+Real Mixed)", fontsize=14)
            else:
                ax.set_title(f"Spatial Error Map ({mode_labels.get(spatial_mode, spatial_mode)})", fontsize=14)

            ax.set_xlabel("X (arcsec)", fontsize=12)
            ax.set_ylabel("Y (arcsec)", fontsize=12)
            ax.set_aspect("equal")
            ax.grid(alpha=0.3)
    else:
        ax.text(0.5, 0.5, "Spatial error map data not available", ha="center", va="center", transform=ax.transAxes, fontsize=12)
        _add_panel_label(ax, "f")
        ax.set_title("Spatial Error Map", fontsize=14)

    plt.figure(fig2.number)
    plt.tight_layout()
    # Reserve space for labels
    top_frac = 0.12
    top_frac_spatial = 0.04
    
    # Adjust axes2[0] (distribution)
    pos0 = axes2[0].get_position()
    axes2[0].set_position([pos0.x0, pos0.y0, pos0.width, pos0.height - (pos0.height * top_frac)])
    
    # Adjust axes2[1] (spatial)
    pos1 = axes2[1].get_position()
    axes2[1].set_position([pos1.x0, pos1.y0, pos1.width, pos1.height - (pos1.height * top_frac_spatial)])

    if specific_run is not None:
        path2 = output_dir / f"comparison_analysis_run{specific_run}.png"
        path2_pdf = output_dir / f"comparison_analysis_run{specific_run}.pdf"
    else:
        path2 = output_dir / "comparison_analysis.png"
        path2_pdf = output_dir / "comparison_analysis.pdf"
    plt.savefig(path2, dpi=300, bbox_inches="tight")
    plt.savefig(path2_pdf, bbox_inches="tight")
    print(f"Saved: {path2} and {path2_pdf}")
    plt.close(fig2)


def create_custom_plot(
    exp_dir: Path,
    modes: List[str],
    n_runs: int,
    output_dir: Path,
    panels: List[str] = ["a", "d"],
    figsize: tuple = (9.6, 7.2),  # 960x720 pixels at 100 DPI
    generate_if_missing: bool = False,
    specific_run: Optional[int] = None,
):
    """
    Create custom plot with selected panels.
    
    Parameters:
    -----------
    panels : List[str]
        List of panel IDs to include: "a" (MAE), "b" (Median), "c" (Euclidean), 
        "d" (P95), "e" (Distribution), "f" (Spatial)
    figsize : tuple
        Figure size in inches (width, height)
    """
    # Load metrics (same as create_comparison_plot)
    all_modes_metrics: Dict[str, Dict] = {}
    all_modes_error_data: Dict[str, Optional[Dict]] = {}
    all_modes_raw_metrics: Dict[str, List[Dict]] = {}

    for mode in modes:
        if specific_run is not None:
            metrics = load_metrics(exp_dir, mode, specific_run)
            if not metrics and generate_if_missing:
                if generate_metrics_from_model(exp_dir, mode, specific_run):
                    metrics = load_metrics(exp_dir, mode, specific_run)
            if metrics:
                all_modes_metrics[mode] = metrics
                all_modes_raw_metrics[mode] = [metrics]
                all_modes_error_data[mode] = load_error_data(exp_dir, mode, specific_run)
        elif n_runs > 1:
            all_metrics: List[Dict] = []
            for run_id in range(n_runs):
                metrics = load_metrics(exp_dir, mode, run_id)
                if metrics:
                    all_metrics.append(metrics)
            if all_metrics:
                all_modes_metrics[mode] = aggregate_metrics(all_metrics)
                all_modes_raw_metrics[mode] = all_metrics
                
                # Load aggregated error data
                agg_errors = aggregate_error_data(exp_dir, mode, n_runs)
                if agg_errors:
                    all_modes_error_data[mode] = agg_errors
                else:
                    all_modes_error_data[mode] = load_error_data(exp_dir, mode, 0)
        else:
            metrics = load_metrics(exp_dir, mode, None)
            if metrics:
                all_modes_metrics[mode] = metrics
                all_modes_raw_metrics[mode] = [metrics]
                all_modes_error_data[mode] = load_error_data(exp_dir, mode, None)

    if not all_modes_metrics:
        print("Error: No metrics found for any mode")
        return

    mode_labels = {
        "syn": "Synthetic\nOnly",
        "mixed": "Syn+Real\nMixed",
        "real": "Real Only\n(Baseline)",
        "finetune": "Syn+Real\nFinetune",
    }

    model_colors = {
        "syn": ["#1e8449", "#58d68d"],
        "mixed": ["#d35400", "#e59866"],  # Changed from Blue to Orange for better contrast
        "real": ["#922b21", "#ec7063"],
        "finetune": ["#7d3c98", "#bb8fce"],
    }

    main_colors = {
        "syn": "#27ae60",
        "mixed": "#e67e22",  # Changed from Blue to Orange
        "real": "#e74c3c",
        "finetune": "#9b59b6",
    }

    available_modes = [m for m in modes if m in all_modes_metrics]
    if not available_modes:
        print("Error: No modes with metrics available")
        return

    # Determine layout based on panels
    n_panels = len(panels)
    if n_panels == 1:
        rows, cols = 1, 1
    elif n_panels == 2:
        rows, cols = 1, 2
    elif n_panels <= 4:
        rows, cols = 2, 2
    else:
        rows, cols = 2, 3

    fig, axes = plt.subplots(rows, cols, figsize=figsize)
    if n_panels == 1:
        axes = [axes]
    elif rows == 1 or cols == 1:
        axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]
    else:
        axes = axes.flatten()

    x = np.arange(len(available_modes))
    w = 0.38
    summary_mode = _pick_summary_mode(available_modes)

    # Map panel IDs to positions (canonical order: a,b,c,d,e,f -> 0..5)
    canonical = ["a", "b", "c", "d", "e", "f"]
    panel_positions = {}
    if n_panels == 2:
        if "a" in panels:
            panel_positions["a"] = 0
            other_panel = [p for p in panels if p != "a"][0]
            panel_positions[other_panel] = 1
        else:
            for i, panel_id in enumerate(panels):
                panel_positions[panel_id] = i
    elif n_panels == 4 and rows == 2 and cols == 2:
        position_map = {"a": 0, "d": 1, "b": 2, "f": 3}
        for panel_id in panels:
            panel_positions[panel_id] = position_map.get(panel_id, panels.index(panel_id))
    elif n_panels == 6 and rows == 2 and cols == 3:
        for panel_id in panels:
            panel_positions[panel_id] = canonical.index(panel_id)
    else:
        for i, panel_id in enumerate(panels):
            panel_positions[panel_id] = i

    # Panel (a): Mean Absolute Error
    if "a" in panels:
        pos = panel_positions.get("a", panels.index("a")) if panel_positions else panels.index("a")
        ax = axes[pos] if n_panels > 1 else axes
        for i, m in enumerate(available_modes):
            metrics = all_modes_metrics[m]
            mae_x = _get_mean(metrics, "mae_x")
            mae_y = _get_mean(metrics, "mae_y")

            bx = ax.bar(i - w/2, mae_x, w, color=model_colors[m][0], edgecolor='white', linewidth=1.5)
            by = ax.bar(i + w/2, mae_y, w, color=model_colors[m][1], edgecolor='white', linewidth=1.5)
            _label_inside_bar(ax, bx[0], f"{mae_x:.1f}")
            _label_inside_bar(ax, by[0], f"{mae_y:.1f}")
            
            if n_runs > 1 and m in all_modes_raw_metrics:
                raw_metrics = all_modes_raw_metrics[m]
                mae_x_values = [rm.get("mae_x", 0) for rm in raw_metrics if "mae_x" in rm]
                mae_y_values = [rm.get("mae_y", 0) for rm in raw_metrics if "mae_y" in rm]
                
                if len(mae_x_values) > 1 and len(mae_y_values) > 1:
                    std_x = np.std(mae_x_values)
                    std_y = np.std(mae_y_values)
                    ax.errorbar(i - w/2, mae_x, yerr=std_x, fmt='none', ecolor='black', 
                               elinewidth=2, capsize=4, capthick=2)
                    ax.errorbar(i + w/2, mae_y, yerr=std_y, fmt='none', ecolor='black', 
                               elinewidth=2, capsize=4, capthick=2)

        _add_panel_label(ax, "a")
        ax.set_title("Mean Absolute Error")
        ax.set_ylabel("MAE (arcsec)")
        ax.set_xticks(x)
        ax.set_xticklabels([mode_labels.get(m, m) for m in available_modes])
        ax.set_ylim(0, 220)
        ax.grid(axis="y", alpha=0.3)

    # Panel (b): Median Absolute Error
    if "b" in panels:
        pos = panel_positions.get("b", panels.index("b")) if panel_positions else panels.index("b")
        ax = axes[pos] if n_panels > 1 else axes
        for i, m in enumerate(available_modes):
            metrics = all_modes_metrics[m]
            med_x = _get_mean(metrics, "median_ae_x")
            med_y = _get_mean(metrics, "median_ae_y")

            bx = ax.bar(i - w/2, med_x, w, color=model_colors[m][0], edgecolor='white', linewidth=1.5)
            by = ax.bar(i + w/2, med_y, w, color=model_colors[m][1], edgecolor='white', linewidth=1.5)
            _label_inside_bar(ax, bx[0], f"{med_x:.1f}")
            _label_inside_bar(ax, by[0], f"{med_y:.1f}")
            
            if n_runs > 1 and m in all_modes_raw_metrics:
                raw_metrics = all_modes_raw_metrics[m]
                med_x_values = [rm.get("median_ae_x", 0) for rm in raw_metrics if "median_ae_x" in rm]
                med_y_values = [rm.get("median_ae_y", 0) for rm in raw_metrics if "median_ae_y" in rm]
                
                if len(med_x_values) > 1 and len(med_y_values) > 1:
                    std_x = np.std(med_x_values)
                    std_y = np.std(med_y_values)
                    ax.errorbar(i - w/2, med_x, yerr=std_x, fmt='none', ecolor='black', 
                               elinewidth=2, capsize=4, capthick=2)
                    ax.errorbar(i + w/2, med_y, yerr=std_y, fmt='none', ecolor='black', 
                               elinewidth=2, capsize=4, capthick=2)

        _add_panel_label(ax, "b")
        ax.set_title("Median Absolute Error")
        ax.set_ylabel("Median AE (arcsec)")
        ax.set_xticks(x)
        ax.set_xticklabels([mode_labels.get(m, m) for m in available_modes])
        ax.set_ylim(0, 70)
        ax.grid(axis="y", alpha=0.3)

    # Panel (c): Euclidean Distance Error (X / Y / Eucl)
    if "c" in panels:
        pos = panel_positions.get("c", panels.index("c")) if panel_positions else panels.index("c")
        ax = axes[pos] if n_panels > 1 else axes
        for i, m in enumerate(available_modes):
            metrics = all_modes_metrics[m]
            euc_x = _get_mean(metrics, "mae_x")
            euc_y = _get_mean(metrics, "mae_y")
            euc_eucl = _get_mean(metrics, "euclidean_mean")

            bx = ax.bar(i - w, euc_x, w, color=model_colors[m][0], edgecolor='white', linewidth=1.5)
            by = ax.bar(i, euc_y, w, color=model_colors[m][1], edgecolor='white', linewidth=1.5)
            be = ax.bar(i + w, euc_eucl, w, color=main_colors[m], edgecolor='white', linewidth=1.5, alpha=0.7)
            _label_inside_bar(ax, bx[0], f"{euc_x:.1f}")
            _label_inside_bar(ax, by[0], f"{euc_y:.1f}")
            _label_inside_bar(ax, be[0], f"{euc_eucl:.1f}")

            if n_runs > 1 and m in all_modes_raw_metrics:
                raw_metrics = all_modes_raw_metrics[m]
                ex = [rm.get("mae_x", 0) for rm in raw_metrics if "mae_x" in rm]
                ey = [rm.get("mae_y", 0) for rm in raw_metrics if "mae_y" in rm]
                ee = [rm.get("euclidean_mean", 0) for rm in raw_metrics if "euclidean_mean" in rm]
                if len(ex) > 1:
                    ax.errorbar(i - w, np.mean(ex), yerr=np.std(ex), fmt='none', ecolor='black', elinewidth=2, capsize=4, capthick=2)
                if len(ey) > 1:
                    ax.errorbar(i, np.mean(ey), yerr=np.std(ey), fmt='none', ecolor='black', elinewidth=2, capsize=4, capthick=2)
                if len(ee) > 1:
                    ax.errorbar(i + w, np.mean(ee), yerr=np.std(ee), fmt='none', ecolor='black', elinewidth=2, capsize=4, capthick=2)

        _add_panel_label(ax, "c")
        ax.set_title("Euclidean Distance Error (X / Y / Eucl)")
        ax.set_ylabel("Error (arcsec)")
        ax.set_xticks(x)
        ax.set_xticklabels([mode_labels.get(m, m) for m in available_modes])
        ax.grid(axis="y", alpha=0.3)

    # Panel (d): 95th Percentile Error (X / Y / Eucl)
    if "d" in panels:
        pos = panel_positions.get("d", panels.index("d")) if panel_positions else panels.index("d")
        ax = axes[pos] if n_panels > 1 else axes
        for i, m in enumerate(available_modes):
            metrics = all_modes_metrics[m]
            p95_x = _get_mean(metrics, "q95_x")
            p95_y = _get_mean(metrics, "q95_y")
            p95_eucl = _get_mean(metrics, "euclidean_q95")

            bx = ax.bar(i - w, p95_x, w, color=model_colors[m][0], edgecolor='white', linewidth=1.5)
            by = ax.bar(i, p95_y, w, color=model_colors[m][1], edgecolor='white', linewidth=1.5)
            be = ax.bar(i + w, p95_eucl, w, color=main_colors[m], edgecolor='white', linewidth=1.5, alpha=0.7)
            _label_inside_bar(ax, bx[0], f"{p95_x:.1f}")
            _label_inside_bar(ax, by[0], f"{p95_y:.1f}")
            _label_inside_bar(ax, be[0], f"{p95_eucl:.1f}")

            if n_runs > 1 and m in all_modes_raw_metrics:
                raw_metrics = all_modes_raw_metrics[m]
                p95_x_values = [rm.get("q95_x", 0) for rm in raw_metrics if "q95_x" in rm]
                p95_y_values = [rm.get("q95_y", 0) for rm in raw_metrics if "q95_y" in rm]
                p95_e_values = [rm.get("euclidean_q95", 0) for rm in raw_metrics if "euclidean_q95" in rm]
                if len(p95_x_values) > 1:
                    ax.errorbar(i - w, np.mean(p95_x_values), yerr=np.std(p95_x_values), fmt='none', ecolor='black', elinewidth=2, capsize=4, capthick=2)
                if len(p95_y_values) > 1:
                    ax.errorbar(i, np.mean(p95_y_values), yerr=np.std(p95_y_values), fmt='none', ecolor='black', elinewidth=2, capsize=4, capthick=2)
                if len(p95_e_values) > 1:
                    ax.errorbar(i + w, np.mean(p95_e_values), yerr=np.std(p95_e_values), fmt='none', ecolor='black', elinewidth=2, capsize=4, capthick=2)

        _add_panel_label(ax, "d")
        ax.set_title("95th Percentile Error (X / Y / Eucl)")
        ax.set_ylabel("P95 Error (arcsec)")
        ax.set_xticks(x)
        ax.set_xticklabels([mode_labels.get(m, m) for m in available_modes])
        ax.grid(axis="y", alpha=0.3)

    # Panel (e): Error Distribution
    if "e" in panels:
        pos = panel_positions.get("e", panels.index("e")) if panel_positions else panels.index("e")
        ax = axes[pos] if n_panels > 1 else axes
        bins = np.linspace(0, 500, 60)
        for m in available_modes:
            ed = all_modes_error_data.get(m) or {}
            errors = ed.get("euclidean_errors", None)
            if errors is None:
                continue
                
            # Make "mixed" appear on top of others
            zorder = 10 if m == "mixed" else 3
            
            ax.hist(
                errors,
                bins=bins,
                alpha=0.4,
                histtype="stepfilled",
                color=main_colors[m],
                edgecolor=main_colors[m],
                linewidth=1.5,
                label=mode_labels.get(m, m),
                zorder=zorder,
            )
            ax.hist(
                errors,
                bins=bins,
                alpha=1.0,
                histtype="step",
                color=main_colors[m],
                edgecolor=main_colors[m],
                linewidth=1.5,
                label=None,
                zorder=zorder,
            )
        ax.set_yscale("log")
        _add_panel_label(ax, "e")
        ax.set_title("Error Distribution")
        ax.set_xlabel("Euclidean Error (arcsec)")
        ax.set_ylabel("Count (log)")
        ax.legend(loc="upper right", fontsize=9, framealpha=0.95, edgecolor="gray")
        ax.grid(alpha=0.3)

    # Panel (f): Spatial Error Map
    if "f" in panels:
        pos = panel_positions.get("f", panels.index("f")) if panel_positions else panels.index("f")
        ax = axes[pos] if n_panels > 1 else axes
        spatial_mode = "mixed" if "mixed" in available_modes else summary_mode
        
        if specific_run is not None:
            spatial_data = load_spatial_data(exp_dir, spatial_mode, specific_run)
            if not spatial_data and generate_if_missing:
                if generate_metrics_from_model(exp_dir, spatial_mode, specific_run):
                    spatial_data = load_spatial_data(exp_dir, spatial_mode, specific_run)
        elif n_runs > 1:
            spatial_data = aggregate_spatial_data(exp_dir, spatial_mode, n_runs)
            if not spatial_data and generate_if_missing:
                for r in range(n_runs):
                    generate_metrics_from_model(exp_dir, spatial_mode, r)
                spatial_data = aggregate_spatial_data(exp_dir, spatial_mode, n_runs)
        else:
            spatial_data = load_spatial_data(exp_dir, spatial_mode, None)
            if not spatial_data and generate_if_missing:
                if generate_metrics_from_model(exp_dir, spatial_mode, None):
                    spatial_data = load_spatial_data(exp_dir, spatial_mode, None)

        if spatial_data:
            true_x = spatial_data.get("true_x")
            true_y = spatial_data.get("true_y")
            euclidean_error = spatial_data.get("euclidean_error")

            if true_x is None or true_y is None or euclidean_error is None:
                ax.text(0.5, 0.5, "Spatial error map data incomplete", ha="center", va="center",
                        transform=ax.transAxes, fontsize=12)
                _add_panel_label(ax, "f")
                ax.set_title("Spatial Error Map")
            else:
                # ax.scatter(true_x, true_y, s=4, alpha=0.05, color="black", label="True")

                hb = ax.hexbin(
                    true_x, true_y,
                    C=euclidean_error,
                    reduce_C_function=np.mean,
                    gridsize=55,
                    mincnt=1,
                    vmin=0,
                    vmax=100,
                )
                cb = fig.colorbar(hb, ax=ax)
                cb.set_label("Mean Euclidean Error (arcsec)")

                _add_panel_label(ax, "f")
                if spatial_mode == "mixed":
                    ax.set_title("Spatial Error Map (Syn+Real Mixed)")
                else:
                    ax.set_title(f"Spatial Error Map ({mode_labels.get(spatial_mode, spatial_mode)})")

                ax.set_xlabel("X (arcsec)")
                ax.set_ylabel("Y (arcsec)")
                ax.set_aspect("equal")
                ax.grid(alpha=0.3)
        else:
            ax.text(0.5, 0.5, "Spatial error map data not available", ha="center", va="center",
                    transform=ax.transAxes, fontsize=12)
            _add_panel_label(ax, "f")
            ax.set_title("Spatial Error Map")

    # Hide unused subplots
    used_positions = set(panel_positions.values()) if panel_positions else set(range(len(panels)))
    for idx in range(len(axes)):
        if idx not in used_positions:
            axes[idx].set_visible(False)

    plt.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)
    panel_str = "_".join(panels)
    if specific_run is not None:
        output_path = output_dir / f"custom_plot_{panel_str}_run{specific_run}.png"
        output_path_pdf = output_dir / f"custom_plot_{panel_str}_run{specific_run}.pdf"
    else:
        output_path = output_dir / f"custom_plot_{panel_str}.png"
        output_path_pdf = output_dir / f"custom_plot_{panel_str}.pdf"
    plt.savefig(output_path, dpi=100, bbox_inches='tight')  # 100 DPI: 9.6" * 100 = 960px
    plt.savefig(output_path_pdf, bbox_inches='tight')
    print(f"Saved: {output_path} and {output_path_pdf}")
    plt.close()


def create_mae_spatial_plot(
    exp_dir: Path,
    modes: List[str],
    n_runs: int,
    output_dir: Path,
    generate_if_missing: bool = False,
    specific_run: Optional[int] = None,
):
    """
    Create a spatial error map using Mean Absolute Error (MAE) instead of Euclidean distance.
    """
    # Only modes that actually have metrics
    available_modes = [m for m in modes if load_metrics(exp_dir, m, 0 if n_runs > 1 else None)]
    if not available_modes:
        print("Error: No modes with metrics available")
        return

    summary_mode = _pick_summary_mode(available_modes)
    spatial_mode = "mixed" if "mixed" in available_modes else summary_mode
    
    # Load data
    if specific_run is not None:
        spatial_data = load_spatial_data(exp_dir, spatial_mode, specific_run)
    elif n_runs > 1:
        spatial_data = aggregate_spatial_data(exp_dir, spatial_mode, n_runs)
    else:
        spatial_data = load_spatial_data(exp_dir, spatial_mode, None)
        
    if not spatial_data:
        print("Error: Could not load spatial data for MAE plot")
        return
        
    true_x = spatial_data.get("true_x")
    true_y = spatial_data.get("true_y")
    pred_x = spatial_data.get("pred_x")
    pred_y = spatial_data.get("pred_y")
    
    if true_x is None or pred_x is None:
        print("Error: Prediction data (pred_x/pred_y) not found in spatial files. Cannot calculate MAE.")
        return
        
    # Calculate MAE per point: (abs(err_x) + abs(err_y)) / 2
    abs_err_x = np.abs(true_x - pred_x)
    abs_err_y = np.abs(true_y - pred_y)
    mae_per_point = (abs_err_x + abs_err_y) / 2
    
    # Create Plot
    fig, ax = plt.subplots(figsize=(8, 7))
    
    hb = ax.hexbin(
        true_x, true_y,
        C=mae_per_point,
        reduce_C_function=np.mean,
        gridsize=55,
        mincnt=1,
        vmin=0,
        vmax=100,
        cmap="plasma",
    )
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label("Mean Absolute Error (arcsec)", fontsize=12)
    
    mode_labels = {
        "syn": "Synthetic Only",
        "mixed": "Syn+Real Mixed",
        "real": "Real Only",
        "finetune": "Syn+Real Finetune",
    }
    
    ax.set_title(f"Spatial MAE Map ({mode_labels.get(spatial_mode, spatial_mode)})", fontsize=14, fontweight="bold")
    ax.set_xlabel("X (arcsec)", fontsize=12)
    ax.set_ylabel("Y (arcsec)", fontsize=12)
    ax.set_aspect("equal")
    ax.grid(alpha=0.3)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    if specific_run is not None:
        out_name = f"spatial_mae_run{specific_run}"
    else:
        out_name = "spatial_mae"
        
    path_png = output_dir / f"{out_name}.png"
    path_pdf = output_dir / f"{out_name}.pdf"
    
    plt.savefig(path_png, dpi=300, bbox_inches="tight")
    plt.savefig(path_pdf, bbox_inches="tight")
    print(f"Saved MAE spatial plot: {path_png} and {path_pdf}")
    plt.close()


def find_available_modes(exp_dir: Path) -> List[str]:
    """Find available training modes in experiment directory."""
    modes: List[str] = []
    for mode in ["real", "syn", "mixed", "finetune"]:
        mode_dir = exp_dir / mode
        if mode_dir.exists():
            # Check for metrics files with appropriate naming
            if mode == "finetune":
                if any(mode_dir.glob("finetune*_extended_metrics.json")):
                    modes.append(mode)
            else:
                if any(mode_dir.glob(f"train_on_{mode}*_extended_metrics.json")):
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
    parser.add_argument(
        "--run",
        type=int,
        default=None,
        help="Plot only a specific run ID (e.g., --run 3 for run 3)",
    )
    parser.add_argument(
        "--panels",
        type=str,
        nargs="+",
        default=None,
        help="Select specific panels: a (MAE), b (Median), c (Euclidean), d (P95), e (Distribution), f (Spatial). Example: --panels a d",
    )
    parser.add_argument(
        "--figsize",
        type=float,
        nargs=2,
        default=None,
        help="Figure size in inches (width height). Example: --figsize 9.6 7.2",
    )

    parser.add_argument(
        "--mae-map",
        action="store_true",
        help="Generate an additional spatial map using Mean Absolute Error instead of Euclidean distance",
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

    # Detect n_runs (only if not plotting a specific run)
    if args.run is not None:
        n_runs = 1
        print(f"\nPlotting experiment: {experiment_name}")
        print(f"Modes: {', '.join(available_modes)}")
        print(f"Run: {args.run} (specific run only)")
    else:
        n_runs = 1
        for run_id in range(10):
            p = exp_dir / available_modes[0] / f"train_on_{available_modes[0]}_run{run_id}_extended_metrics.json"
            if p.exists():
                n_runs = max(n_runs, run_id + 1)
        print(f"\nPlotting experiment: {experiment_name}")
        print(f"Modes: {', '.join(available_modes)}")
        print(f"Runs: {n_runs}")

    output_dir = exp_dir / "plots"
    
    # Use custom plot if panels are specified
    if args.panels:
        # Validate panels
        valid_panels = ["a", "b", "c", "d", "e", "f"]
        panels = [p.lower() for p in args.panels if p.lower() in valid_panels]
        if not panels:
            print("Error: No valid panels specified. Valid panels: a, b, c, d, e, f")
            sys.exit(1)
        
        # Set default figsize if not provided (960x720 pixels = 9.6x7.2 inches at 100 DPI)
        if args.figsize:
            figsize = tuple(args.figsize)
        else:
            figsize = (9.6, 7.2)
        
        create_custom_plot(
            exp_dir,
            available_modes,
            n_runs,
            output_dir,
            panels=panels,
            figsize=figsize,
            generate_if_missing=args.generate_metrics,
            specific_run=args.run,
        )
        
        panel_str = "_".join(panels)
        if args.run is not None:
            print(f"\nDone! Plot saved to: {output_dir}/custom_plot_{panel_str}_run{args.run}.png")
        else:
            print(f"\nDone! Plot saved to: {output_dir}/custom_plot_{panel_str}.png")
    else:
        # Use standard comparison plot (all panels a–f, labels on the sides)
        figsize = tuple(args.figsize) if args.figsize else None
        create_comparison_plot(
            exp_dir,
            available_modes,
            n_runs,
            output_dir,
            generate_if_missing=args.generate_metrics,
            specific_run=args.run,
            figsize=figsize,
        )

        if args.run is not None:
            print(f"\nDone! Plot saved to: {output_dir}/comparison_performance_run{args.run}.png")
            
        else:
            print(f"\nDone! Plot saved to: {output_dir}/comparison_performance.png")

    if args.mae_map:
        create_mae_spatial_plot(
            exp_dir,
            available_modes,
            n_runs,
            output_dir,
            generate_if_missing=args.generate_metrics,
            specific_run=args.run,
        )


if __name__ == "__main__":
    main()