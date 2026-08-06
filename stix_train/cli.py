#!/usr/bin/env python3
"""STIX localization training CLI (`stix-train`)."""

import os

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
os.environ["TF_ENABLE_ONEDNN_OPTS"] = "0"

import argparse
from pathlib import Path

import numpy as np

from .config import TrainConfig, EXPERIMENTS_DIR, PROJECT_ROOT
from .data import load_data, ensure_synthetic_data
from .training import train_single_mode, save_experiment_config


def parse_args() -> TrainConfig:
    parser = argparse.ArgumentParser(description="Train STIX localization models")
    parser.add_argument("--experiment", type=str, required=True)
    parser.add_argument("--n-samples", type=int, default=1_000_000)
    parser.add_argument("--x-min", type=float, default=-2295.0)
    parser.add_argument("--x-max", type=float, default=2295.0)
    parser.add_argument("--y-min", type=float, default=-1878.0)
    parser.add_argument("--y-max", type=float, default=1878.0)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--epochs", type=int, default=1000)
    parser.add_argument("--batch", type=int, default=512)
    parser.add_argument("--patience", type=int, default=50)
    parser.add_argument(
        "--mode",
        type=str,
        default="all",
        choices=["real", "syn", "mixed", "finetune", "all"],
    )
    parser.add_argument("--n-runs", type=int, default=1)
    parser.add_argument("--col-count", type=int, default=12)
    parser.add_argument("--sidelobes-threshold", type=float, default=0.84)
    parser.add_argument("--hidden-dims", type=str, default="1024,512,256,128,64")
    parser.add_argument(
        "--backend",
        type=str,
        default="tensorflow",
        choices=["tensorflow", "pytorch"],
    )
    parser.add_argument(
        "--arch",
        type=str,
        default="mlp",
        choices=["mlp", "res", "dualhead_y"],
        help="PyTorch architecture",
    )
    parser.add_argument("--loss", type=str, default="mse", choices=["mse", "huber"])
    parser.add_argument("--huber-delta", type=float, default=0.05)
    parser.add_argument(
        "--synthetic-data-path",
        type=str,
        default=None,
        help="Explicit NPZ path (otherwise derived from n_samples/x/y bounds)",
    )
    args = parser.parse_args()
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
        backend=args.backend,
        arch=args.arch,
        loss=args.loss,
        huber_delta=args.huber_delta,
        synthetic_data_path=args.synthetic_data_path,
    )


def main():
    print("[DEBUG] Starting stix-train", flush=True)
    config = parse_args()
    if config.backend == "pytorch":
        import torch

        print(f"[GPU][PyTorch] cuda={torch.cuda.is_available()}", flush=True)
    print(f"[DEBUG] Parsed config: {config}", flush=True)

    modes_to_train = (
        ["real", "syn", "mixed", "finetune"] if config.mode == "all" else [config.mode]
    )

    needs_synthetic = any(m in ["syn", "mixed", "finetune"] for m in modes_to_train)
    syn_data_path = None
    if needs_synthetic:
        if config.synthetic_data_path:
            syn_path = Path(config.synthetic_data_path)
            if not syn_path.is_absolute():
                syn_path = (PROJECT_ROOT / syn_path).resolve()
            if not syn_path.exists():
                raise FileNotFoundError(f"Synthetic data not found: {syn_path}")
            syn_data_path = str(syn_path)
            print(f"Using explicit synthetic data: {syn_data_path}")
        else:
            syn_data_path = ensure_synthetic_data(
                n_samples=config.n_samples,
                min_x=int(config.x_min),
                max_x=int(config.x_max),
                min_y=int(config.y_min),
                max_y=int(config.y_max),
            )

    save_experiment_config(config, syn_data_path)

    print("[DEBUG] Loading data...", flush=True)
    data = load_data(
        syn_data_path=syn_data_path,
        x_min=config.x_min,
        x_max=config.x_max,
        y_min=config.y_min,
        y_max=config.y_max,
        sidelobes_threshold=config.sidelobes_threshold,
        col_count=config.col_count,
    )
    print("[DEBUG] Data loaded", flush=True)

    results = {}
    syn_model_path = None
    for mode in modes_to_train:
        pretrained_path = syn_model_path if mode == "finetune" else None
        results[mode] = train_single_mode(mode, data, config, pretrained_path)
        if mode == "syn":
            syn_dir = EXPERIMENTS_DIR / config.experiment / "syn"
            suffix = ".pt" if config.backend == "pytorch" else ".keras"
            if config.n_runs > 1:
                syn_model_path = str(syn_dir / f"model_run0{suffix}")
            else:
                syn_model_path = str(syn_dir / f"model{suffix}")

    if len(modes_to_train) > 1:
        print("\n" + "=" * 60)
        print("FINAL COMPARISON (MAE on real test set)")
        print("=" * 60)
        for name, mae in results.items():
            if isinstance(mae, np.ndarray) and mae.ndim == 1:
                print(
                    f"{name:10s}: X={mae[0]:7.2f}, Y={mae[1]:7.2f}, Mean={np.mean(mae):7.2f}"
                )
            else:
                print(f"{name:10s}: {mae}")


if __name__ == "__main__":
    main()
