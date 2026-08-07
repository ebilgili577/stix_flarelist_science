#!/usr/bin/env python3
"""Interactive SLURM submission for `stix-train` (`stix-submit`)."""

import subprocess
from pathlib import Path

from .config import TrainConfig, DEFAULT_HIDDEN_DIMS, EXPERIMENTS_DIR


def ask(prompt: str, default=None) -> str:
    suffix = f" [{default}]" if default is not None else ""
    response = input(f"{prompt}{suffix}: ").strip()
    return response if response else str(default) if default is not None else ""


def main():
    package_root = Path(__file__).resolve().parent.parent
    exp = ask("experiment name", "exp1")
    mode = ask("mode (real/syn/mixed/finetune/all)", "all")
    col_count = ask("how many subcols?", TrainConfig.col_count)
    sidelobes_threshold = ask(
        "sidelobes threshold (1 for no filtering)", TrainConfig.sidelobes_threshold
    )
    n_samples = (
        ask("synthetic n_samples", TrainConfig.n_samples)
        if mode != "real"
        else str(TrainConfig.n_samples)
    )
    x_min = ask("x_min", TrainConfig.x_min)
    x_max = ask("x_max", TrainConfig.x_max)
    y_min = ask("y_min", TrainConfig.y_min)
    y_max = ask("y_max", TrainConfig.y_max)
    lr = ask("learning rate", TrainConfig.learning_rate)
    epochs = ask("epochs", TrainConfig.epochs)
    batch = ask("batch size", TrainConfig.batch_size)
    patience = ask("early stop patience", TrainConfig.patience)
    n_runs = ask("n_runs", TrainConfig.n_runs)
    hidden_dims = ask(
        "hidden layer sizes (comma-separated)",
        ",".join(map(str, DEFAULT_HIDDEN_DIMS)),
    )
    backend = ask("backend (tensorflow/pytorch)", TrainConfig.backend)
    arch = TrainConfig.arch
    loss = TrainConfig.loss
    huber_delta = TrainConfig.huber_delta
    synthetic_data_path = ""
    if backend == "pytorch":
        arch = ask("arch (mlp/res/dualhead_y)", TrainConfig.arch)
        loss = ask("loss (mse/huber)", TrainConfig.loss)
        if loss == "huber":
            huber_delta = ask("huber delta", TrainConfig.huber_delta)
    if mode != "real":
        synthetic_data_path = ask("synthetic NPZ path (empty = auto)", "")

    description = ask("description of experiment", "")
    log_dir = EXPERIMENTS_DIR / exp / mode
    log_dir.mkdir(parents=True, exist_ok=True)
    if description:
        desc_path = EXPERIMENTS_DIR / exp / "description.txt"
        desc_path.parent.mkdir(parents=True, exist_ok=True)
        desc_path.write_text(description, encoding="utf-8")
        print(f"Description saved to {desc_path}")

    cmd = [
        "sbatch",
        "-J",
        exp,
        "-o",
        str(EXPERIMENTS_DIR / exp / mode / "train_%j.log"),
        str(package_root / "slurm" / "send_job.sh"),
        "--experiment",
        exp,
        "--n-samples",
        str(n_samples),
        "--x-min",
        str(x_min),
        "--x-max",
        str(x_max),
        "--y-min",
        str(y_min),
        "--y-max",
        str(y_max),
        "--lr",
        str(lr),
        "--epochs",
        str(epochs),
        "--batch",
        str(batch),
        "--patience",
        str(patience),
        "--mode",
        mode,
        "--n-runs",
        str(n_runs),
        "--col-count",
        str(col_count),
        "--sidelobes-threshold",
        str(sidelobes_threshold),
        "--hidden-dims",
        hidden_dims,
        "--backend",
        backend,
        "--arch",
        str(arch),
        "--loss",
        str(loss),
        "--huber-delta",
        str(huber_delta),
    ]
    if synthetic_data_path:
        cmd.extend(["--synthetic-data-path", synthetic_data_path])

    print("\nSubmitting:", " ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
