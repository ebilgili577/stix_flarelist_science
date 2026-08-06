import json
from datetime import datetime
from pathlib import Path
from typing import Any, List, Optional, Tuple

import numpy as np

from .config import TrainConfig, EXPERIMENTS_DIR, RANDOM_STATE
from .data import DataDict, denormalize_locations


def _get_torch_backend():
    from . import torch_backend

    return torch_backend


def evaluate_model(
    model,
    X_test: np.ndarray,
    y_test: np.ndarray,
    model_name: str,
) -> np.ndarray:
    preds = denormalize_locations(model.predict(X_test, verbose=0))
    y_true = denormalize_locations(y_test)
    mae = np.mean(np.abs(preds - y_true), axis=0)
    print(f"\n{model_name} - MAE on real test set:")
    print(f"  X: {mae[0]:.2f}, Y: {mae[1]:.2f}, Mean: {np.mean(mae):.2f}")
    return mae


def evaluate_model_torch(
    model,
    X_test: np.ndarray,
    y_test: np.ndarray,
    model_name: str,
) -> np.ndarray:
    torch_backend = _get_torch_backend()
    device = torch_backend.get_device()
    preds = denormalize_locations(torch_backend.predict(model, X_test, device))
    y_true = denormalize_locations(y_test)
    mae = np.mean(np.abs(preds - y_true), axis=0)
    print(f"\n{model_name} - MAE on real test set:")
    print(f"  X: {mae[0]:.2f}, Y: {mae[1]:.2f}, Mean: {np.mean(mae):.2f}")
    return mae


def _split_train_val(
    X: np.ndarray,
    y: np.ndarray,
    val_fraction: float = 0.1,
    random_seed: Optional[int] = None,
    real_head: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n = len(X)
    seed = RANDOM_STATE if random_seed is None else random_seed
    rng = np.random.default_rng(seed)
    if real_head is None or real_head >= n:
        n_val = max(1, min(n - 1, int(round(n * val_fraction))))
        indices = rng.permutation(n)
        val_idx, train_idx = indices[:n_val], indices[n_val:]
    else:
        if real_head < 2:
            raise ValueError(f"real_head={real_head}: need >= 2 real rows")
        n_val = max(1, min(real_head - 1, int(round(real_head * val_fraction))))
        real_indices = rng.permutation(real_head)
        val_idx = real_indices[:n_val]
        train_idx = np.concatenate([real_indices[n_val:], np.arange(real_head, n)])
        print(
            f"[VAL] Validation drawn from real rows only: "
            f"{n_val} real val / {real_head - n_val} real train / {n - real_head} synthetic train",
            flush=True,
        )
    return X[train_idx], y[train_idx], X[val_idx], y[val_idx]


def train_model(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    y_test: np.ndarray,
    model_name: str,
    learning_rate: float,
    epochs: int,
    batch_size: int,
    patience: int,
    save_dir: Path,
    hidden_dims: List[int],
    run_id: Optional[int] = None,
    random_seed: Optional[int] = None,
) -> Tuple[Any, np.ndarray]:
    from tensorflow import keras
    import tensorflow as tf
    from .model import build_model

    print(f"\n{'='*60}\nTraining {model_name}\n{'='*60}", flush=True)
    if random_seed is not None:
        tf.random.set_seed(random_seed)
        np.random.seed(random_seed)

    model = build_model(X_train.shape[1], hidden_dims, learning_rate)
    early_stop = keras.callbacks.EarlyStopping(
        monitor="val_loss", patience=patience, restore_best_weights=True
    )
    model.fit(
        X_train,
        y_train,
        validation_split=0.1,
        epochs=epochs,
        batch_size=batch_size,
        callbacks=[early_stop],
        verbose=0,
    )
    mae = evaluate_model(model, X_test, y_test, model_name)
    save_dir.mkdir(parents=True, exist_ok=True)
    save_path = (
        save_dir / f"model_run{run_id}.keras"
        if run_id is not None
        else save_dir / "model.keras"
    )
    model.save(save_path)
    print(f"Saved model to {save_path}")
    return model, mae


def train_model_pytorch(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    y_test: np.ndarray,
    model_name: str,
    learning_rate: float,
    epochs: int,
    batch_size: int,
    patience: int,
    save_dir: Path,
    hidden_dims: List[int],
    run_id: Optional[int] = None,
    random_seed: Optional[int] = None,
    config: Optional[TrainConfig] = None,
    val_real_head: Optional[int] = None,
    X_val_external: Optional[np.ndarray] = None,
    y_val_external: Optional[np.ndarray] = None,
) -> Tuple[Any, np.ndarray]:
    import torch

    torch_backend = _get_torch_backend()
    print(f"\n{'='*60}\nTraining {model_name} [PyTorch]\n{'='*60}", flush=True)
    if random_seed is not None:
        np.random.seed(random_seed)
        torch.manual_seed(random_seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(random_seed)

    if X_val_external is not None and y_val_external is not None:
        X_tr, y_tr = X_train, y_train
        X_val, y_val = X_val_external, y_val_external
        print(
            f"[VAL] External real validation: {len(y_val)} real val / {len(y_tr)} train",
            flush=True,
        )
    else:
        X_tr, y_tr, X_val, y_val = _split_train_val(
            X_train, y_train, random_seed=random_seed, real_head=val_real_head
        )

    device = torch_backend.get_device()
    print(f"[GPU][PyTorch] Device: {device}")
    arch = config.arch if config is not None else "mlp"
    loss = config.loss if config is not None else "mse"
    huber_delta = config.huber_delta if config is not None else 0.05

    model = torch_backend.build_model(
        input_shape=(X_train.shape[1],),
        hidden_dims=hidden_dims,
        arch=arch,
    )
    model = torch_backend.train_model(
        model=model,
        X_train=X_tr,
        y_train=y_tr,
        X_val=X_val,
        y_val=y_val,
        learning_rate=learning_rate,
        epochs=epochs,
        batch_size=batch_size,
        patience=patience,
        device=device,
        loss=loss,
        huber_delta=huber_delta,
    )
    mae = evaluate_model_torch(model, X_test, y_test, model_name)
    save_dir.mkdir(parents=True, exist_ok=True)
    save_path = (
        save_dir / f"model_run{run_id}.pt" if run_id is not None else save_dir / "model.pt"
    )
    torch_backend.save_model(
        save_path,
        model,
        input_shape=(X_train.shape[1],),
        hidden_dims=hidden_dims,
        arch=arch,
    )
    print(f"Saved model to {save_path}")
    return model, mae


def train_finetune(
    X_train_syn: np.ndarray,
    y_train_syn: np.ndarray,
    X_train_real: np.ndarray,
    y_train_real: np.ndarray,
    X_test: np.ndarray,
    y_test: np.ndarray,
    learning_rate: float,
    epochs: int,
    batch_size: int,
    patience: int,
    save_dir: Path,
    hidden_dims: List[int],
    run_id: Optional[int] = None,
    random_seed: Optional[int] = None,
    pretrained_model_path: Optional[str] = None,
) -> Tuple[Any, np.ndarray]:
    from tensorflow import keras
    import tensorflow as tf
    from .model import build_model

    print(f"\n{'='*60}\nFINETUNING MODE\n{'='*60}")
    if random_seed is not None:
        tf.random.set_seed(random_seed)
        np.random.seed(random_seed)

    if pretrained_model_path is not None and Path(pretrained_model_path).exists():
        print(f"\nStep 1: Loading pretrained model from {pretrained_model_path}")
        model = keras.models.load_model(pretrained_model_path)
    else:
        print("\nStep 1: Pretraining on synthetic data")
        model = build_model(X_train_syn.shape[1], hidden_dims, learning_rate)
        early_stop_pretrain = keras.callbacks.EarlyStopping(
            monitor="val_loss", patience=patience, restore_best_weights=True
        )
        model.fit(
            X_train_syn,
            y_train_syn,
            validation_split=0.1,
            epochs=epochs,
            batch_size=batch_size,
            callbacks=[early_stop_pretrain],
            verbose=0,
        )

    print("\nStep 2: Finetuning on real data")
    num_layers = len(model.layers)
    for i, layer in enumerate(model.layers):
        layer.trainable = i >= num_layers - 2
    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=learning_rate * 0.1),
        loss="mse",
        metrics=["mae"],
    )
    early_stop_finetune = keras.callbacks.EarlyStopping(
        monitor="val_loss", patience=patience, restore_best_weights=True
    )
    model.fit(
        X_train_real,
        y_train_real,
        validation_split=0.1,
        epochs=epochs,
        batch_size=batch_size,
        callbacks=[early_stop_finetune],
        verbose=0,
    )
    mae = evaluate_model(model, X_test, y_test, "finetune")
    save_dir.mkdir(parents=True, exist_ok=True)
    save_path = (
        save_dir / f"model_run{run_id}.keras"
        if run_id is not None
        else save_dir / "model.keras"
    )
    model.save(save_path)
    print(f"Saved model to {save_path}")
    return model, mae


def train_single_mode(
    mode: str,
    data: DataDict,
    config: TrainConfig,
    pretrained_syn_model_path: Optional[str] = None,
) -> np.ndarray:
    X_test, y_test = data["test_real"]
    models_dir = EXPERIMENTS_DIR / config.experiment / mode
    models_dir.mkdir(parents=True, exist_ok=True)
    all_maes = []

    for run_id in range(config.n_runs):
        random_seed = run_id if config.n_runs > 1 else None
        run_id_arg = run_id if config.n_runs > 1 else None
        use_torch = config.backend == "pytorch"

        if mode == "real":
            X_train, y_train = data["train_real"]
            if use_torch:
                _, mae = train_model_pytorch(
                    X_train,
                    y_train,
                    X_test,
                    y_test,
                    "train_on_real",
                    config.learning_rate,
                    config.epochs,
                    config.batch_size,
                    config.patience,
                    models_dir,
                    config.hidden_dims,
                    run_id_arg,
                    random_seed,
                    config,
                )
            else:
                _, mae = train_model(
                    X_train,
                    y_train,
                    X_test,
                    y_test,
                    "train_on_real",
                    config.learning_rate,
                    config.epochs,
                    config.batch_size,
                    config.patience,
                    models_dir,
                    config.hidden_dims,
                    run_id_arg,
                    random_seed,
                )
        elif mode == "syn":
            X_train, y_train = data["train_syn"]
            X_val_ext = y_val_ext = None
            if use_torch and config.val_on_real:
                X_val_ext, y_val_ext = data["train_real"]
            if use_torch:
                _, mae = train_model_pytorch(
                    X_train,
                    y_train,
                    X_test,
                    y_test,
                    "train_on_syn",
                    config.learning_rate,
                    config.epochs,
                    config.batch_size,
                    config.patience,
                    models_dir,
                    config.hidden_dims,
                    run_id_arg,
                    random_seed,
                    config,
                    X_val_external=X_val_ext,
                    y_val_external=y_val_ext,
                )
            else:
                _, mae = train_model(
                    X_train,
                    y_train,
                    X_test,
                    y_test,
                    "train_on_syn",
                    config.learning_rate,
                    config.epochs,
                    config.batch_size,
                    config.patience,
                    models_dir,
                    config.hidden_dims,
                    run_id_arg,
                    random_seed,
                )
        elif mode == "mixed":
            X_train, y_train = data["train_mixed"]
            val_real_head = len(data["train_real"][0]) if config.val_on_real else None
            if use_torch:
                _, mae = train_model_pytorch(
                    X_train,
                    y_train,
                    X_test,
                    y_test,
                    "train_on_mixed",
                    config.learning_rate,
                    config.epochs,
                    config.batch_size,
                    config.patience,
                    models_dir,
                    config.hidden_dims,
                    run_id_arg,
                    random_seed,
                    config,
                    val_real_head=val_real_head,
                )
            else:
                _, mae = train_model(
                    X_train,
                    y_train,
                    X_test,
                    y_test,
                    "train_on_mixed",
                    config.learning_rate,
                    config.epochs,
                    config.batch_size,
                    config.patience,
                    models_dir,
                    config.hidden_dims,
                    run_id_arg,
                    random_seed,
                )
        elif mode == "finetune":
            if use_torch:
                raise NotImplementedError(
                    "finetune with --backend pytorch is not implemented yet; use tensorflow or mode=mixed/syn"
                )
            X_train_syn, y_train_syn = data["train_syn"]
            X_train_real, y_train_real = data["train_real"]
            syn_model_path = None
            if pretrained_syn_model_path is not None:
                if config.n_runs > 1:
                    syn_model_path = str(
                        Path(pretrained_syn_model_path).parent / f"model_run{run_id}.keras"
                    )
                else:
                    syn_model_path = pretrained_syn_model_path
            _, mae = train_finetune(
                X_train_syn,
                y_train_syn,
                X_train_real,
                y_train_real,
                X_test,
                y_test,
                config.learning_rate,
                config.epochs,
                config.batch_size,
                config.patience,
                models_dir,
                config.hidden_dims,
                run_id_arg,
                random_seed,
                syn_model_path,
            )
        else:
            raise ValueError(f"Unknown mode: {mode}")
        all_maes.append(mae)

    if config.n_runs > 1:
        mae_array = np.array(all_maes)
        mae_mean = np.mean(mae_array, axis=0)
        mae_std = np.std(mae_array, axis=0)
        print(f"\n{'='*60}")
        print(f"STATISTICS FOR {mode.upper()} MODE ({config.n_runs} runs)")
        print(f"  X: {mae_mean[0]:.2f} ± {mae_std[0]:.2f}")
        print(f"  Y: {mae_mean[1]:.2f} ± {mae_std[1]:.2f}")
        print(f"{'='*60}")
        return mae_mean
    return all_maes[0]


def save_experiment_config(config: TrainConfig, syn_data_path: Optional[str] = None) -> None:
    exp_dir = EXPERIMENTS_DIR / config.experiment
    exp_dir.mkdir(parents=True, exist_ok=True)
    config_path = exp_dir / "config.json"
    config_dict = config.to_dict()
    config_dict["timestamp"] = datetime.now().isoformat()
    config_dict["synthetic_data_path"] = syn_data_path
    with open(config_path, "w") as f:
        json.dump(config_dict, f, indent=2)
    print(f"\nConfiguration saved to {config_path}")
