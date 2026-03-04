import json
from datetime import datetime
from pathlib import Path
from typing import Optional, Tuple, List

import numpy as np
from tensorflow import keras

from .config import TrainConfig, EXPERIMENTS_DIR
from .model import build_model
from .data import DataDict, denormalize_locations


def compute_extended_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    save_dir: Path,
    metrics_name: str,
) -> dict:
    """
    Compute comprehensive evaluation metrics and save to disk.
    
    Saves three files:
      - {metrics_name}_extended_metrics.json  (scalar metrics)
      - {metrics_name}_errors.npz             (error distributions)
      - {metrics_name}_spatial_errors.npz     (per-sample spatial data)
    
    Args:
        y_true: Ground-truth locations (N, 2), denormalized
        y_pred: Predicted locations (N, 2), denormalized
        save_dir: Directory to save output files
        metrics_name: Base name for output files
        
    Returns:
        Dictionary of scalar metrics
    """
    errors = y_pred - y_true
    abs_err_x = np.abs(errors[:, 0])
    abs_err_y = np.abs(errors[:, 1])
    eucl_dist = np.sqrt(np.sum(errors**2, axis=1))

    metrics = {
        "mae_x": float(np.mean(abs_err_x)),
        "mae_y": float(np.mean(abs_err_y)),
        "mae_mean": float(np.mean(np.abs(errors))),
        "median_ae_x": float(np.median(abs_err_x)),
        "median_ae_y": float(np.median(abs_err_y)),
        "median_ae_mean": float(np.median(np.abs(errors))),
        "mse_x": float(np.mean(errors[:, 0] ** 2)),
        "mse_y": float(np.mean(errors[:, 1] ** 2)),
        "rmse_x": float(np.sqrt(np.mean(errors[:, 0] ** 2))),
        "rmse_y": float(np.sqrt(np.mean(errors[:, 1] ** 2))),
        "q90_x": float(np.percentile(abs_err_x, 90)),
        "q90_y": float(np.percentile(abs_err_y, 90)),
        "q95_x": float(np.percentile(abs_err_x, 95)),
        "q95_y": float(np.percentile(abs_err_y, 95)),
        "q99_x": float(np.percentile(abs_err_x, 99)),
        "q99_y": float(np.percentile(abs_err_y, 99)),
        "max_x": float(np.max(abs_err_x)),
        "max_y": float(np.max(abs_err_y)),
        "euclidean_mean": float(np.mean(eucl_dist)),
        "euclidean_median": float(np.median(eucl_dist)),
        "euclidean_std": float(np.std(eucl_dist)),
        "euclidean_q50": float(np.percentile(eucl_dist, 50)),
        "euclidean_q90": float(np.percentile(eucl_dist, 90)),
        "euclidean_q95": float(np.percentile(eucl_dist, 95)),
        "euclidean_q99": float(np.percentile(eucl_dist, 99)),
        "euclidean_max": float(np.max(eucl_dist)),
        "n_test_samples": len(y_true),
    }

    save_dir.mkdir(parents=True, exist_ok=True)

    with open(save_dir / f"{metrics_name}_extended_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)

    np.savez(
        save_dir / f"{metrics_name}_errors.npz",
        euclidean_errors=eucl_dist,
    )

    np.savez(
        save_dir / f"{metrics_name}_spatial_errors.npz",
        true_x=y_true[:, 0],
        true_y=y_true[:, 1],
        pred_x=y_pred[:, 0],
        pred_y=y_pred[:, 1],
        euclidean_error=eucl_dist,
    )

    print(f"  Extended metrics saved to {save_dir / metrics_name}_*.json/npz")
    return metrics


def evaluate_model(
    model: keras.Model,
    X_test: np.ndarray,
    y_test: np.ndarray,
    model_name: str,
    save_dir: Optional[Path] = None,
    metrics_name: Optional[str] = None,
) -> dict:
    """
    Evaluate model on test set, compute extended metrics, and optionally save.
    
    Args:
        model: Trained Keras model
        X_test: Test features
        y_test: Test targets (normalized)
        model_name: Name for logging
        save_dir: Directory to save metrics (if None, metrics are not saved)
        metrics_name: Base filename for saved metrics
        
    Returns:
        Dictionary of extended metrics
    """
    y_pred = denormalize_locations(model.predict(X_test, verbose=0))
    y_true = denormalize_locations(y_test)

    if save_dir is not None and metrics_name is not None:
        metrics = compute_extended_metrics(y_true, y_pred, save_dir, metrics_name)
    else:
        errors = y_pred - y_true
        eucl_dist = np.sqrt(np.sum(errors**2, axis=1))
        metrics = {
            "mae_x": float(np.mean(np.abs(errors[:, 0]))),
            "mae_y": float(np.mean(np.abs(errors[:, 1]))),
            "mae_mean": float(np.mean(np.abs(errors))),
            "euclidean_mean": float(np.mean(eucl_dist)),
        }

    print(f"\n{model_name} - Results on real test set:")
    print(f"  MAE  X: {metrics['mae_x']:.2f}, Y: {metrics['mae_y']:.2f}, Mean: {metrics['mae_mean']:.2f}")
    print(f"  Eucl Mean: {metrics.get('euclidean_mean', 0):.2f}, "
          f"Q95: {metrics.get('euclidean_q95', 0):.2f}")

    return metrics


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
) -> Tuple[keras.Model, dict]:
    """
    Train a single model and return it along with evaluation metrics.
    
    Args:
        X_train: Training features
        y_train: Training targets
        X_test: Test features
        y_test: Test targets
        model_name: Name for the model
        learning_rate: Learning rate
        epochs: Maximum epochs
        batch_size: Batch size
        patience: Early stopping patience
        save_dir: Directory to save model
        hidden_dims: Hidden layer sizes
        run_id: Optional run ID for multiple runs
        random_seed: Optional random seed
        
    Returns:
        Tuple of (trained model, metrics dict)
    """
    import tensorflow as tf
    
    print(f"\n{'='*60}")
    print(f"Training {model_name}")
    if run_id is not None:
        print(f"Run {run_id}")
    print(f"{'='*60}", flush=True)
    
    # Set random seed if provided
    if random_seed is not None:
        tf.random.set_seed(random_seed)
        np.random.seed(random_seed)
    
    model = build_model(
        input_dim=X_train.shape[1],
        hidden_dims=hidden_dims,
        learning_rate=learning_rate
    )
    
    early_stop = keras.callbacks.EarlyStopping(
        monitor="val_loss",
        patience=patience,
        restore_best_weights=True,
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
    
    # Build metrics filename matching plot_results.py expectations
    if run_id is not None:
        metrics_name = f"{model_name}_run{run_id}"
    else:
        metrics_name = model_name

    metrics = evaluate_model(model, X_test, y_test, model_name, save_dir, metrics_name)
    
    # Save model
    save_dir.mkdir(parents=True, exist_ok=True)
    if run_id is not None:
        save_path = save_dir / f"model_run{run_id}.keras"
    else:
        save_path = save_dir / "model.keras"
    model.save(save_path)
    print(f"Saved model to {save_path}")
    
    return model, metrics


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
) -> Tuple[keras.Model, dict]:
    """
    Pretrain on synthetic data, then finetune on real data.
    
    Args:
        X_train_syn: Synthetic training features
        y_train_syn: Synthetic training targets
        X_train_real: Real training features
        y_train_real: Real training targets
        X_test: Test features
        y_test: Test targets
        learning_rate: Learning rate
        epochs: Maximum epochs
        batch_size: Batch size
        patience: Early stopping patience
        save_dir: Directory to save model
        hidden_dims: Hidden layer sizes
        run_id: Optional run ID
        random_seed: Optional random seed
        pretrained_model_path: Optional path to pretrained model
        
    Returns:
        Tuple of (trained model, metrics dict)
    """
    import tensorflow as tf
    
    print(f"\n{'='*60}")
    print("FINETUNING MODE")
    if run_id is not None:
        print(f"Run {run_id}")
    print(f"{'='*60}")
    
    # Set random seed if provided
    if random_seed is not None:
        tf.random.set_seed(random_seed)
        np.random.seed(random_seed)
    
    # Step 1: Load pretrained model or train on synthetic data
    if pretrained_model_path is not None and Path(pretrained_model_path).exists():
        print(f"\nStep 1: Loading pretrained model from {pretrained_model_path}")
        model = keras.models.load_model(pretrained_model_path)
    else:
        print("\nStep 1: Pretraining on synthetic data")
        model = build_model(
            input_dim=X_train_syn.shape[1],
            hidden_dims=hidden_dims,
            learning_rate=learning_rate
        )
        
        early_stop_pretrain = keras.callbacks.EarlyStopping(
            monitor="val_loss",
            patience=patience,
            restore_best_weights=True,
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
    
    # Step 2: Finetune on real data (freeze all except last 2 layers)
    print("\nStep 2: Finetuning on real data")
    
    num_layers = len(model.layers)
    for i, layer in enumerate(model.layers):
        layer.trainable = i >= num_layers - 2
    
    # Use lower learning rate for finetuning
    finetune_lr = learning_rate * 0.1
    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=finetune_lr),
        loss="mse",
        metrics=["mae"],
    )
    
    early_stop_finetune = keras.callbacks.EarlyStopping(
        monitor="val_loss",
        patience=patience,
        restore_best_weights=True,
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
    
    if run_id is not None:
        metrics_name = f"finetune_run{run_id}"
    else:
        metrics_name = "finetune"

    metrics = evaluate_model(model, X_test, y_test, "finetune", save_dir, metrics_name)
    
    # Save model
    save_dir.mkdir(parents=True, exist_ok=True)
    if run_id is not None:
        save_path = save_dir / f"model_run{run_id}.keras"
    else:
        save_path = save_dir / "model.keras"
    model.save(save_path)
    print(f"Saved model to {save_path}")
    
    return model, metrics


def train_single_mode(
    mode: str,
    data: DataDict,
    config: TrainConfig,
    pretrained_syn_model_path: Optional[str] = None,
) -> dict:
    """
    Train a single mode (real, syn, mixed, or finetune) with optional multiple runs.
    
    Args:
        mode: Training mode ('real', 'syn', 'mixed', 'finetune')
        data: Dictionary with train/test data
        config: Training configuration
        pretrained_syn_model_path: Optional path to pretrained syn model for finetune
        
    Returns:
        Dictionary of extended metrics (averaged over runs if n_runs > 1)
    """
    
    X_test, y_test = data["test_real"]
    
    # Create model directory
    models_dir = EXPERIMENTS_DIR / config.experiment / mode
    models_dir.mkdir(parents=True, exist_ok=True)
    
    all_metrics: List[dict] = []
    
    for run_id in range(config.n_runs):
        random_seed = run_id if config.n_runs > 1 else None
        run_id_arg = run_id if config.n_runs > 1 else None
        
        if mode == "real":
            X_train, y_train = data["train_real"]
            model, metrics = train_model(
                X_train, y_train, X_test, y_test,
                "train_on_real", config.learning_rate, config.epochs,
                config.batch_size, config.patience, models_dir,
                config.hidden_dims, run_id_arg, random_seed
            )
        elif mode == "syn":
            X_train, y_train = data["train_syn"]
            model, metrics = train_model(
                X_train, y_train, X_test, y_test,
                "train_on_syn", config.learning_rate, config.epochs,
                config.batch_size, config.patience, models_dir,
                config.hidden_dims, run_id_arg, random_seed
            )
        elif mode == "mixed":
            X_train, y_train = data["train_mixed"]
            model, metrics = train_model(
                X_train, y_train, X_test, y_test,
                "train_on_mixed", config.learning_rate, config.epochs,
                config.batch_size, config.patience, models_dir,
                config.hidden_dims, run_id_arg, random_seed
            )
        elif mode == "finetune":
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
            
            model, metrics = train_finetune(
                X_train_syn, y_train_syn, X_train_real, y_train_real,
                X_test, y_test, config.learning_rate, config.epochs,
                config.batch_size, config.patience, models_dir,
                config.hidden_dims, run_id_arg, random_seed, syn_model_path
            )
        else:
            raise ValueError(f"Unknown mode: {mode}")
        
        all_metrics.append(metrics)
    
    if config.n_runs > 1:
        numeric_keys = [k for k in all_metrics[0] if isinstance(all_metrics[0][k], (int, float))]
        mean_metrics = {}
        for key in numeric_keys:
            values = [m[key] for m in all_metrics]
            mean_metrics[key] = float(np.mean(values))
            mean_metrics[f"{key}_std"] = float(np.std(values))

        print(f"\n{'='*60}")
        print(f"STATISTICS FOR {mode.upper()} MODE ({config.n_runs} runs)")
        print(f"{'='*60}")
        print(f"MAE (mean +/- std):")
        print(f"  X: {mean_metrics['mae_x']:.2f} +/- {mean_metrics['mae_x_std']:.2f}")
        print(f"  Y: {mean_metrics['mae_y']:.2f} +/- {mean_metrics['mae_y_std']:.2f}")
        print(f"  Mean: {mean_metrics['mae_mean']:.2f} +/- {mean_metrics['mae_mean_std']:.2f}")
        print(f"Euclidean distance (mean +/- std):")
        print(f"  Mean: {mean_metrics['euclidean_mean']:.2f} +/- {mean_metrics['euclidean_mean_std']:.2f}")
        print(f"  Q95:  {mean_metrics['euclidean_q95']:.2f} +/- {mean_metrics['euclidean_q95_std']:.2f}")
        print(f"{'='*60}")

        return mean_metrics
    
    return all_metrics[0]


def save_experiment_config(config: TrainConfig, syn_data_path: Optional[str] = None) -> None:
    """
    Save experiment configuration to JSON file.
    
    Args:
        config: Training configuration
        syn_data_path: Path to synthetic data (for logging)
    """
    
    exp_dir = EXPERIMENTS_DIR / config.experiment
    exp_dir.mkdir(parents=True, exist_ok=True)
    
    config_path = exp_dir / "config.json"
    
    config_dict = config.to_dict()
    config_dict["timestamp"] = datetime.now().isoformat()
    config_dict["synthetic_data_path"] = syn_data_path
    
    with open(config_path, "w") as f:
        json.dump(config_dict, f, indent=2)
    
    print(f"\nConfiguration saved to {config_path}")
