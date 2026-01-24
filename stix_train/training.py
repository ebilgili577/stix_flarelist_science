import json
from datetime import datetime
from pathlib import Path
from typing import Optional, Tuple, List, Dict

import numpy as np
from tensorflow import keras
from sklearn.metrics import r2_score

from .config import TrainConfig, EXPERIMENTS_DIR
from .model import build_model
from .data import DataDict, denormalize_locations


def compute_extended_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    save_dir: Path,
    model_name: str
) -> Dict:
    """
    Compute extended metrics including q95, error distribution, and spatial error map.
    
    Args:
        y_true: True values (denormalized, in arcsec)
        y_pred: Predicted values (denormalized, in arcsec)
        save_dir: Directory to save metrics and plots
        model_name: Name for file naming
        
    Returns:
        Dictionary with all metrics
    """
    errors = y_pred - y_true
    abs_errors = np.abs(errors)
    
    # Basic metrics
    mae_x = np.mean(abs_errors[:, 0])
    mae_y = np.mean(abs_errors[:, 1])
    median_ae_x = np.median(abs_errors[:, 0])
    median_ae_y = np.median(abs_errors[:, 1])
    p95_x = np.percentile(abs_errors[:, 0], 95)
    p95_y = np.percentile(abs_errors[:, 1], 95)
    r2_x = r2_score(y_true[:, 0], y_pred[:, 0])
    r2_y = r2_score(y_true[:, 1], y_pred[:, 1])
    
    # Euclidean distance error
    euclidean = np.sqrt(errors[:, 0]**2 + errors[:, 1]**2)
    euclidean_mae = np.mean(euclidean)
    euclidean_median = np.median(euclidean)
    euclidean_p95 = np.percentile(euclidean, 95)
    
    metrics = {
        'mae_x': float(mae_x),
        'mae_y': float(mae_y),
        'mae_mean': float((mae_x + mae_y) / 2),
        'median_ae_x': float(median_ae_x),
        'median_ae_y': float(median_ae_y),
        'p95_x': float(p95_x),
        'p95_y': float(p95_y),
        'r2_x': float(r2_x),
        'r2_y': float(r2_y),
        'euclidean_mae': float(euclidean_mae),
        'euclidean_median': float(euclidean_median),
        'euclidean_p95': float(euclidean_p95),
    }
    
    # Save error distribution data
    error_data = {
        'errors_x': errors[:, 0],
        'errors_y': errors[:, 1],
        'abs_errors_x': abs_errors[:, 0],
        'abs_errors_y': abs_errors[:, 1],
        'euclidean_errors': euclidean,
        'true_xy': y_true,
        'pred_xy': y_pred,
    }
    
    save_dir.mkdir(parents=True, exist_ok=True)
    error_file = save_dir / f"{model_name}_errors.npz"
    np.savez_compressed(error_file, **error_data)
    print(f"Saved error distribution data to {error_file}")
    
    # Save spatial error map data (for plotting)
    spatial_data = {
        'true_x': y_true[:, 0],
        'true_y': y_true[:, 1],
        'pred_x': y_pred[:, 0],
        'pred_y': y_pred[:, 1],
        'euclidean_error': euclidean,
    }
    spatial_file = save_dir / f"{model_name}_spatial_errors.npz"
    np.savez_compressed(spatial_file, **spatial_data)
    print(f"Saved spatial error map data to {spatial_file}")
    
    # Save metrics JSON
    metrics_file = save_dir / f"{model_name}_extended_metrics.json"
    with open(metrics_file, 'w') as f:
        json.dump(metrics, f, indent=2)
    print(f"Saved extended metrics to {metrics_file}")
    
    return metrics


def evaluate_model(
    model: keras.Model,
    X_test: np.ndarray,
    y_test: np.ndarray,
    model_name: str,
    extended_metrics: bool = False,
    save_dir: Optional[Path] = None,
    metrics_name: Optional[str] = None
) -> np.ndarray:
    """
    Evaluate model on test set and print MAE.
    
    Args:
        model: Trained Keras model
        X_test: Test features
        y_test: Test targets (normalized)
        model_name: Name for logging
        extended_metrics: If True, compute and save extended metrics
        save_dir: Directory to save extended metrics (required if extended_metrics=True)
        
    Returns:
        MAE array [x_mae, y_mae]
    """
    preds = denormalize_locations(model.predict(X_test, verbose=0))
    y_true = denormalize_locations(y_test)
    mae = np.mean(np.abs(preds - y_true), axis=0)
    
    print(f"\n{model_name} - MAE on real test set:", flush=True)
    print(f"  X: {mae[0]:.2f}, Y: {mae[1]:.2f}, Mean: {np.mean(mae):.2f}", flush=True)
    
    if extended_metrics and save_dir is not None:
        metrics_file_name = metrics_name if metrics_name is not None else model_name
        print(f"\nComputing extended metrics for {model_name}...", flush=True)
        extended = compute_extended_metrics(y_true, preds, save_dir, metrics_file_name)
        print(f"  MAE:       X={extended['mae_x']:7.2f}, Y={extended['mae_y']:7.2f}, Mean={extended['mae_mean']:7.2f}", flush=True)
        print(f"  Median AE: X={extended['median_ae_x']:7.2f}, Y={extended['median_ae_y']:7.2f}", flush=True)
        print(f"  P95:       X={extended['p95_x']:7.2f}, Y={extended['p95_y']:7.2f}", flush=True)
        print(f"  R²:        X={extended['r2_x']:7.4f}, Y={extended['r2_y']:7.4f}", flush=True)
        print(f"  Euclidean: MAE={extended['euclidean_mae']:7.2f}, Median={extended['euclidean_median']:7.2f}, P95={extended['euclidean_p95']:7.2f}", flush=True)
    
    return mae


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
    extended_metrics: bool = False,
) -> Tuple[keras.Model, np.ndarray]:
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
        run_id: Optional run ID for multiple runs
        random_seed: Optional random seed
        
    Returns:
        Tuple of (trained model, MAE array)
    """
    import tensorflow as tf
    
    print(f"\n{'='*60}", flush=True)
    print(f"Training {model_name}", flush=True)
    if run_id is not None:
        print(f"Run {run_id}", flush=True)
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
    
    print(f"Starting training: {X_train.shape[0]} samples, {epochs} max epochs, batch_size={batch_size}", flush=True)
    print("Training in progress... (this may take a while)", flush=True)
    
    history = model.fit(
        X_train,
        y_train,
        validation_split=0.1,
        epochs=epochs,
        batch_size=batch_size,
        callbacks=[early_stop],
        verbose=1,  # Show progress
    )
    
    epochs_trained = len(history.history['loss'])
    print(f"Training completed: {epochs_trained} epochs trained", flush=True)
    
    # Determine metrics name for extended metrics
    if run_id is not None:
        metrics_name = f"{model_name}_run{run_id}"
    else:
        metrics_name = model_name
    
    mae = evaluate_model(model, X_test, y_test, model_name, extended_metrics, save_dir, metrics_name)
    
    # Save model
    save_dir.mkdir(parents=True, exist_ok=True)
    if run_id is not None:
        save_path = save_dir / f"model_run{run_id}.keras"
    else:
        save_path = save_dir / "model.keras"
    model.save(save_path)
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
    extended_metrics: bool = False,
) -> Tuple[keras.Model, np.ndarray]:
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
        run_id: Optional run ID
        random_seed: Optional random seed
        pretrained_model_path: Optional path to pretrained model
        
    Returns:
        Tuple of (trained model, MAE array)
    """
    import tensorflow as tf
    
    print(f"\n{'='*60}", flush=True)
    print("FINETUNING MODE", flush=True)
    if run_id is not None:
        print(f"Run {run_id}", flush=True)
    print(f"{'='*60}", flush=True)
    
    # Set random seed if provided
    if random_seed is not None:
        tf.random.set_seed(random_seed)
        np.random.seed(random_seed)
    
    # Step 1: Load pretrained model or train on synthetic data
    if pretrained_model_path is not None and Path(pretrained_model_path).exists():
        print(f"\nStep 1: Loading pretrained model from {pretrained_model_path}", flush=True)
        model = keras.models.load_model(pretrained_model_path)
    else:
        print("\nStep 1: Pretraining on synthetic data", flush=True)
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
        
        print(f"Pretraining: {X_train_syn.shape[0]} samples, {epochs} max epochs", flush=True)
        history_pretrain = model.fit(
            X_train_syn,
            y_train_syn,
            validation_split=0.1,
            epochs=epochs,
            batch_size=batch_size,
            callbacks=[early_stop_pretrain],
            verbose=1,  # Show progress
        )
        print(f"Pretraining completed: {len(history_pretrain.history['loss'])} epochs", flush=True)
    
    # Step 2: Finetune on real data (freeze all except last 2 layers)
    print("\nStep 2: Finetuning on real data", flush=True)
    
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
    
    print(f"Finetuning: {X_train_real.shape[0]} samples, {epochs} max epochs", flush=True)
    history_finetune = model.fit(
        X_train_real,
        y_train_real,
        validation_split=0.1,
        epochs=epochs,
        batch_size=batch_size,
        callbacks=[early_stop_finetune],
        verbose=1,  # Show progress
    )
    print(f"Finetuning completed: {len(history_finetune.history['loss'])} epochs", flush=True)
    
    # Determine metrics name for extended metrics
    if run_id is not None:
        metrics_name = f"finetune_run{run_id}"
    else:
        metrics_name = "finetune"
    
    mae = evaluate_model(model, X_test, y_test, "finetune", extended_metrics, save_dir, metrics_name)
    
    # Save model
    save_dir.mkdir(parents=True, exist_ok=True)
    if run_id is not None:
        save_path = save_dir / f"model_run{run_id}.keras"
    else:
        save_path = save_dir / "model.keras"
    model.save(save_path)
    print(f"Saved model to {save_path}")
    
    return model, mae


def train_single_mode(
    mode: str,
    data: DataDict,
    config: TrainConfig,
    pretrained_syn_model_path: Optional[str] = None,
    extended_metrics: bool = False,
) -> np.ndarray:
    """
    Train a single mode (real, syn, mixed, or finetune) with optional multiple runs.
    
    Args:
        mode: Training mode ('real', 'syn', 'mixed', 'finetune')
        data: Dictionary with train/test data
        config: Training configuration
        pretrained_syn_model_path: Optional path to pretrained syn model for finetune
        
    Returns:
        MAE array (or mean MAE if multiple runs)
    """
    
    X_test, y_test = data["test_real"]
    
    # Create model directory
    models_dir = EXPERIMENTS_DIR / config.experiment / mode
    models_dir.mkdir(parents=True, exist_ok=True)
    
    all_maes = []
    
    for run_id in range(config.n_runs):
        random_seed = run_id if config.n_runs > 1 else None
        run_id_arg = run_id if config.n_runs > 1 else None
        
        if mode == "real":
            X_train, y_train = data["train_real"]
            model, mae = train_model(
                X_train, y_train, X_test, y_test,
                "train_on_real", config.learning_rate, config.epochs,
                config.batch_size, config.patience, models_dir,
                config.hidden_dims, run_id_arg, random_seed, extended_metrics
            )
        elif mode == "syn":
            X_train, y_train = data["train_syn"]
            model, mae = train_model(
                X_train, y_train, X_test, y_test,
                "train_on_syn", config.learning_rate, config.epochs,
                config.batch_size, config.patience, models_dir,
                config.hidden_dims, run_id_arg, random_seed, extended_metrics
            )
        elif mode == "mixed":
            X_train, y_train = data["train_mixed"]
            model, mae = train_model(
                X_train, y_train, X_test, y_test,
                "train_on_mixed", config.learning_rate, config.epochs,
                config.batch_size, config.patience, models_dir,
                config.hidden_dims, run_id_arg, random_seed, extended_metrics
            )
        elif mode == "finetune":
            X_train_syn, y_train_syn = data["train_syn"]
            X_train_real, y_train_real = data["train_real"]
            
            # Determine pretrained model path for this run
            syn_model_path = None
            if pretrained_syn_model_path is not None:
                if config.n_runs > 1:
                    syn_model_path = str(
                        Path(pretrained_syn_model_path).parent / f"model_run{run_id}.keras"
                    )
                else:
                    syn_model_path = pretrained_syn_model_path
            
            model, mae = train_finetune(
                X_train_syn, y_train_syn, X_train_real, y_train_real,
                X_test, y_test, config.learning_rate, config.epochs,
                config.batch_size, config.patience, models_dir,
                config.hidden_dims, run_id_arg, random_seed, syn_model_path, extended_metrics
            )
        else:
            raise ValueError(f"Unknown mode: {mode}")
        
        all_maes.append(mae)
    
    # Compute statistics if multiple runs
    if config.n_runs > 1:
        mae_array = np.array(all_maes)
        mae_mean = np.mean(mae_array, axis=0)
        mae_std = np.std(mae_array, axis=0)
        
        print(f"\n{'='*60}", flush=True)
        print(f"STATISTICS FOR {mode.upper()} MODE ({config.n_runs} runs)", flush=True)
        print(f"{'='*60}", flush=True)
        print(f"MAE (mean ± std):", flush=True)
        print(f"  X: {mae_mean[0]:.2f} ± {mae_std[0]:.2f}", flush=True)
        print(f"  Y: {mae_mean[1]:.2f} ± {mae_std[1]:.2f}", flush=True)
        print(f"  Mean: {np.mean(mae_mean):.2f} ± {np.mean(mae_std):.2f}", flush=True)
        print(f"{'='*60}", flush=True)
        
        return mae_mean
    
    return all_maes[0]


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
