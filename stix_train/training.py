import json
from datetime import datetime
from pathlib import Path
from typing import Optional, Tuple, List

import numpy as np
from tensorflow import keras

from .config import TrainConfig, EXPERIMENTS_DIR
from .model import build_model
from .data import DataDict, denormalize_locations


def evaluate_model(
    model: keras.Model,
    X_test: np.ndarray,
    y_test: np.ndarray,
    model_name: str
) -> np.ndarray:
    """
    Evaluate model on test set and print MAE.
    
    Args:
        model: Trained Keras model
        X_test: Test features
        y_test: Test targets (normalized)
        model_name: Name for logging
        
    Returns:
        MAE array [x_mae, y_mae]
    """
    preds = denormalize_locations(model.predict(X_test, verbose=0))
    y_true = denormalize_locations(y_test)
    mae = np.mean(np.abs(preds - y_true), axis=0)
    
    print(f"\n{model_name} - MAE on real test set:")
    print(f"  X: {mae[0]:.2f}, Y: {mae[1]:.2f}, Mean: {np.mean(mae):.2f}")
    
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
    
    mae = evaluate_model(model, X_test, y_test, model_name)
    
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
    
    mae = evaluate_model(model, X_test, y_test, "finetune")
    
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
                config.hidden_dims, run_id_arg, random_seed
            )
        elif mode == "syn":
            X_train, y_train = data["train_syn"]
            model, mae = train_model(
                X_train, y_train, X_test, y_test,
                "train_on_syn", config.learning_rate, config.epochs,
                config.batch_size, config.patience, models_dir,
                config.hidden_dims, run_id_arg, random_seed
            )
        elif mode == "mixed":
            X_train, y_train = data["train_mixed"]
            model, mae = train_model(
                X_train, y_train, X_test, y_test,
                "train_on_mixed", config.learning_rate, config.epochs,
                config.batch_size, config.patience, models_dir,
                config.hidden_dims, run_id_arg, random_seed
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
                config.hidden_dims, run_id_arg, random_seed, syn_model_path
            )
        else:
            raise ValueError(f"Unknown mode: {mode}")
        
        all_maes.append(mae)
    
    # Compute statistics if multiple runs
    if config.n_runs > 1:
        mae_array = np.array(all_maes)
        mae_mean = np.mean(mae_array, axis=0)
        mae_std = np.std(mae_array, axis=0)
        
        print(f"\n{'='*60}")
        print(f"STATISTICS FOR {mode.upper()} MODE ({config.n_runs} runs)")
        print(f"{'='*60}")
        print(f"MAE (mean ± std):")
        print(f"  X: {mae_mean[0]:.2f} ± {mae_std[0]:.2f}")
        print(f"  Y: {mae_mean[1]:.2f} ± {mae_std[1]:.2f}")
        print(f"  Mean: {np.mean(mae_mean):.2f} ± {np.mean(mae_std):.2f}")
        print(f"{'='*60}")
        
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
