from typing import Dict, Tuple, Optional

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split


from .config import (
    REAL_DATA_PATH,
    SYNTHETIC_DATA_DIR,
    NORMALIZATION_FACTOR,
    TEST_SPLIT,
    RANDOM_STATE,
    get_feature_columns,
    TARGET_COLUMNS,
)

# Type alias for data dictionary
DataDict = Dict[str, Tuple[np.ndarray, np.ndarray]]


def normalize_counts(X: np.ndarray) -> np.ndarray:
    """Normalize detector counts by event maximum count."""
    return X / X.max(axis=1, keepdims=True)


def normalize_locations(y: np.ndarray) -> np.ndarray:
    """Normalize locations by normalization factor."""
    return y / NORMALIZATION_FACTOR


def denormalize_locations(y: np.ndarray) -> np.ndarray:
    """Denormalize locations back to original scale."""
    return y * NORMALIZATION_FACTOR


def ensure_synthetic_data(n_samples: int, fov_big: int) -> str:
    """
    Ensure synthetic data exists. Expects data to be generated externally using
    the stix-data-generator repository: https://github.com/i4Ds/stix-data-generator
    
    Args:
        n_samples: Number of synthetic samples
        fov_big: FOV range for synthetic data generation
        
    Returns:
        Path to synthetic data file
        
    Raises:
        FileNotFoundError: If the synthetic data file does not exist
    """
    
    filepath = SYNTHETIC_DATA_DIR / f"sim_{n_samples}_{fov_big}.npz"
    
    if not filepath.exists():
        raise FileNotFoundError(
            f"Synthetic data file not found: {filepath}\n"
            f"Please generate synthetic data using the stix-data-generator repository:\n"
            f"https://github.com/i4Ds/stix-data-generator\n"
            f"Expected file: sim_{n_samples}_{fov_big}.npz in data/synthetic/"
        )
    
    print(f"Using existing synthetic data: {filepath}")
    return str(filepath)


def load_data(syn_data_path: Optional[str], x_fov: float, y_fov: float, sidelobes_threshold: float, col_count: int = None) -> DataDict:
    """
    Load and prepare datasets.
    
    Args:
        syn_data_path: Path to synthetic data npz file (None if not needed)
        x_fov: X FOV bound for filtering (absolute value)
        y_fov: Y FOV bound for filtering (absolute value)
        sidelobes_threshold: Threshold for sidelobes ratio filtering
        col_count: Number of detectors to use (e.g., 12 or 24). If None, uses all 24.
    Returns:
        Dictionary with train/test splits for real, synthetic, and mixed data
    """
    # Load real data
    df_raw = pd.read_csv(REAL_DATA_PATH)

    print(f"Loaded real data: {len(df_raw)} events")
    print(f"Filtering by sidelobes ratio threshold: {sidelobes_threshold}")
    # filter by sidelobes ratio threshold
    df = df_raw[df_raw['sidelobes_ratio'] < sidelobes_threshold]
    print(f"Filtered sidelobes data: {len(df)} events")
    # Filter by FOV
    print(f"Filtering by FOV: x_fov={x_fov}, y_fov={y_fov}")
    df_filtered = df[
        (df["loc_x_stix"] > -x_fov) & (df["loc_x_stix"] < x_fov) &
        (df["loc_y_stix"] > -y_fov) & (df["loc_y_stix"] < y_fov)
    ]
    print(f"Filtered FOV data: {len(df_filtered)} events")
    feature_cols = get_feature_columns(col_count)
  
    
    X_real = df_filtered[feature_cols].values.astype(float)
    y_real = df_filtered[TARGET_COLUMNS].values.astype(float)
    
    print("Normalizing real data...")
    X_real_norm = normalize_counts(X_real)
    y_real_norm = normalize_locations(y_real)
    
    print("Splitting real data into train and test sets...")
    X_train_r, X_test_r, y_train_r, y_test_r = train_test_split(
        X_real_norm, y_real_norm, test_size=TEST_SPLIT, random_state=RANDOM_STATE
    )
    
    result: DataDict = {
        "train_real": (X_train_r, y_train_r),
        "test_real": (X_test_r, y_test_r),
    }
    
    print(f"Real training set X shape: {X_train_r.shape}")
    print(f"Real training set y shape: {y_train_r.shape}")

    print(f"Real test set shape X: {X_test_r.shape}")
    print(f"Real test set shape y: {y_test_r.shape}")
    
    # Load synthetic data if path is provided
    if syn_data_path is not None:
        print(f"Loading synthetic data from: {syn_data_path}")
        data_syn = np.load(syn_data_path)
        print(f"Getting synthetic data with {len(feature_cols)} features")
        X_syn = data_syn['X'][:, :len(feature_cols)]  # Slice to match col_count detectors
        y_syn = data_syn['Y']
        
        print(f"Filtering synthetic data by FOV: x_fov={x_fov}, y_fov={y_fov}")
        # Filter by FOV
        mask = (
            (y_syn[:, 0] > -x_fov) & (y_syn[:, 0] < x_fov) &
            (y_syn[:, 1] > -y_fov) & (y_syn[:, 1] < y_fov)
        )
        X_syn = X_syn[mask]
        y_syn = y_syn[mask]
        
        print("Normalizing synthetic data...")
        X_syn_norm = normalize_counts(X_syn)
        y_syn_norm = normalize_locations(y_syn)
        
        print(f"Synthetic data X shape {X_syn_norm.shape}")
        print(f"Synthetic data Y shape {y_syn_norm.shape}")
        # Create mixed dataset
        print("Mixing real and synthetic training sets...")

        X_train_mixed = np.concatenate((X_train_r, X_syn_norm), axis=0)
        y_train_mixed = np.concatenate((y_train_r, y_syn_norm), axis=0)
        print(f"Mixed training set shape: {X_train_mixed.shape}")
        
        result["train_syn"] = (X_syn_norm, y_syn_norm)
        result["train_mixed"] = (X_train_mixed, y_train_mixed)
    
    return result
