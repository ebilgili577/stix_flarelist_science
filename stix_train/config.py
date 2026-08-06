from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import List, Optional

PROJECT_ROOT = Path(__file__).resolve().parent.parent
REAL_DATA_PATH = PROJECT_ROOT / "data/6_final/stix_flarelist_final_20210214_20250830.csv"
SYNTHETIC_DATA_DIR = PROJECT_ROOT / "data/synthetic"
EXPERIMENTS_DIR = PROJECT_ROOT / "experiments"
UNSEEN_DATA_PATH = PROJECT_ROOT / "data/unseen/stix_flarelist_final_20250901_20260121.csv"

NORMALIZATION_FACTOR = 4000.0
TEST_SPLIT = 0.20
RANDOM_STATE = 123

DEFAULT_HIDDEN_DIMS = [1024, 512, 256, 128, 64]

DETECTOR_ORDER: List[str] = [
    "3", "20", "22", "16", "14", "32", "21", "26", "4", "24", "8", "28",
    "15", "27", "31", "6", "30", "2", "25", "5", "23", "7", "29", "1",
]
POSITIONS: List[str] = ["a", "b", "c", "d"]
SIDES: List[str] = ["top", "bot"]


def get_feature_columns(col_count: int = None) -> List[str]:
    detectors = DETECTOR_ORDER[:col_count] if col_count else DETECTOR_ORDER
    return [
        f"{det}_{pos}_{side}"
        for det in detectors
        for pos in POSITIONS
        for side in SIDES
    ]


TARGET_COLUMNS: List[str] = ["loc_x_stix", "loc_y_stix"]

detector_dict = dict(
    zip(DETECTOR_ORDER, [f"{i}{j}" for i in range(10, 2, -1) for j in ["a", "b", "c", "d"]])
)


@dataclass
class TrainConfig:
    experiment: str
    mode: str = "mixed"
    col_count: int = 12
    sidelobes_threshold: float = 0.84
    n_samples: int = 1_000_000
    x_min: float = -2295.0
    x_max: float = 2295.0
    y_min: float = -1878.0
    y_max: float = 1878.0
    learning_rate: float = 1e-3
    epochs: int = 1000
    batch_size: int = 512
    patience: int = 50
    n_runs: int = 1
    hidden_dims: List[int] = field(default_factory=lambda: DEFAULT_HIDDEN_DIMS.copy())
    # PyTorch-only options
    backend: str = "tensorflow"
    arch: str = "mlp"
    loss: str = "mse"
    huber_delta: float = 0.05
    # Mixed: validate on real train rows only (synthetic stays in train)
    val_on_real: bool = True
    synthetic_data_path: Optional[str] = None

    def to_dict(self) -> dict:
        return asdict(self)
