from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from time import perf_counter
from typing import List, Optional, Sequence, Tuple, Union

import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


TorchInput = Union[np.ndarray, Sequence[np.ndarray]]


class TorchMLP(nn.Module):
    def __init__(self, input_dim: int, hidden_dims: List[int]):
        super().__init__()
        layers: List[nn.Module] = []
        prev = input_dim
        for dim in hidden_dims:
            layers.append(nn.Linear(prev, dim))
            layers.append(nn.ReLU())
            prev = dim
        layers.append(nn.Linear(prev, 2))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class ResidualBlock(nn.Module):
    def __init__(self, width: int):
        super().__init__()
        self.fc1 = nn.Linear(width, width)
        self.norm = nn.LayerNorm(width)
        self.act = nn.GELU()
        self.fc2 = nn.Linear(width, width)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.fc2(self.act(self.norm(self.fc1(x))))


class TorchResidualMLP(nn.Module):
    """Residual blocks at constant width; `hidden_dims` = one width per block."""

    def __init__(self, input_dim: int, hidden_dims: List[int]):
        super().__init__()
        if len(set(hidden_dims)) != 1:
            raise ValueError(f"Residual arch needs a constant width, got {hidden_dims}")
        width = hidden_dims[0]
        self.stem = nn.Sequential(nn.Linear(input_dim, width), nn.GELU())
        self.blocks = nn.Sequential(*[ResidualBlock(width) for _ in hidden_dims])
        self.head = nn.Linear(width, 2)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.head(self.blocks(self.stem(x)))


class TorchDualHeadMLP(nn.Module):
    """Shared trunk with a wider Y head."""

    def __init__(self, input_dim: int, hidden_dims: List[int]):
        super().__init__()
        trunk: List[nn.Module] = []
        prev = input_dim
        for dim in hidden_dims:
            trunk.append(nn.Linear(prev, dim))
            trunk.append(nn.ReLU())
            prev = dim
        self.trunk_net = nn.Sequential(*trunk)
        w = prev
        self.x_head = nn.Sequential(nn.Linear(w, w // 4), nn.ReLU(), nn.Linear(w // 4, 1))
        self.y_head = nn.Sequential(
            nn.Linear(w, w // 2),
            nn.ReLU(),
            nn.Linear(w // 2, w // 4),
            nn.ReLU(),
            nn.Linear(w // 4, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        t = self.trunk_net(x)
        return torch.cat([self.x_head(t), self.y_head(t)], dim=1)


ARCH_CLASSES = {
    "mlp": TorchMLP,
    "res": TorchResidualMLP,
    "dualhead_y": TorchDualHeadMLP,
}


def make_loss_fn(loss: str, huber_delta: float):
    if loss == "mse":
        return nn.functional.mse_loss
    if loss == "huber":

        def huber(pred: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
            return nn.functional.huber_loss(pred, target, delta=huber_delta)

        return huber
    raise ValueError(f"Unknown loss: {loss!r}")


def get_device() -> torch.device:
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def _primary_array(X: TorchInput) -> np.ndarray:
    if isinstance(X, (list, tuple)):
        return np.asarray(X[0], dtype=np.float32)
    return np.asarray(X, dtype=np.float32)


def _make_loader(X: TorchInput, y: np.ndarray, batch_size: int, shuffle: bool) -> DataLoader:
    ds = TensorDataset(
        torch.from_numpy(_primary_array(X)),
        torch.from_numpy(np.asarray(y, dtype=np.float32)),
    )
    return DataLoader(ds, batch_size=batch_size, shuffle=shuffle, drop_last=False)


def build_model(
    input_shape: Tuple[int, ...],
    hidden_dims: List[int],
    arch: str = "mlp",
) -> nn.Module:
    try:
        cls = ARCH_CLASSES[arch]
    except KeyError as exc:
        raise ValueError(f"Unknown arch: {arch!r}") from exc
    return cls(input_dim=int(input_shape[0]), hidden_dims=hidden_dims)


def train_model(
    model: nn.Module,
    X_train: TorchInput,
    y_train: np.ndarray,
    X_val: TorchInput,
    y_val: np.ndarray,
    learning_rate: float,
    epochs: int,
    batch_size: int,
    patience: int,
    device: torch.device,
    loss: str = "mse",
    huber_delta: float = 0.05,
) -> nn.Module:
    model = model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    loss_fn = make_loss_fn(loss, huber_delta)
    train_loader = _make_loader(X_train, y_train, batch_size=batch_size, shuffle=True)
    val_loader = _make_loader(X_val, y_val, batch_size=batch_size, shuffle=False)

    best_state = deepcopy(model.state_dict())
    best_val = float("inf")
    wait = 0
    log_every = 25 if epochs >= 100 else max(5, epochs // 10) if epochs > 10 else 1

    print(
        f"[TRAIN][PyTorch] device={device} "
        f"train_samples={len(train_loader.dataset)} val_samples={len(val_loader.dataset)} "
        f"batch_size={batch_size} loss={loss}"
        + (f" (delta={huber_delta})" if loss == "huber" else ""),
        flush=True,
    )

    for epoch in range(1, epochs + 1):
        epoch_start = perf_counter()
        model.train()
        train_loss = 0.0
        train_n = 0
        for xb, yb in train_loader:
            xb = xb.to(device, non_blocking=True)
            yb = yb.to(device, non_blocking=True)
            optimizer.zero_grad(set_to_none=True)
            batch_loss = loss_fn(model(xb), yb)
            batch_loss.backward()
            optimizer.step()
            bs = yb.shape[0]
            train_loss += batch_loss.item() * bs
            train_n += bs
        train_loss /= max(train_n, 1)

        model.eval()
        val_loss = 0.0
        val_n = 0
        with torch.no_grad():
            for xb, yb in val_loader:
                xb = xb.to(device, non_blocking=True)
                yb = yb.to(device, non_blocking=True)
                batch_loss = loss_fn(model(xb), yb)
                bs = yb.shape[0]
                val_loss += batch_loss.item() * bs
                val_n += bs
        val_loss /= max(val_n, 1)

        improved = val_loss < best_val
        if improved:
            best_val = val_loss
            best_state = deepcopy(model.state_dict())
            wait = 0
        else:
            wait += 1

        if epoch == 1 or epoch == epochs or epoch % log_every == 0 or improved or wait >= patience:
            marker = "*" if improved else " "
            print(
                f"[TRAIN][PyTorch]{marker} epoch {epoch:4d}/{epochs} "
                f"train_loss={train_loss:.6f} val_loss={val_loss:.6f} "
                f"best_val={best_val:.6f} wait={wait}/{patience} "
                f"time={perf_counter() - epoch_start:.1f}s",
                flush=True,
            )

        if wait >= patience:
            print(
                f"[TRAIN][PyTorch] Early stopping at epoch {epoch} (best_val={best_val:.6f})",
                flush=True,
            )
            break

    model.load_state_dict(best_state)
    return model


def predict(
    model: nn.Module,
    X: TorchInput,
    device: torch.device,
    batch_size: int = 8192,
) -> np.ndarray:
    model.eval()
    loader = DataLoader(
        TensorDataset(torch.from_numpy(_primary_array(X))),
        batch_size=batch_size,
        shuffle=False,
    )
    preds = []
    with torch.no_grad():
        for (xb,) in loader:
            preds.append(model(xb.to(device, non_blocking=True)).cpu().numpy())
    return np.concatenate(preds, axis=0)


def save_model(
    path: Path,
    model: nn.Module,
    input_shape: Tuple[int, ...],
    hidden_dims: List[int],
    arch: str = "mlp",
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "state_dict": model.state_dict(),
            "input_shape": tuple(int(v) for v in input_shape),
            "hidden_dims": list(hidden_dims),
            "arch": arch,
        },
        path,
    )


def load_model(path: Path, map_location: Optional[torch.device] = None) -> nn.Module:
    payload = torch.load(path, map_location=map_location or "cpu", weights_only=False)
    input_shape = tuple(payload.get("input_shape", (payload.get("input_dim"),)))
    model = build_model(
        input_shape=input_shape,
        hidden_dims=payload["hidden_dims"],
        arch=payload.get("arch", "mlp"),
    )
    model.load_state_dict(payload["state_dict"])
    return model
