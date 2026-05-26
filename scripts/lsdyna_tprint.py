"""Small helpers for LS-DYNA thermal ``tprint`` files."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

RADIUS_MM = 0.500
THICKNESS_MM = 0.200
V1_NR = 100
V1_NZ = 80

_TIME_RE = re.compile(r"\btime\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)\s+time\s+step")


def parse_tprint(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Parse LS-DYNA ``tprint`` into ``(times_ms, temperatures)``.

    ``temperatures`` has shape ``(n_frames, n_nodes)`` and is indexed by
    ``node_id - 1``. Only the node temperature column is read; flux columns are
    intentionally ignored.
    """
    times: list[float] = []
    frames: list[list[float]] = []
    current: list[float] | None = None
    reading_nodes = False

    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            m = _TIME_RE.search(line)
            if m:
                if current is not None:
                    frames.append(current)
                times.append(float(m.group(1)))
                current = []
                reading_nodes = False
                continue

            if current is None:
                continue

            if line.lstrip().startswith("*TEMPERATURE_NODE"):
                reading_nodes = True
                continue

            if not reading_nodes:
                continue

            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                node_id = int(parts[0])
                temp = float(parts[1])
            except ValueError:
                continue

            if node_id <= 0:
                continue
            if len(current) < node_id:
                current.extend([np.nan] * (node_id - len(current)))
            current[node_id - 1] = temp

    if current is not None:
        frames.append(current)

    if not times or not frames:
        raise ValueError(f"no temperature frames found in {path}")
    if len(times) != len(frames):
        raise ValueError(f"time/frame count mismatch in {path}: {len(times)} vs {len(frames)}")

    n_nodes = max(len(frame) for frame in frames)
    arr = np.full((len(frames), n_nodes), np.nan, dtype=float)
    for k, frame in enumerate(frames):
        arr[k, : len(frame)] = frame
    return np.asarray(times, dtype=float), arr


def reshape_to_grid(T_vec: np.ndarray, nr: int = V1_NR, nz: int = V1_NZ) -> np.ndarray:
    """Reshape a V1-style row-major node vector into ``(nz + 1, nr + 1)``."""
    expected = (nr + 1) * (nz + 1)
    if len(T_vec) < expected:
        raise ValueError(f"temperature vector has {len(T_vec)} nodes, expected {expected}")
    return np.asarray(T_vec[:expected], dtype=float).reshape((nz + 1, nr + 1))
