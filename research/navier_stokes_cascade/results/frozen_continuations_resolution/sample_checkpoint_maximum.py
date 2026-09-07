#!/usr/bin/env python3
"""Check a saved field's vorticity maximum on nested grids; no evolution."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

source = Path(__file__).resolve().parents[1] / "frozen_continuations/analyze_critical_budget.py"
spec = importlib.util.spec_from_file_location("budget", source)
budget = importlib.util.module_from_spec(spec)
spec.loader.exec_module(budget)


def maximum(coefficients, sampling_grid):
    n = coefficients.shape[0]
    m = sampling_grid
    if m < n or m % n or m & (m - 1):
        raise ValueError("Use a nested power-of-two sampling grid")
    frequencies = np.rint(np.fft.fftfreq(n) * n).astype(int)
    k = (frequencies[:, None, None], frequencies[None, :, None], np.arange(n // 2 + 1)[None, None, :])
    positive = coefficients[:, :, :n // 2 + 1, :]
    locations = np.ix_(frequencies % m, frequencies % m, np.arange(n // 2 + 1))
    squared = np.zeros((m, m, m))
    for component in range(3):
        a, b = (component + 1) % 3, (component + 2) % 3
        padded = np.zeros((m, m, m // 2 + 1), dtype=complex)
        padded[locations] = 1j * (k[a] * positive[..., b] - k[b] * positive[..., a])
        values = np.fft.irfftn(padded, s=(m, m, m), axes=(0, 1, 2)) * m**3
        del padded
        squared += values**2
        del values
    return float(np.sqrt(squared.max()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--initial-state", required=True)
    parser.add_argument("--sampling-grid", type=int, default=256)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    output = Path(args.output)
    if output.exists():
        raise ValueError("Refusing to replace existing sampling evidence")
    field, viscosity, observed, metadata = budget.load_checkpoint(Path(args.checkpoint), args.sha256)
    n, m = metadata["grid"], args.sampling_grid
    native = maximum(field, n)
    if abs(native - observed[6]) / observed[6] > 1e-10:
        raise ValueError("Independent real FFT samples disagree with FFTW at the native grid")
    dense = maximum(field, m)
    initial_path = Path(args.initial_state)
    initial = budget.load_initial(initial_path, n)
    initial_native, initial_dense = maximum(initial, n), maximum(initial, m)
    if dense < native * (1 - 1e-12) or initial_dense < initial_native * (1 - 1e-12):
        raise ValueError("Nested samples unexpectedly reduced the maximum")
    record = {
        "time": observed[0], "evolution_grid": n, "cutoff": metadata["cutoff"], "viscosity": viscosity,
        "sampling_grids": [n, m], "final_maxima": [native, dense],
        "initial_maxima": [initial_native, initial_dense],
        "amplification_ratios": [native / initial_native, dense / initial_dense],
        "final_relative_sampling_shift": abs(dense - native) / dense,
        "native_numpy_vs_fftw_relative_error": abs(native - observed[6]) / observed[6],
        "checkpoint_sha256": args.sha256,
        "initial_state_sha256": hashlib.sha256(initial_path.read_bytes()).hexdigest(),
        "analysis_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "reader_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "numpy_version": np.__version__, "command": sys.argv,
        "scope": "Endpoint samples of the same truncated field; no additional evolution modes or uniform continuum bound",
    }
    output.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    print(json.dumps(record, indent=2), flush=True)


if __name__ == "__main__":
    main()
