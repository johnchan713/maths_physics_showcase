#!/usr/bin/env python3
"""Reproduce a bounded, frozen-field continuation; never optimize the holdouts.

Run one grid at a time. Cross-grid promotion is assessed separately, so a
completed trajectory is not itself a passed resolution or blow-up test.
"""
import argparse
import gzip
import importlib.util
import json
from pathlib import Path
import struct
import subprocess
import time

import numpy as np

RESEARCH = Path(__file__).resolve().parents[2]
INITIAL = RESEARCH / "candidates/wave_k3_late_rate_t004_screening.csv"
INITIAL_SHA = "893ef61450d670936dc851d28ab3cf7c1aa315c6ff23934646e74137c40ac6ef"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


driver = module("continuation_driver", RESEARCH / "scripts/continue_candidates.py")
analysis = module("critical_budget", RESEARCH / "results/frozen_continuations/analyze_critical_budget.py")


def save(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def snapshot(path, row):
    """Check the emitted checkpoint against its CSV, then use NumPy's own FFTs."""
    raw = path.read_bytes()
    if raw[:8] != b"NSCONT1\n":
        raise ValueError("Unknown checkpoint format")
    n, cutoff, sampling, dense, backend = struct.unpack_from("<5Q", raw, 24)
    length = struct.unpack_from("<Q", raw, 8)[0]
    state_size = struct.unpack_from("<Q", raw, 488)[0]
    if length != len(raw) - 24 or state_size != n**3 or len(raw) != 496 + 48 * state_size:
        raise ValueError("Checkpoint length mismatch")
    if 3 * cutoff >= n or backend != 1 or n != row["grid"]:
        raise ValueError("Expected a dealiased FFTW checkpoint on the recorded grid")
    viscosity = struct.unpack_from("<d", raw, 64)[0]
    observed = struct.unpack_from("<19d", raw, 336)
    if observed[0] != row["time"] or observed[4] != row["h_half"]:
        raise ValueError("Checkpoint and CSV endpoint differ")
    field = np.frombuffer(raw, dtype="<c16", offset=496).reshape(n, n, n, 3)
    if not np.isfinite(field).all():
        raise ValueError("Nonfinite field")
    result = analysis.budget(field, viscosity)
    names = {"energy": "energy", "enstrophy": "enstrophy", "h_half": "h_half",
             "enstrophy_stretching": "stretching",
             "enstrophy_viscous_destruction": "viscous_destruction"}
    if sampling == n:
        names["sampled_vorticity_native"] = "vorticity_sample"
    errors = {key: driver.relative(result[key], row[value]) for key, value in names.items()}
    if max(errors.values()) > 1e-10:
        raise ValueError("Independent snapshot diagnostics disagree: " + str(errors))
    return {"time": row["time"], "grid": n, "cutoff": cutoff,
            "checkpoint_sha256": driver.sha(path), "checkpoint_bytes": len(raw),
            "numpy_vs_csv_relative_errors": errors, **result}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--grid", type=int, choices=[32, 64, 128], required=True)
    parser.add_argument("--max-dt", type=float, default=0.000125)
    parser.add_argument("--observation-interval", type=float, default=0.01)
    args = parser.parse_args()
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)
    solver = str(Path(args.solver).resolve(strict=True))
    if driver.sha(INITIAL) != INITIAL_SHA:
        raise ValueError("Frozen initial field has changed")
    analysis.test_exact_shear()
    sampling = max(64, args.grid)
    manifest = {
        "status": "running", "initial_state": str(INITIAL.relative_to(RESEARCH)),
        "initial_state_sha256": INITIAL_SHA, "solver_sha256": driver.sha(solver),
        "source_sha256": {str(p.relative_to(RESEARCH)): driver.sha(p) for p in
                          [Path(__file__).resolve(), RESEARCH / "src/continue.cpp",
                           RESEARCH / "scripts/continue_candidates.py",
                           RESEARCH / "results/frozen_continuations/analyze_critical_budget.py",
                           *sorted((RESEARCH / "include/ns_cascade").glob("*.hpp"))]},
        "settings": {"grid": args.grid, "sampling_grid": sampling, "dense_sampling_grid": 128,
                     "maximum_dt": args.max_dt, "observation_interval": args.observation_interval,
                     "viscosity": 0.02, "cutoff_limit": 0.008, "cfl": 0.4,
                     "endpoints": [0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10]},
        "numpy_version": np.__version__, "analytic_shear_check": "passed",
        "runs": [], "budgets": [], "checkpoints": [],
    }
    save(output / "manifest.json", manifest)
    histories = []
    checkpoint = output / "current.chk"
    for endpoint in manifest["settings"]["endpoints"]:
        stem = "t" + str(round(endpoint * 1000)).zfill(3)
        evidence, log = output / (stem + ".csv"), output / (stem + ".log")
        command = [solver, "--final-time", str(endpoint), "--output", str(evidence),
                   "--checkpoint-output", str(checkpoint)]
        restart_sha = None
        if histories:
            restart_sha = driver.sha(checkpoint)
            command += ["--restart", str(checkpoint)]
        else:
            command += ["--state-input", str(INITIAL), "--grid", str(args.grid),
                        "--sampling-grid", str(sampling), "--dense-sampling-grid", "128",
                        "--viscosity", "0.02", "--dt", str(args.max_dt),
                        "--observation-interval", str(args.observation_interval),
                        "--cutoff-limit", "0.008", "--backend", "fftw"]
        record = {"endpoint": endpoint, "command": command, "restart_sha256": restart_sha}
        manifest["runs"].append(record)
        save(output / "manifest.json", manifest)
        started = time.monotonic()
        with log.open("w") as stream:
            completed = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
        record.update(returncode=completed.returncode, wall_seconds=time.monotonic() - started,
                      log_sha256=driver.sha(log))
        save(output / "manifest.json", manifest)
        if completed.returncode != 0:
            manifest["status"] = "cutoff-stopped" if completed.returncode == 2 else "solver-failed"
            save(output / "manifest.json", manifest)
            raise RuntimeError("Continuation stopped; failure and last checkpoint retained")
        rows = driver.read_evidence(evidence)
        record["evidence_sha256"] = driver.sha(evidence)
        if histories and histories[-1] != rows[0]:
            raise ValueError("Restart changed boundary evidence")
        histories.extend(rows[1:] if histories else rows)
        expected = round(endpoint / args.observation_interval) + 1
        if len(histories) != expected or abs(histories[-1]["time"] - endpoint) > 1e-12:
            raise ValueError("Incomplete physical observation clock")
        budget = snapshot(checkpoint, rows[-1])
        manifest["budgets"].append(budget)
        if endpoint in (0.08, 0.10):
            packed = output / (stem + ".chk.gz")
            packed.write_bytes(gzip.compress(checkpoint.read_bytes(), compresslevel=6, mtime=0))
            manifest["checkpoints"].append({"time": endpoint, "path": packed.name,
                "sha256": driver.sha(packed), "uncompressed_sha256": driver.sha(checkpoint),
                "uncompressed_bytes": checkpoint.stat().st_size})
        record["checkpoint_sha256"] = budget["checkpoint_sha256"]
        save(output / "manifest.json", manifest)
        print(args.grid, endpoint, "H_ratio", rows[-1]["h_half_ratio"],
              "gamma_H", budget["h_half_logarithmic_rate"],
              "cutoff", rows[-1]["peak_cutoff_fraction"], flush=True)
    if driver.sha(INITIAL) != INITIAL_SHA:
        raise ValueError("Initial coefficients changed during experiment")
    manifest["status"] = "completed-trajectory-awaiting-crosschecks"
    save(output / "manifest.json", manifest)
    checkpoint.unlink()


if __name__ == "__main__":
    main()
