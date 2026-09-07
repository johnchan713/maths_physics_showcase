#!/usr/bin/env python3
"""Continue frozen finalists on a fixed observation clock, preserving final states.

These are empirical screening gates for finite amplification, not error bounds
for the infinite-dimensional PDE. A resolution failure pauses this grid pair;
it does not reject the underlying initial field as a singularity mechanism.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def read_evidence(path):
    with Path(path).open() as source:
        raw = list(csv.DictReader(source))
    if not raw:
        raise ValueError("Empty continuation evidence")
    rows = []
    for item in raw:
        if None in item or None in item.values():
            raise ValueError("Malformed continuation evidence row")
        row = {key: value if key in ("backend", "cutoff_stopped") else float(value)
               for key, value in item.items()}
        if row["backend"] not in ("fftw", "radix2") or row["cutoff_stopped"] not in ("true", "false"):
            raise ValueError("Invalid backend/stop identifier")
        if any(not math.isfinite(v) for v in row.values() if isinstance(v, float)):
            raise ValueError("Nonfinite continuation evidence")
        if rows and (row["time"] <= rows[-1]["time"] or row["step"] <= rows[-1]["step"]):
            raise ValueError("Continuation clock did not advance")
        rows.append(row)
    return rows


def assess_pair(coarse, fine, endpoint, interval=0.01):
    failures = []
    expected = round(endpoint / interval) + 1
    if len(coarse) != expected or len(fine) != expected:
        return {"status": "resolution-or-run-failure", "failures": ["incomplete observation clock"]}
    for rows in (coarse, fine):
        for index, row in enumerate(rows):
            if abs(row["time"] - index * interval) > 1e-12:
                failures.append("inconsistent observation clock")
            if row["cutoff_stopped"] != "false" or row["peak_cutoff_fraction"] > 0.008:
                failures.append("cutoff gate")
            if row["divergence_defect"] > 1e-9 or row["reality_defect"] > 1e-9:
                failures.append("Fourier invariant gate")
    gaps = {key: max(relative(a[key], b[key]) for a, b in zip(coarse, fine)) for key in
            ("h_half_ratio", "l3_dense_ratio", "enstrophy_ratio", "k_rms_ratio", "vorticity_dense_ratio")}
    for key, tolerance in (("h_half_ratio", 0.02), ("l3_dense_ratio", 0.02),
                           ("enstrophy_ratio", 0.05), ("k_rms_ratio", 0.05),
                           ("vorticity_dense_ratio", 0.10)):
        if gaps[key] > tolerance:
            failures.append("cross-resolution " + key)
    late_start = max(interval, endpoint - 0.02)
    paired_late = [(a, b) for a, b in zip(coarse, fine) if a["time"] >= late_start - 1e-12]
    gaps["late_stretching"] = max(relative(a["stretching"], b["stretching"]) for a, b in paired_late)
    if gaps["late_stretching"] > 0.10:
        failures.append("cross-resolution late stretching")
    numerical_failures = sorted(set(failures))
    late_gains = []
    for rows in (coarse, fine):
        late = [row for row in rows if row["time"] >= late_start - 1e-12]
        late_gains.append(late[-1]["h_half"] / late[0]["h_half"])
        if late_gains[-1] < 1.001 or any(b["h_half"] < a["h_half"] * (1 - 1e-10) for a, b in zip(late, late[1:])):
            failures.append("late critical-norm growth stalled or reversed")
        if min(row["production_to_dissipation"] for row in late) <= 1.0:
            failures.append("late stretching no longer exceeds viscous destruction")
    return {
        "status": "resolution-or-run-failure" if numerical_failures else
                  "finite-growth-stalled" if failures else "continue-finite-amplification",
        "failures": sorted(set(failures)), "numerical_failures": numerical_failures,
        "maximum_trajectory_relative_gaps": gaps,
        "late_h_half_gains": late_gains,
        "sampling_shift": {key: max(relative(row[key], row[dense]) for rows in (coarse, fine) for row in rows)
                           for key, dense in (("l3_sample", "l3_dense"), ("vorticity_sample", "vorticity_dense"))},
        "endpoints": [{key: rows[-1][key] for key in
                       ("grid", "time", "step", "h_half_ratio", "l3_dense_ratio", "vorticity_dense_ratio",
                        "enstrophy_ratio", "k_rms_ratio", "production_to_dissipation", "peak_cutoff_fraction")}
                      for rows in (coarse, fine)],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--candidates", nargs="+", default=["checkpoint", "wave-packets-0", "vortex-tubes-0"],
                        choices=["checkpoint", "wave-packets-0", "vortex-tubes-0"])
    args = parser.parse_args()
    research = Path(__file__).resolve().parents[1]
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise SystemExit("Refusing to replace an existing continuation campaign")
    solver = str(Path(args.solver).resolve())
    manifest = {
        "purpose": "Frozen-field finite-amplification continuation; no PDE singularity claim",
        "solver_sha256": sha(solver),
        "source_sha256": {str(p.relative_to(research)): sha(p) for p in
                          [research / "src/continue.cpp", Path(__file__).resolve(),
                           *sorted((research / "include/ns_cascade").glob("*.hpp"))]},
        "settings": {"grids": [32, 64], "viscosity": 0.02, "maximum_dt": 0.000125,
                     "cfl": 0.4, "observation_interval": 0.01, "cutoff_limit": 0.008,
                     "sampling_grid": 64, "dense_sampling_grid": 128,
                     "stage_endpoints": [0.08, 0.10, 0.12, 0.16]},
        "gates": {"h_half_relative_gap": 0.02, "l3_relative_gap": 0.02,
                  "enstrophy_relative_gap": 0.05, "k_rms_relative_gap": 0.05,
                  "vorticity_relative_gap": 0.10, "late_stretching_relative_gap": 0.10,
                  "late_h_half_gain_minimum": 1.001, "late_production_to_dissipation_minimum": 1.0,
                  "late_window_duration": 0.02,
                  "sampling_shift_is_reported_separately": True},
        "candidates": [],
    }

    def persist():
        temporary = output / "manifest.json.tmp"
        temporary.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
        temporary.replace(output / "manifest.json")

    persist()
    for name in args.candidates:
        source = research / "results/robust_amplification_verified" / (name + "-discovery") / "state.csv"
        candidate = {"name": name, "initial_state": str(source.relative_to(research)),
                     "initial_state_sha256": sha(source), "runs": [], "stages": []}
        manifest["candidates"].append(candidate)
        histories = {32: [], 64: []}
        for endpoint in manifest["settings"]["stage_endpoints"]:
            calls = []
            for grid in (32, 64):
                directory = output / name / ("grid-" + str(grid))
                directory.mkdir(parents=True, exist_ok=True)
                checkpoint = directory / "current.chk"
                stem = "t" + str(round(endpoint * 1000)).zfill(3)
                evidence, log = directory / (stem + ".csv"), directory / (stem + ".log")
                command = [solver, "--final-time", str(endpoint), "--output", str(evidence),
                           "--checkpoint-output", str(checkpoint)]
                if histories[grid]:
                    command += ["--restart", str(checkpoint)]
                else:
                    command += ["--state-input", str(source), "--grid", str(grid),
                                "--sampling-grid", "64", "--dense-sampling-grid", "128",
                                "--viscosity", "0.02", "--dt", "0.000125",
                                "--observation-interval", "0.01", "--cutoff-limit", "0.008"]
                run = {"grid": grid, "endpoint": endpoint, "command": command,
                       "evidence": str(evidence.relative_to(output)), "log": str(log.relative_to(output))}
                candidate["runs"].append(run)
                calls.append((run, evidence, log))
            persist()

            def execute(call):
                run, evidence, log = call
                started = time.monotonic()
                with log.open("w") as stream:
                    result = subprocess.run(run["command"], stdout=stream, stderr=subprocess.STDOUT, check=False)
                return result.returncode, time.monotonic() - started, read_evidence(evidence) if evidence.exists() else []

            print(f"Continuing {name} on 32/64 to {endpoint}", flush=True)
            with ThreadPoolExecutor(max_workers=2) as pool:
                results = list(pool.map(execute, calls))
            for (run, evidence, log), (returncode, seconds, rows) in zip(calls, results):
                run.update(returncode=returncode, wall_seconds=seconds, log_sha256=sha(log))
                if evidence.exists():
                    run["evidence_sha256"] = sha(evidence)
                grid = run["grid"]
                if histories[grid] and rows:
                    if histories[grid][-1] != rows[0]:
                        raise RuntimeError("Restart changed its boundary evidence")
                    rows = rows[1:]
                histories[grid].extend(rows)
            assessment = assess_pair(histories[32], histories[64], endpoint)
            if any(r[0] != 0 for r in results):
                assessment["status"] = "resolution-or-run-failure"
                assessment["failures"].append("solver exit gate")
            candidate["stages"].append({"endpoint": endpoint, **assessment})
            persist()
            print(name, endpoint, assessment["status"], assessment["failures"], flush=True)
            if assessment["status"] != "continue-finite-amplification":
                break
        candidate["checkpoints"] = []
        for grid in (32, 64):
            checkpoint = output / name / ("grid-" + str(grid)) / "current.chk"
            if checkpoint.exists():
                data = checkpoint.read_bytes()
                packed = checkpoint.with_suffix(".chk.gz")
                packed.write_bytes(gzip.compress(data, compresslevel=9, mtime=0))
                candidate["checkpoints"].append({"grid": grid, "path": str(packed.relative_to(output)),
                                                  "sha256": sha(packed), "uncompressed_sha256": sha(checkpoint),
                                                  "uncompressed_bytes": len(data)})
                checkpoint.unlink()
        persist()
    print("Continuation campaign recorded with checksummed final states", flush=True)


if __name__ == "__main__":
    main()
