#!/usr/bin/env python3
"""Extend a frozen field from its own published evolved state, without refitting."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
PREVIOUS = ROOT.parent / "late_growth_continuation"
sys.path.insert(0, str(PREVIOUS))
from run import INITIAL, INITIAL_SHA, RESEARCH, analysis, driver, module, save, snapshot

comparison = module("previous_comparison", PREVIOUS / "assess.py")
restorer = module("published_restorer", ROOT.parent / "frozen_continuations_resolution/restore_checkpoint.py")


def stitch(base_rows, segments, endpoint, interval=0.01):
    """Keep original baselines and reject changed, missing or duplicated clock rows."""
    rows = list(base_rows)
    if not rows:
        raise ValueError("Missing earlier trajectory")
    for segment in segments:
        if len(segment) < 2 or segment[0] != rows[-1]:
            raise ValueError("Restart changed or omitted its boundary evidence")
        rows.extend(segment[1:])
    if len(rows) != round(endpoint / interval) + 1:
        raise ValueError("Incomplete physical observation clock")
    for index, row in enumerate(rows):
        if abs(row["time"] - index * interval) > 1e-12:
            raise ValueError("Physical observation clock changed")
        if index and row["step"] <= rows[index - 1]["step"]:
            raise ValueError("Accepted-step clock did not advance")
    return rows


def provenance(grid, solver=None):
    protocol = json.loads((ROOT / "protocol.json").read_text())
    source = RESEARCH / protocol["source_histories"][str(grid)]
    manifest, rows = comparison.history(source)
    if driver.sha(INITIAL) != INITIAL_SHA or manifest["initial_state_sha256"] != INITIAL_SHA:
        raise ValueError("Frozen initial field changed")
    if rows[-1]["time"] != protocol["start_time"] or manifest["settings"]["grid"] != grid:
        raise ValueError("Wrong source trajectory")
    settings = manifest["settings"]
    for key, expected in {"maximum_dt": protocol["maximum_time_step"],
                          "viscosity": protocol["viscosity"], "cfl": protocol["cfl_target"],
                          "observation_interval": protocol["observation_interval"],
                          "dense_sampling_grid": protocol["common_physical_sampling_grid"],
                          "cutoff_limit": protocol["retained_gates"]["maximum_cutoff_energy_fraction_each_step"]}.items():
        if settings[key] != expected:
            raise ValueError("Different source setting: " + key)
    for name, digest in manifest["source_sha256"].items():
        if driver.sha(RESEARCH / name) != digest:
            raise ValueError("Validated source changed: " + name)
    if solver is not None and driver.sha(solver) != manifest["solver_sha256"]:
        raise ValueError("Executable differs from the verified preceding run")
    return protocol, source, manifest, rows


def pack_checkpoint(path, output, budget):
    with output.open("xb") as stream:
        stream.write(gzip.compress(path.read_bytes(), compresslevel=6, mtime=0))
    return {"time": budget["time"], "path": output.name, "sha256": driver.sha(output),
            "uncompressed_sha256": budget["checkpoint_sha256"],
            "uncompressed_bytes": budget["checkpoint_bytes"]}


def restore_source(protocol, original, grid, output):
    archive = protocol["source_checkpoints"][str(grid)]
    source = RESEARCH / archive["path"]
    if archive["format"] == "split-gzip":
        restored = restorer.restore(source, output)
        if restored["uncompressed_sha256"] != original["uncompressed_sha256"] or restored["grid"] != grid:
            raise ValueError("Published parts point at a different restart state")
    elif archive["format"] == "gzip":
        if driver.sha(source) != original["sha256"]:
            raise ValueError("Published compressed state changed")
        raw = gzip.decompress(source.read_bytes())
        if len(raw) != original["uncompressed_bytes"] or hashlib.sha256(raw).hexdigest() != original["uncompressed_sha256"]:
            raise ValueError("Restored state identity changed")
        with output.open("xb") as destination:
            destination.write(raw)
    else:
        raise ValueError("Unknown published checkpoint encoding")
    return {**archive, "sha256": driver.sha(source)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", required=True)
    parser.add_argument("--grid", type=int, choices=(64, 128), required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    solver = Path(args.solver).resolve(strict=True)
    protocol, source, previous, earlier_rows = provenance(args.grid, solver)
    analysis.test_exact_shear()
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)
    checkpoint = output / "current.chk"
    original = next(item for item in previous["checkpoints"] if item["time"] == protocol["start_time"])
    archive = restore_source(protocol, original, args.grid, checkpoint)
    restart_budget = snapshot(checkpoint, earlier_rows[-1])
    if restart_budget["checkpoint_sha256"] != original["uncompressed_sha256"]:
        raise ValueError("Restored field differs from the verified .10 state")
    manifest = {
        "status": "running", "initial_state_sha256": INITIAL_SHA,
        "solver_sha256": driver.sha(solver), "protocol_sha256": driver.sha(ROOT / "protocol.json"),
        "extension_source_sha256": driver.sha(Path(__file__)),
        "restorer_source_sha256": driver.sha(Path(restorer.__file__)),
        "comparison_source_sha256": driver.sha(Path(comparison.__file__)),
        "source_sha256": previous["source_sha256"],
        "source_history": {"path": str(source.relative_to(RESEARCH)),
                           "manifest_sha256": driver.sha(source / "manifest.json"),
                           "checkpoint_archive": archive,
                           "checkpoint_sha256": original["uncompressed_sha256"]},
        "settings": {**previous["settings"], "start_time": protocol["start_time"],
                     "endpoints": protocol["endpoints"]},
        "numpy_version": analysis.np.__version__,
        "analytic_shear_check": "passed", "restart_verification": restart_budget,
        "runs": [], "budgets": [], "checkpoints": [],
    }
    save(output / "manifest.json", manifest)
    segments = []
    try:
        for endpoint in protocol["endpoints"]:
            stem = "t%03d" % round(endpoint * 1000)
            evidence, log = output / (stem + ".csv"), output / (stem + ".log")
            command = [str(solver), "--restart", str(checkpoint), "--final-time", str(endpoint),
                       "--output", str(evidence), "--checkpoint-output", str(checkpoint)]
            run = {"endpoint": endpoint, "command": command, "restart_sha256": driver.sha(checkpoint)}
            manifest["runs"].append(run)
            save(output / "manifest.json", manifest)
            print("Starting grid", args.grid, "to", endpoint, flush=True)
            started = time.monotonic()
            with log.open("w") as stream:
                completed = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
            run.update(returncode=completed.returncode, wall_seconds=time.monotonic() - started,
                       log_sha256=driver.sha(log))
            if evidence.exists():
                run["evidence_sha256"] = driver.sha(evidence)
            if completed.returncode:
                manifest["status"] = "cutoff-stopped" if completed.returncode == 2 else "solver-failed"
                save(output / "manifest.json", manifest)
                raise RuntimeError("Continuation stopped; last checkpoint and failure evidence retained")
            rows = driver.read_evidence(evidence)
            segments.append(rows)
            stitch(earlier_rows, segments, endpoint, protocol["observation_interval"])
            budget = snapshot(checkpoint, rows[-1])
            run["checkpoint_sha256"] = budget["checkpoint_sha256"]
            manifest["budgets"].append(budget)
            manifest["checkpoints"].append(pack_checkpoint(checkpoint, output / (stem + ".chk.gz"), budget))
            save(output / "manifest.json", manifest)
            print(args.grid, endpoint, "H_ratio", rows[-1]["h_half_ratio"],
                  "L3_ratio", rows[-1]["l3_dense_ratio"], "vorticity_ratio", rows[-1]["vorticity_dense_ratio"],
                  "gamma_H", budget["h_half_logarithmic_rate"], "cutoff", rows[-1]["peak_cutoff_fraction"], flush=True)
        if driver.sha(INITIAL) != INITIAL_SHA:
            raise ValueError("Initial coefficients changed during the extension")
        manifest["status"] = "completed-trajectory-awaiting-crosschecks"
        save(output / "manifest.json", manifest)
        checkpoint.unlink()
    except Exception as error:
        if manifest["status"] == "running":
            manifest["status"] = "verification-failed"
        manifest["error"] = str(error)
        save(output / "manifest.json", manifest)
        raise


if __name__ == "__main__":
    main()
