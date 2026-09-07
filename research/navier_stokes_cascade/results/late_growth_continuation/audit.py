#!/usr/bin/env python3
"""Verify archived run identity, completeness, and the stated numerical verdict."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import struct

from assess import history, paired, timestep
from run import INITIAL, INITIAL_SHA, RESEARCH, driver, save

ROOT = Path(__file__).resolve().parent


def main():
    if driver.sha(INITIAL) != INITIAL_SHA:
        raise ValueError("Candidate identity changed")
    report = {"initial_state_sha256": INITIAL_SHA, "checkpoints": [], "baseline_reproduction": {},
              "scope": "Archive integrity and floating-point screening gates; no rigorous PDE certificate"}
    data = {}
    for name in ("reference32", "reference64", "refined64"):
        folder = ROOT / name
        manifest, rows = data[name] = history(folder)
        if manifest["initial_state_sha256"] != INITIAL_SHA:
            raise ValueError("A run used a different candidate")
        for path, expected in manifest["source_sha256"].items():
            if driver.sha(RESEARCH / path) != expected:
                raise ValueError("Recorded source changed: " + path)
        times = [row["time"] for row in manifest["budgets"]]
        if times != manifest["settings"]["endpoints"]:
            raise ValueError("Critical budgets incomplete")
        for run, budget in zip(manifest["runs"], manifest["budgets"]):
            log = folder / ("t%03d.log" % round(run["endpoint"] * 1000))
            if driver.sha(log) != run["log_sha256"] or run["checkpoint_sha256"] != budget["checkpoint_sha256"]:
                raise ValueError("Run provenance changed")
            if max(budget["numpy_vs_csv_relative_errors"].values()) > 1e-10:
                raise ValueError("Snapshot crosscheck failed")
        if {p["time"] for p in manifest["checkpoints"]} != {0.08, 0.10}:
            raise ValueError("Missing evolved holdout state")
        for checkpoint in manifest["checkpoints"]:
            packed = folder / checkpoint["path"]
            if driver.sha(packed) != checkpoint["sha256"]:
                raise ValueError("Compressed checkpoint changed")
            raw = gzip.decompress(packed.read_bytes())
            if len(raw) != checkpoint["uncompressed_bytes"] or hashlib.sha256(raw).hexdigest() != checkpoint["uncompressed_sha256"]:
                raise ValueError("Restored checkpoint changed")
            grid = struct.unpack_from("<Q", raw, 24)[0]
            time = struct.unpack_from("<d", raw, 336)[0]
            if raw[:8] != b"NSCONT1\n" or grid != manifest["settings"]["grid"] or time != checkpoint["time"]:
                raise ValueError("Checkpoint metadata differ")
            if len(raw) != 496 + 48 * grid**3:
                raise ValueError("Checkpoint schema length differs")
            report["checkpoints"].append({"path": str(packed.relative_to(ROOT)), "grid": grid, "time": time,
                "restored_sha256": checkpoint["uncompressed_sha256"], "bytes": len(raw)})
        if name.startswith("reference"):
            grid = manifest["settings"]["grid"]
            baseline = 1.089587921450797 if grid == 32 else 1.0895966877726906
            measured = next(row["h_half_ratio"] for row in rows if abs(row["time"] - 0.06) < 1e-12)
            gap = driver.relative(baseline, measured)
            if gap > 1e-10:
                raise ValueError("Earlier .06 result was not reproduced")
            report["baseline_reproduction"][str(grid)] = {"relative_h_half_ratio_gap": gap}
    assessments = {str(t): paired(data["reference32"], data["reference64"], t) for t in (0.08, 0.10)}
    refinements = {str(t): timestep(data["reference64"], data["refined64"], t) for t in (0.08, 0.10)}
    if any(item["failures"] for item in [*assessments.values(), *refinements.values()]):
        raise ValueError("Numerical screening failed; do not label the campaign passed")
    oracle = ROOT / "independent64"
    record = json.loads((oracle / "manifest.json").read_text())
    if record["status"] != "passed" or record["returncode"] or record["initial_state_sha256"] != INITIAL_SHA:
        raise ValueError("Independent trajectory failed or used another field")
    for path, expected in record["source_sha256"].items():
        if driver.sha(RESEARCH / path) != expected:
            raise ValueError("Independent solver source changed: " + path)
    if driver.sha(oracle / "run.log") != record["log_sha256"] or driver.sha(oracle / "trajectory.csv") != record["evidence_sha256"]:
        raise ValueError("Independent evidence changed")
    log = (oracle / "run.log").read_text()
    if "verdict: pass" not in log:
        raise ValueError("Missing independent verdict")
    with (oracle / "trajectory.csv").open() as source:
        rows = list(csv.DictReader(source))
    if not rows or abs(float(rows[-1]["time"]) - 0.10) > 1e-12:
        raise ValueError("Independent comparison did not reach the holdout")
    report["independent64"] = {
        "steps": int(rows[-1]["step"]),
        "peak_relative_state_difference": float(re.search(r"peak relative state difference: (\S+)", log)[1]),
        "peak_diagnostic_scaled_difference": float(re.search(r"peak diagnostic scaled difference: (\S+)", log)[1]),
    }
    report["status"] = "passed-declared-finite-amplification-gates"
    report["material_limitation"] = "32/64 vorticity and late-stretching gaps approach 10%; no 128 dynamics in this campaign"
    report["maximum_snapshot_diagnostic_relative_gap"] = max(max(b["numpy_vs_csv_relative_errors"].values())
        for manifest, _ in data.values() for b in manifest["budgets"])
    report["files"] = [{"path": str(p.relative_to(ROOT)), "bytes": p.stat().st_size, "sha256": driver.sha(p)}
        for p in sorted(ROOT.rglob("*")) if p.is_file() and p.name != "audit.json" and
        "__pycache__" not in p.parts and not p.name.endswith((".tmp", ".chk"))]
    save(ROOT / "audit.json", report)
    print(json.dumps({key: value for key, value in report.items() if key != "files"}, indent=2))


if __name__ == "__main__":
    main()
