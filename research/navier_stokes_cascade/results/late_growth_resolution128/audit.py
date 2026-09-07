#!/usr/bin/env python3
"""Audit archive provenance and report the actual 64/128 outcome, including failure."""
import gzip
import hashlib
import json
import math
from pathlib import Path
import struct

from analyze import ROOT, PREVIOUS, INITIAL, INITIAL_SHA, RESEARCH, comparison, driver, save


def main():
    if driver.sha(INITIAL) != INITIAL_SHA:
        raise ValueError("Frozen initial coefficients changed")
    coarse, fine = comparison.history(PREVIOUS / "reference64"), comparison.history(ROOT / "reference128")
    recorded = json.loads((ROOT / "assessment.json").read_text())
    expected_pairs = {str(t): comparison.paired(coarse, fine, t) for t in (0.08, 0.10)}
    if recorded["pairs"] != expected_pairs:
        raise ValueError("Recorded spatial assessment differs from archived evidence")
    if recorded["coarse_manifest_sha256"] != driver.sha(PREVIOUS / "reference64/manifest.json"):
        raise ValueError("Coarse provenance changed")
    if recorded["fine_manifest_sha256"] != driver.sha(ROOT / "reference128/manifest.json"):
        raise ValueError("Fine provenance changed")
    if recorded["protocol_sha256"] != driver.sha(ROOT / "protocol.json"):
        raise ValueError("Declared protocol changed")
    if recorded["analysis_sha256"] != driver.sha(ROOT / "analyze.py"):
        raise ValueError("Analysis source changed")
    sources = {
        "comparison_source_sha256": PREVIOUS / "assess.py",
        "sampler_source_sha256": ROOT.parent / "frozen_continuations_resolution/sample_checkpoint_maximum.py",
        "checkpoint_reader_sha256": ROOT.parent / "frozen_continuations/analyze_critical_budget.py",
    }
    for key, path in sources.items():
        if recorded[key] != driver.sha(path):
            raise ValueError("Analysis dependency changed: " + key)
    if recorded["status"] == "analysis-in-progress":
        raise ValueError("Sampling analysis incomplete")
    if len(recorded["sampling"]) != 4 or {(x["grid"], x["time"]) for x in recorded["sampling"]} != {(n, t) for n in (64, 128) for t in (0.08, 0.10)}:
        raise ValueError("Sampling observations incomplete")
    protocol = json.loads((ROOT / "protocol.json").read_text())
    limits = protocol["additional_sampling"]
    sampling_failures = []
    shared_baselines = recorded["sampling"][0]["initial_maxima"]
    for item in recorded["sampling"]:
        values, baselines = item["final_maxima"], item["initial_maxima"]
        if item["sampling_grids"] != limits["sample_grids"] or len(values) != 2 or len(baselines) != 2 or baselines != shared_baselines:
            raise ValueError("Sampling grids or common initial normalization changed")
        if any(not math.isfinite(v) or v <= 0 for v in values + baselines):
            raise ValueError("Invalid sampled maximum")
        if values[1] < values[0] * (1 - 1e-12) or baselines[1] < baselines[0] * (1 - 1e-12):
            raise ValueError("Nested spatial samples lost an existing maximum")
        ratios = [v / b for v, b in zip(values, baselines)]
        if ratios != item["amplification_ratios"] or driver.relative(*values) != item["raw_relative_shift"] or driver.relative(*ratios) != item["ratio_relative_shift"]:
            raise ValueError("Sampling ratios or shifts differ from recorded maxima")
        original, history_rows = coarse if item["grid"] == 64 else fine
        endpoint = next(row for row in history_rows if row["time"] == item["time"])
        errors = {"maximum": driver.relative(values[0], endpoint["vorticity_dense"]),
                  "ratio": driver.relative(ratios[0], endpoint["vorticity_dense_ratio"])}
        if errors != item["numpy_vs_fftw_relative_errors"] or max(errors.values()) > 1e-10 or driver.relative(baselines[0], history_rows[0]["vorticity_dense"]) > 1e-10:
            raise ValueError("Independent spatial sampling crosscheck failed")
        checkpoint = next(p for p in original["checkpoints"] if p["time"] == item["time"])
        budget = next(p for p in original["budgets"] if p["time"] == item["time"])
        if item["cutoff"] != budget["cutoff"] or item["viscosity"] != original["settings"]["viscosity"]:
            raise ValueError("Sampling metadata differs from the evolved state")
        if item["checkpoint_sha256"] != checkpoint["uncompressed_sha256"]:
            raise ValueError("Samples belong to a different state")
        if item["raw_relative_shift"] > limits["maximum_raw_vorticity_sampling_shift"] or item["ratio_relative_shift"] > limits["maximum_vorticity_amplification_ratio_sampling_shift"]:
            sampling_failures.append("N=%d at T=%g: endpoint spatial sampling shift" % (item["grid"], item["time"]))
    for time in (0.08, 0.10):
        selected = [item for item in recorded["sampling"] if item["time"] == time]
        gap = driver.relative(*(item["amplification_ratios"][1] for item in selected))
        if gap != recorded["cross_grid_vorticity_gap_on_256_samples"][str(time)]:
            raise ValueError("Cross-grid sampling gap differs")
        if gap > limits["maximum_cross_evolution_vorticity_ratio_gap_on_256_samples"]:
            sampling_failures.append("T=%g: cross-evolution vorticity ratio on 256 samples" % time)
    status = ("resolution-or-sampling-gate-failed" if any(p["numerical_failures"] for p in expected_pairs.values()) or sampling_failures else
              "finite-growth-stalled" if any(p["failures"] for p in expected_pairs.values()) else
              "passed-preliminary-finite-amplification-resolution-gates")
    if recorded["status"] != status or recorded["sampling_failures"] != sampling_failures:
        raise ValueError("Scientific verdict differs from the declared gates")
    report = {"status": "archive-integrity-passed", "scientific_status": recorded["status"],
              "initial_state_sha256": INITIAL_SHA, "checkpoints": [],
              "scope": "Archive identity and declared empirical outcomes; an integrity pass does not convert a scientific gate failure into a pass"}
    manifest, rows = fine
    if [b["time"] for b in manifest["budgets"]] != manifest["settings"]["endpoints"]:
        raise ValueError("Fine critical budgets incomplete")
    for path, digest in manifest["source_sha256"].items():
        if driver.sha(RESEARCH / path) != digest:
            raise ValueError("Fine solver source changed: " + path)
    for run, budget in zip(manifest["runs"], manifest["budgets"]):
        log = ROOT / "reference128" / ("t%03d.log" % round(run["endpoint"] * 1000))
        if driver.sha(log) != run["log_sha256"] or run["checkpoint_sha256"] != budget["checkpoint_sha256"]:
            raise ValueError("Fine run identity changed")
        if max(budget["numpy_vs_csv_relative_errors"].values()) > 1e-10:
            raise ValueError("Snapshot budget crosscheck failed")
    for metadata in manifest["checkpoints"]:
        root = ROOT / "reference128"
        index_path = root / ("t%03d_checkpoint.json" % round(metadata["time"] * 1000))
        index = json.loads(index_path.read_text())
        if index["compressed_sha256"] != metadata["sha256"] or index["uncompressed_sha256"] != metadata["uncompressed_sha256"] or index["uncompressed_bytes"] != metadata["uncompressed_bytes"]:
            raise ValueError("Archive index differs from the original checkpoint manifest")
        pieces = []
        for part in index["parts"]:
            path = (root / part["path"]).resolve()
            path.relative_to(root.resolve())
            data = path.read_bytes()
            if len(data) != part["bytes"] or hashlib.sha256(data).hexdigest() != part["sha256"]:
                raise ValueError("Checkpoint part changed")
            pieces.append(data)
        packed = b"".join(pieces)
        if len(packed) != index["compressed_bytes"] or hashlib.sha256(packed).hexdigest() != metadata["sha256"]:
            raise ValueError("Compressed checkpoint identity changed")
        raw = gzip.decompress(packed)
        if len(raw) != metadata["uncompressed_bytes"] or hashlib.sha256(raw).hexdigest() != metadata["uncompressed_sha256"]:
            raise ValueError("Restored checkpoint identity changed")
        grid = struct.unpack_from("<Q", raw, 24)[0]
        observed = struct.unpack_from("<19d", raw, 336)
        row = next(row for row in rows if row["time"] == metadata["time"])
        if raw[:8] != b"NSCONT1\n" or grid != 128 or observed[0] != row["time"] or observed[4] != row["h_half"]:
            raise ValueError("Checkpoint metadata and trajectory disagree")
        if len(raw) != 496 + 48 * grid**3:
            raise ValueError("Checkpoint numerical schema length differs")
        report["checkpoints"].append({"index": str(index_path.relative_to(ROOT)), "grid": grid, "time": row["time"],
            "restored_bytes": len(raw), "restored_sha256": metadata["uncompressed_sha256"], "parts": len(index["parts"])})
    if {item["time"] for item in report["checkpoints"]} != {0.08, 0.10}:
        raise ValueError("A holdout restart state is missing")
    report["maximum_snapshot_diagnostic_relative_error"] = max(max(b["numpy_vs_csv_relative_errors"].values()) for b in manifest["budgets"])
    report["files"] = [{"path": str(p.relative_to(ROOT)), "bytes": p.stat().st_size, "sha256": driver.sha(p)}
        for p in sorted(ROOT.rglob("*")) if p.is_file() and p.name != "audit.json" and
        "__pycache__" not in p.parts and not p.name.endswith((".tmp", ".chk", ".chk.gz"))]
    save(ROOT / "audit.json", report)
    print(json.dumps({key: value for key, value in report.items() if key != "files"}, indent=2))


if __name__ == "__main__":
    main()
