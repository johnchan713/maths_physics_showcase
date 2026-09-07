#!/usr/bin/env python3
"""Check archive identity separately from the finite-growth scientific verdict."""
import gzip
import hashlib
import json
import math
import struct

from analyze import (ROOT, INITIAL_SHA, RESEARCH, analysis, comparison, driver,
                     history, sampler, save, scientific_status)
from archive import verify_parts


def main():
    protocol = json.loads((ROOT / "protocol.json").read_text())
    recorded = json.loads((ROOT / "assessment.json").read_text())
    if recorded["status"] == "analysis-in-progress":
        raise ValueError("Sampling assessment incomplete")
    if recorded["initial_state_sha256"] != INITIAL_SHA:
        raise ValueError("Different initial coefficients")
    sources = {"protocol_sha256": ROOT / "protocol.json", "analysis_sha256": ROOT / "analyze.py",
               "sampler_source_sha256": sampler.__file__, "checkpoint_reader_sha256": analysis.__file__}
    for key, path in sources.items():
        if recorded[key] != driver.sha(path):
            raise ValueError("Assessment dependency changed: " + key)
    data = {grid: history(grid) for grid in protocol["evolution_grids"]}
    expected_pairs = {str(t): comparison.paired(data[64][:2], data[128][:2], t) for t in protocol["endpoints"]}
    if expected_pairs != recorded["pairs"]:
        raise ValueError("Recorded assessment differs from the archived histories")
    report = {"status": "archive-integrity-passed", "scientific_status": recorded["status"],
              "initial_state_sha256": INITIAL_SHA, "checkpoints": [],
              "scope": "Identity and declared empirical outcomes; an archive-integrity pass does not override a scientific failure"}
    for grid, (_, rows, manifest) in data.items():
        folder = ROOT / ("reference%d" % grid)
        if driver.sha(folder / "manifest.json") != recorded["extension_manifests"][str(grid)]:
            raise ValueError("Extension manifest changed")
        restart = manifest["restart_verification"]
        if restart["time"] != protocol["start_time"] or restart["checkpoint_sha256"] != manifest["source_history"]["checkpoint_sha256"]:
            raise ValueError("Restored .10 field identity differs")
        for budget in [restart] + manifest["budgets"]:
            if max(budget["numpy_vs_csv_relative_errors"].values()) > 1e-10:
                raise ValueError("Independent critical-budget check failed")
            if budget["grid"] != grid or budget["cutoff"] != (grid - 1) // 3:
                raise ValueError("Budget is from a different grid")
        if [c["time"] for c in manifest["checkpoints"]] != protocol["endpoints"]:
            raise ValueError("A declared restart target is missing")
        for metadata in manifest["checkpoints"]:
            index_path = folder / ("t%03d_checkpoint.json" % round(metadata["time"] * 1000))
            index = json.loads(index_path.read_text())
            if index["grid"] != grid or index["time"] != metadata["time"] or index["initial_state_sha256"] != INITIAL_SHA:
                raise ValueError("Archive metadata differs")
            if index["compressed_sha256"] != metadata["sha256"] or index["uncompressed_sha256"] != metadata["uncompressed_sha256"] or index["uncompressed_bytes"] != metadata["uncompressed_bytes"]:
                raise ValueError("Archive index differs from the evolved-state manifest")
            verify_parts(folder, index)
            packed = b"".join((folder / part["path"]).read_bytes() for part in index["parts"])
            raw = gzip.decompress(packed)
            if len(raw) != metadata["uncompressed_bytes"] or hashlib.sha256(raw).hexdigest() != metadata["uncompressed_sha256"]:
                raise ValueError("Restored evolved state differs")
            n, cutoff = struct.unpack_from("<2Q", raw, 24)
            observed = struct.unpack_from("<19d", raw, 336)
            row = next(row for row in rows if row["time"] == metadata["time"])
            if raw[:8] != b"NSCONT1\n" or n != grid or cutoff != index["cutoff"] or observed[0] != row["time"] or observed[4] != row["h_half"]:
                raise ValueError("Restored checkpoint and evidence disagree")
            if len(raw) != 496 + 48 * grid**3:
                raise ValueError("Restored checkpoint numerical layout differs")
            report["checkpoints"].append({"index": str(index_path.relative_to(ROOT)), "grid": grid,
                "time": row["time"], "restored_bytes": len(raw), "restored_sha256": metadata["uncompressed_sha256"],
                "parts": len(index["parts"])})
    expected_cases = {(n, t) for n in protocol["evolution_grids"] for t in protocol["endpoints"]}
    if len(recorded["sampling"]) != len(expected_cases) or {(x["grid"], x["time"]) for x in recorded["sampling"]} != expected_cases:
        raise ValueError("Sampling cases missing or duplicated")
    limits = protocol["additional_sampling"]
    shared_baselines = recorded["sampling"][0]["initial_maxima"]
    failures = []
    for item in recorded["sampling"]:
        values, baselines = item["final_maxima"], item["initial_maxima"]
        if item["sampling_grids"] != limits["sample_grids"] or len(values) != 2 or len(baselines) != 2 or baselines != shared_baselines:
            raise ValueError("Sampling grids or shared normalization differ")
        if any(not math.isfinite(v) or v <= 0 for v in values + baselines):
            raise ValueError("Invalid spatial maximum")
        if values[1] < values[0] * (1 - 1e-12) or baselines[1] < baselines[0] * (1 - 1e-12):
            raise ValueError("Nested sampling lost a retained maximum")
        ratios = [v / b for v, b in zip(values, baselines)]
        if ratios != item["amplification_ratios"] or driver.relative(*values) != item["raw_relative_shift"] or driver.relative(*ratios) != item["ratio_relative_shift"]:
            raise ValueError("Sampling arithmetic differs")
        _, rows, manifest = data[item["grid"]]
        endpoint = next(row for row in rows if row["time"] == item["time"])
        errors = {"maximum": driver.relative(values[0], endpoint["vorticity_dense"]),
                  "ratio": driver.relative(ratios[0], endpoint["vorticity_dense_ratio"])}
        if errors != item["numpy_vs_fftw_relative_errors"] or max(errors.values()) > 1e-10 or driver.relative(baselines[0], rows[0]["vorticity_dense"]) > 1e-10:
            raise ValueError("Independent physical sampling check failed")
        checkpoint = next(p for p in manifest["checkpoints"] if p["time"] == item["time"])
        budget = next(p for p in manifest["budgets"] if p["time"] == item["time"])
        if item["checkpoint_sha256"] != checkpoint["uncompressed_sha256"] or item["cutoff"] != budget["cutoff"] or item["viscosity"] != protocol["viscosity"]:
            raise ValueError("Sampling field identity differs")
        if item["raw_relative_shift"] > limits["maximum_raw_vorticity_sampling_shift"] or item["ratio_relative_shift"] > limits["maximum_vorticity_amplification_ratio_sampling_shift"]:
            failures.append("N=%d at T=%g: endpoint spatial sampling shift" % (item["grid"], item["time"]))
    for time in protocol["endpoints"]:
        selected = [p for p in recorded["sampling"] if p["time"] == time]
        gap = driver.relative(*(p["amplification_ratios"][1] for p in selected))
        if gap != recorded["cross_grid_vorticity_gap_on_256_samples"][str(time)]:
            raise ValueError("Cross-grid spatial sampling gap differs")
        if gap > limits["maximum_cross_evolution_vorticity_ratio_gap_on_256_samples"]:
            failures.append("T=%g: cross-evolution vorticity ratio on 256 samples" % time)
    status = scientific_status(expected_pairs, failures)
    if status != recorded["status"] or failures != recorded["sampling_failures"]:
        raise ValueError("Recorded scientific verdict differs from the declared gates")
    report["maximum_snapshot_diagnostic_relative_error"] = max(
        max(b["numpy_vs_csv_relative_errors"].values()) for _, _, m in data.values()
        for b in [m["restart_verification"]] + m["budgets"])
    report["files"] = [{"path": str(p.relative_to(ROOT)), "bytes": p.stat().st_size, "sha256": driver.sha(p)}
        for p in sorted(ROOT.rglob("*")) if p.is_file() and p.name != "audit.json" and
        "__pycache__" not in p.parts and not p.name.endswith((".tmp", ".chk", ".chk.gz"))]
    save(ROOT / "audit.json", report)
    print(json.dumps({k: v for k, v in report.items() if k != "files"}, indent=2))


if __name__ == "__main__":
    main()
