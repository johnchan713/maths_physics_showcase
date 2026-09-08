#!/usr/bin/env python3
"""Assess the frozen .14 extension with the unchanged preceding comparison gates."""
import json
from pathlib import Path
import tempfile

from extend import (ROOT, PREVIOUS, INITIAL, INITIAL_SHA, RESEARCH, analysis,
                    comparison, driver, module, provenance, restorer, save, stitch)

sampler = module("extension_sampler", ROOT.parent / "frozen_continuations_resolution/sample_checkpoint_maximum.py")


def history(grid, endpoint=None):
    protocol, source, previous, earlier_rows = provenance(grid)
    final = protocol["endpoints"][-1] if endpoint is None else endpoint
    if final not in protocol["endpoints"]:
        raise ValueError("Undeclared continuation horizon")
    folder = ROOT / ("reference%d" % grid)
    manifest = json.loads((folder / "manifest.json").read_text())
    allowed = ("completed-trajectory-awaiting-crosschecks",) if final == protocol["endpoints"][-1] else (
        "running", "completed-trajectory-awaiting-crosschecks")
    if manifest["status"] not in allowed:
        raise ValueError("Extension incomplete or stopped: " + str(folder))
    if manifest["initial_state_sha256"] != INITIAL_SHA or manifest["protocol_sha256"] != driver.sha(ROOT / "protocol.json"):
        raise ValueError("Extension initial field or protocol changed")
    sources = {"extension_source_sha256": ROOT / "extend.py",
               "restorer_source_sha256": Path(restorer.__file__),
               "comparison_source_sha256": Path(comparison.__file__)}
    for key, path in sources.items():
        if manifest[key] != driver.sha(path):
            raise ValueError("Extension dependency changed: " + key)
    if manifest["source_sha256"] != previous["source_sha256"] or manifest["solver_sha256"] != previous["solver_sha256"]:
        raise ValueError("Extension used different solver provenance")
    parent = manifest["source_history"]
    if parent["path"] != str(source.relative_to(RESEARCH)) or parent["manifest_sha256"] != driver.sha(source / "manifest.json"):
        raise ValueError("Earlier trajectory provenance changed")
    archive = protocol["source_checkpoints"][str(grid)]
    if parent["checkpoint_archive"] != {**archive, "sha256": driver.sha(RESEARCH / archive["path"])}:
        raise ValueError("Earlier checkpoint archive changed")
    original = next(p for p in previous["checkpoints"] if p["time"] == protocol["start_time"])
    if parent["checkpoint_sha256"] != original["uncompressed_sha256"]:
        raise ValueError("Wrong restart field identity")
    settings = {**previous["settings"], "start_time": protocol["start_time"], "endpoints": protocol["endpoints"]}
    if manifest["settings"] != settings:
        raise ValueError("Restored settings changed")
    segments, budgets = [], []
    last_sha = original["uncompressed_sha256"]
    for run in manifest["runs"]:
        if run["endpoint"] > final + 1e-12:
            continue
        stem = "t%03d" % round(run["endpoint"] * 1000)
        path = folder / (stem + ".csv")
        if run.get("returncode") != 0 or run["evidence_sha256"] != driver.sha(path) or run["log_sha256"] != driver.sha(folder / (stem + ".log")):
            raise ValueError("Extension segment failed or changed")
        if run["restart_sha256"] != last_sha:
            raise ValueError("Restart hash chain changed")
        budget = next(b for b in manifest["budgets"] if b["time"] == run["endpoint"])
        if run["checkpoint_sha256"] != budget["checkpoint_sha256"]:
            raise ValueError("Snapshot and evolved checkpoint identities differ")
        last_sha = run["checkpoint_sha256"]
        budgets.append(budget)
        segments.append(driver.read_evidence(path))
    rows = stitch(earlier_rows, segments, final, protocol["observation_interval"])
    expected_times = [t for t in protocol["endpoints"] if t <= final + 1e-12]
    if [b["time"] for b in budgets] != expected_times:
        raise ValueError("Incomplete critical budget clock")
    joined = {**previous, "budgets": previous["budgets"] + budgets}
    return joined, rows, manifest


def scientific_status(pairs, sampling_failures):
    if not pairs:
        raise ValueError("Missing paired trajectory assessments")
    if any(p["numerical_failures"] for p in pairs.values()) or sampling_failures:
        return "resolution-or-sampling-gate-failed"
    if any(p["failures"] for p in pairs.values()):
        return "finite-growth-stalled"
    return "passed-preliminary-finite-amplification-continuation-gates"


def checkpoint_file(folder, metadata, temporary):
    path = folder / metadata["path"]
    if path.exists():
        if driver.sha(path) != metadata["sha256"]:
            raise ValueError("Compressed extension checkpoint changed")
        return path
    output = temporary / (folder.name + "-t%03d.chk" % round(metadata["time"] * 1000))
    restored = restorer.restore(folder / ("t%03d_checkpoint.json" % round(metadata["time"] * 1000)), output)
    if restored["uncompressed_sha256"] != metadata["uncompressed_sha256"]:
        raise ValueError("Parts describe a different evolved state")
    return output


def main():
    protocol = json.loads((ROOT / "protocol.json").read_text())
    data = {grid: history(grid) for grid in protocol["evolution_grids"]}
    report = {
        "status": "analysis-in-progress", "initial_state_sha256": INITIAL_SHA,
        "protocol_sha256": driver.sha(ROOT / "protocol.json"),
        "analysis_sha256": driver.sha(Path(__file__)),
        "sampler_source_sha256": driver.sha(Path(sampler.__file__)),
        "checkpoint_reader_sha256": driver.sha(Path(analysis.__file__)),
        "extension_manifests": {str(n): driver.sha(ROOT / ("reference%d/manifest.json" % n)) for n in data},
        "pairs": {str(t): comparison.paired(data[64][:2], data[128][:2], t) for t in protocol["endpoints"]},
        "sampling": [], "sampling_failures": [],
        "scope": "Frozen-field empirical continuation; no new independent or refined-timestep evolution beyond .10 is claimed",
    }
    save(ROOT / "assessment.json", report)
    settings = protocol["additional_sampling"]
    sample_grids = settings["sample_grids"]
    initial = analysis.load_initial(INITIAL, 64)
    baselines = [sampler.maximum(initial, n) for n in sample_grids]
    with tempfile.TemporaryDirectory(prefix="sampling-", dir=ROOT) as directory:
        temporary = Path(directory)
        for grid, (_, rows, manifest) in data.items():
            if driver.relative(baselines[0], rows[0]["vorticity_dense"]) > 1e-10:
                raise ValueError("Different physical initial normalization")
            folder = ROOT / ("reference%d" % grid)
            for metadata in manifest["checkpoints"]:
                path = checkpoint_file(folder, metadata, temporary)
                field, viscosity, observed, layout = analysis.load_checkpoint(path, metadata["uncompressed_sha256"])
                endpoint = next(row for row in rows if row["time"] == metadata["time"])
                if observed[0] != metadata["time"] or layout["grid"] != grid:
                    raise ValueError("Wrong evolved sampling field")
                values = [sampler.maximum(field, n) for n in sample_grids]
                ratios = [v / b for v, b in zip(values, baselines)]
                errors = {"maximum": driver.relative(values[0], endpoint["vorticity_dense"]),
                          "ratio": driver.relative(ratios[0], endpoint["vorticity_dense_ratio"])}
                if max(errors.values()) > 1e-10 or values[1] < values[0] * (1 - 1e-12):
                    raise ValueError("Independent or nested samples disagree")
                item = {"grid": grid, "cutoff": layout["cutoff"], "time": observed[0], "viscosity": viscosity,
                        "sampling_grids": sample_grids, "initial_maxima": baselines, "final_maxima": values,
                        "amplification_ratios": ratios, "raw_relative_shift": driver.relative(*values),
                        "ratio_relative_shift": driver.relative(*ratios), "numpy_vs_fftw_relative_errors": errors,
                        "checkpoint_sha256": metadata["uncompressed_sha256"]}
                report["sampling"].append(item)
                if item["raw_relative_shift"] > settings["maximum_raw_vorticity_sampling_shift"] or item["ratio_relative_shift"] > settings["maximum_vorticity_amplification_ratio_sampling_shift"]:
                    report["sampling_failures"].append("N=%d at T=%g: endpoint spatial sampling shift" % (grid, observed[0]))
                save(ROOT / "assessment.json", report)
                print("Samples", grid, observed[0], "ratios", ratios, "shift", item["ratio_relative_shift"], flush=True)
                del field
    report["cross_grid_vorticity_gap_on_256_samples"] = {}
    for endpoint in protocol["endpoints"]:
        selected = [row for row in report["sampling"] if row["time"] == endpoint]
        if len(selected) != 2 or {row["grid"] for row in selected} != {64, 128}:
            raise ValueError("Missing sampling leg")
        gap = driver.relative(*(row["amplification_ratios"][1] for row in selected))
        report["cross_grid_vorticity_gap_on_256_samples"][str(endpoint)] = gap
        if gap > settings["maximum_cross_evolution_vorticity_ratio_gap_on_256_samples"]:
            report["sampling_failures"].append("T=%g: cross-evolution vorticity ratio on 256 samples" % endpoint)
    report["status"] = scientific_status(report["pairs"], report["sampling_failures"])
    save(ROOT / "assessment.json", report)
    print("Continuation assessment:", report["status"], flush=True)


if __name__ == "__main__":
    main()
