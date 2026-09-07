#!/usr/bin/env python3
"""Compare frozen 64/128 trajectories and endpoint samples without refitting."""
import json
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
PREVIOUS = ROOT.parent / "late_growth_continuation"
sys.path.insert(0, str(PREVIOUS))
import assess as comparison
from run import INITIAL, INITIAL_SHA, RESEARCH, analysis, driver, module, save

sampler = module("vorticity_sampler", ROOT.parent / "frozen_continuations_resolution/sample_checkpoint_maximum.py")
restorer = module("checkpoint_restorer", ROOT.parent / "frozen_continuations_resolution/restore_checkpoint.py")


def checkpoint_file(folder, metadata, temporary):
    path = folder / metadata["path"]
    if path.exists():
        if driver.sha(path) != metadata["sha256"]:
            raise ValueError("Compressed checkpoint changed")
        return path
    index = folder / ("t%03d_checkpoint.json" % round(metadata["time"] * 1000))
    output = temporary / (folder.name + "-t%03d.chk" % round(metadata["time"] * 1000))
    restored = restorer.restore(index, output)
    if restored["uncompressed_sha256"] != metadata["uncompressed_sha256"]:
        raise ValueError("Archive index points at another evolved field")
    return output


def main():
    if driver.sha(INITIAL) != INITIAL_SHA:
        raise ValueError("Initial coefficients changed")
    protocol = json.loads((ROOT / "protocol.json").read_text())
    coarse = comparison.history(PREVIOUS / "reference64")
    fine = comparison.history(ROOT / "reference128")
    pairs = {str(t): comparison.paired(coarse, fine, t) for t in (0.08, 0.10)}
    report = {
        "status": "analysis-in-progress", "initial_state_sha256": INITIAL_SHA,
        "protocol_sha256": driver.sha(ROOT / "protocol.json"),
        "analysis_sha256": driver.sha(Path(__file__)),
        "comparison_source_sha256": driver.sha(PREVIOUS / "assess.py"),
        "sampler_source_sha256": driver.sha(Path(sampler.__file__)),
        "checkpoint_reader_sha256": driver.sha(Path(analysis.__file__)),
        "coarse_manifest_sha256": driver.sha(PREVIOUS / "reference64/manifest.json"),
        "fine_manifest_sha256": driver.sha(ROOT / "reference128/manifest.json"),
        "pairs": pairs, "sampling": [], "sampling_failures": [],
        "scope": "Separate finite-amplification resolution check; endpoint samples are not certified continuous maxima",
    }
    save(ROOT / "assessment.json", report)
    # The identical initial field supplies a common physical normalization.
    initial = analysis.load_initial(INITIAL, 64)
    baselines = [sampler.maximum(initial, n) for n in (128, 256)]
    for manifest, rows in (coarse, fine):
        if manifest["initial_state_sha256"] != INITIAL_SHA or driver.relative(rows[0]["vorticity_dense"], baselines[0]) > 1e-10:
            raise ValueError("Different initial field or physical normalization")
    settings = protocol["additional_sampling"]
    with tempfile.TemporaryDirectory(prefix="sampling-", dir=ROOT) as directory:
        temporary = Path(directory)
        for grid, folder, history in ((64, PREVIOUS / "reference64", coarse), (128, ROOT / "reference128", fine)):
            manifest, rows = history
            for metadata in manifest["checkpoints"]:
                path = checkpoint_file(folder, metadata, temporary)
                field, viscosity, observed, layout = analysis.load_checkpoint(path, metadata["uncompressed_sha256"])
                values = [sampler.maximum(field, n) for n in (128, 256)]
                endpoint = next(row for row in rows if row["time"] == metadata["time"])
                ratios = [value / baseline for value, baseline in zip(values, baselines)]
                errors = {"maximum": driver.relative(values[0], endpoint["vorticity_dense"]),
                          "ratio": driver.relative(ratios[0], endpoint["vorticity_dense_ratio"])}
                if max(errors.values()) > 1e-10 or values[1] < values[0] * (1 - 1e-12):
                    raise ValueError("Independent or nested spatial samples disagree")
                raw_shift, ratio_shift = driver.relative(*values), driver.relative(*ratios)
                item = {"grid": grid, "cutoff": layout["cutoff"], "time": observed[0], "viscosity": viscosity,
                        "sampling_grids": [128, 256], "initial_maxima": baselines, "final_maxima": values,
                        "amplification_ratios": ratios, "raw_relative_shift": raw_shift,
                        "ratio_relative_shift": ratio_shift, "numpy_vs_fftw_relative_errors": errors,
                        "checkpoint_sha256": metadata["uncompressed_sha256"]}
                report["sampling"].append(item)
                if raw_shift > settings["maximum_raw_vorticity_sampling_shift"] or ratio_shift > settings["maximum_vorticity_amplification_ratio_sampling_shift"]:
                    report["sampling_failures"].append("N=%d at T=%g: endpoint spatial sampling shift" % (grid, observed[0]))
                save(ROOT / "assessment.json", report)
                print("Samples", grid, observed[0], "ratios", ratios, "shift", ratio_shift, flush=True)
                del field
    report["cross_grid_vorticity_gap_on_256_samples"] = {}
    for time in (0.08, 0.10):
        selected = [row for row in report["sampling"] if row["time"] == time]
        if len(selected) != 2 or {row["grid"] for row in selected} != {64, 128}:
            raise ValueError("Missing endpoint sampling leg")
        gap = driver.relative(*(row["amplification_ratios"][1] for row in selected))
        report["cross_grid_vorticity_gap_on_256_samples"][str(time)] = gap
        if gap > settings["maximum_cross_evolution_vorticity_ratio_gap_on_256_samples"]:
            report["sampling_failures"].append("T=%g: cross-evolution vorticity ratio on 256 samples" % time)
    if any(row["numerical_failures"] for row in pairs.values()) or report["sampling_failures"]:
        report["status"] = "resolution-or-sampling-gate-failed"
    elif any(row["failures"] for row in pairs.values()):
        report["status"] = "finite-growth-stalled"
    else:
        report["status"] = "passed-preliminary-finite-amplification-resolution-gates"
    save(ROOT / "assessment.json", report)
    print("Resolution follow-up:", report["status"], flush=True)


if __name__ == "__main__":
    main()
