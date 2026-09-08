#!/usr/bin/env python3
"""Verify input identity, spectral arithmetic and recorded vectors by direct sums.

The audit recomputes spectral norms and Fourier values at recorded locations.
It does not repeat every 256-grid maximum search or certify continuous maxima.
"""
import argparse
import json
import math
from pathlib import Path
import tempfile

import numpy as np

from analyze import (ROOT, RESEARCH, direct_curl, embed, frequencies, low_pass, projection_fractions,
                     reader, relative, restorer, save, sha, spectral_norms)


def close(a, b, tolerance=1e-10, absolute_tolerance=0.):
    if any(not math.isfinite(v) for v in (a, b, tolerance, absolute_tolerance)) or min(tolerance, absolute_tolerance) < 0:
        raise ValueError("Comparison needs finite values and nonnegative tolerances")
    if relative(a, b) > tolerance and abs(a - b) > absolute_tolerance:
        raise ValueError("Recorded arithmetic differs: %r versus %r" % (a, b))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "audit.json")
    args = parser.parse_args()
    if args.output.exists() and args.output.resolve() != (ROOT / "audit.json").resolve():
        raise ValueError("Refusing to overwrite a noncanonical audit output")
    protocol = json.loads((ROOT / "protocol.json").read_text())
    diagnosis = json.loads((ROOT / "diagnosis.json").read_text())
    probe = json.loads((ROOT / "probe.json").read_text())
    controls = json.loads((ROOT / "probe_controls.json").read_text())
    build = json.loads((ROOT / "build.json").read_text())
    if build["source_sha256"] != sha(ROOT / "benchmark.cpp") or build["binary_sha256"] != probe["binary_sha256"]:
        raise ValueError("Probe build provenance differs")
    for name, digest in protocol["source_sha256"].items():
        if sha(RESEARCH / name) != digest:
            raise ValueError("Frozen dependency changed: " + name)
    initial_sha = sha(RESEARCH / protocol["initial_state"])
    if initial_sha != protocol["initial_state_sha256"]:
        raise ValueError("Initial field changed")
    if diagnosis["status"] != "diagnosis-complete-original-resolution-failure-retained":
        raise ValueError("Diagnosis incomplete or improperly promoted")
    if diagnosis["analysis_sha256"] != sha(ROOT / "analyze.py"):
        raise ValueError("Diagnosis source changed")
    if any(item["protocol_sha256"] != sha(ROOT / "protocol.json") or
           item["initial_state_sha256"] != initial_sha for item in (diagnosis, probe)):
        raise ValueError("Protocol or initial identity differs")
    if probe["settings"] != protocol["cost_probe"] or probe["status"] != "short-probe-passed" or probe["returncode"] != 0:
        raise ValueError("Capacity probe did not pass its bounded specification")
    if probe["stdout_sha256"] != sha(ROOT / "probe.stdout") or probe["log_sha256"] != sha(ROOT / "probe.log"):
        raise ValueError("Capacity-probe raw evidence changed")
    if probe["measurements"] != json.loads((ROOT / "probe.stdout").read_text()):
        raise ValueError("Capacity-probe measurements differ from stdout")
    measurements = probe["measurements"]
    if measurements["relative_state_l2_gap"] > protocol["cost_probe"]["maximum_relative_state_l2_gap"]:
        raise ValueError("Short-probe agreement failed")
    if [p["grid"] for p in measurements["probes"]] != protocol["cost_probe"]["grids"]:
        raise ValueError("Cost-probe grids differ")
    for p, cutoff in zip(measurements["probes"], protocol["cost_probe"]["cutoffs"]):
        if p["steps"] != 2 or p["final_time"] != .00002 or p["cutoff"] != cutoff or len(p["rk4_seconds"]) != 2:
            raise ValueError("Cost-probe clock or cutoff differs")
        if any(not math.isfinite(t) or t <= 0 for t in p["rk4_seconds"]):
            raise ValueError("Invalid step timing")
        close(p["initial_energy"], 10.)
        if p["final_energy"] > p["initial_energy"] * (1 + 1e-10) or max(p["divergence_defect"], p["reality_defect"]) > 1e-9:
            raise ValueError("Probe invariant failure")
    if controls["status"] != "passed" or controls["binary_sha256"] != probe["binary_sha256"] or controls["source_sha256"] != sha(ROOT / "test_probe.py"):
        raise ValueError("Probe controls are from different code")
    if [c["returncode"] for c in controls["controls"]] != [0, 1, 1]:
        raise ValueError("Control results changed")
    if controls["controls"][0]["measurements"]["relative_state_l2_gap"] > 1e-10 or controls["controls"][1]["measurements"]["relative_state_l2_gap"] <= 1e-10:
        raise ValueError("Positive or negative agreement control is ineffective")
    settings = protocol["snapshot_audit"]
    if [c["time"] for c in diagnosis["cases"]] != settings["times"]:
        raise ValueError("Snapshot clock changed")
    report = {"status": "input-and-diagnostic-audit-passed",
              "scientific_status": "original-t014-resolution-failure-retained",
              "capacity_scope": "Two steps only; no completed .14 grid-256 evolution",
              "verification_scope": __doc__, "cases": [], "prior_assessment_sha256": {}}
    with tempfile.TemporaryDirectory(prefix="audit-inputs-", dir=ROOT) as directory:
        for case in diagnosis["cases"]:
            fields = {}
            for item in [c for c in protocol["checkpoints"] if c["time"] == case["time"]]:
                if case["sources"][str(item["grid"])] != item:
                    raise ValueError("Diagnosis source link changed")
                index = RESEARCH / item["index"]
                if sha(index) != item["index_sha256"]:
                    raise ValueError("Input index changed")
                path = Path(directory) / ("n%d-t%03d.chk" % (item["grid"], round(item["time"] * 1000)))
                restorer.restore(index, path)
                field, nu, observed, layout = reader.load_checkpoint(path, item["state_sha256"])
                if layout["grid"] != item["grid"] or nu != .02 or observed[0] != case["time"]:
                    raise ValueError("Input state metadata changed")
                fields[item["grid"]] = field
                path.unlink()
            coarse, fine = fields[64], fields[128]
            lifted = embed(coarse, 128)
            low = low_pass(fine, settings["shared_cutoff"])
            shared, extra = low - lifted, fine - low
            actual = {"coarse": coarse, "fine": fine, "shared_error": shared, "extra": extra,
                      "total_error": fine - lifted}
            for name, coefficients in actual.items():
                recomputed = spectral_norms(coefficients)
                for key, value in recomputed.items():
                    close(value, case["norms"][name][key])
                f = frequencies(coefficients.shape[0])
                k = (f[:, None, None], f[None, :, None], f[None, None, :])
                curl_squared = 0.
                for c in range(3):
                    a, b = (c + 1) % 3, (c + 2) % 3
                    curl_squared += float(np.sum(np.abs(k[a] * coefficients[..., b] - k[b] * coefficients[..., a])**2))
                close(curl_squared, recomputed["vorticity_l2_squared"])
            for key in case["norms"]["fine"]:
                close(case["norms"]["total_error"][key], case["norms"]["shared_error"][key] + case["norms"]["extra"][key])
                close(case["high_mode_fractions_of_fine_squared_norms"][key], case["norms"]["extra"][key] / case["norms"]["fine"][key])
                close(case["shared_mode_fractions_of_squared_error"][key], case["norms"]["shared_error"][key] / case["norms"]["total_error"][key])
            points = [p["index"] for p in case["pointwise"]]
            if any(len(p) != 3 or any(not isinstance(i, int) or not 0 <= i < settings["sampling_grid"] for i in p) for p in points):
                raise ValueError("Invalid sample location")
            direct_error = 0.
            scale = case["peaks"]["fine"]["maximum"]
            for name in ("coarse", "fine", "shared_error", "extra"):
                values = direct_curl(actual[name], points, settings["sampling_grid"])
                recorded = np.array([p["vectors"][name] for p in case["pointwise"]])
                direct_error = max(direct_error, float(np.max(np.abs(values - recorded))) / scale)
            if direct_error > settings["identity_tolerance"]:
                raise ValueError("Direct Fourier values disagree with the recorded vectors")
            for point in case["pointwise"]:
                vectors = point["vectors"]
                np.testing.assert_allclose(np.array(vectors["fine"]) - vectors["coarse"],
                                           np.array(vectors["shared_error"]) + vectors["extra"], rtol=0., atol=scale * 1e-10)
                fractions = projection_fractions(vectors["shared_error"], vectors["extra"])
                recorded = point["projected_error_fractions"]
                if fractions is None or recorded is None:
                    if fractions != recorded:
                        raise ValueError("Signed pointwise projection changed")
                else:
                    for name in ("shared", "extra"):
                        close(fractions[name], recorded[name])
            close(float(np.linalg.norm(case["pointwise"][0]["vectors"]["coarse"])), case["peaks"]["coarse"]["maximum"])
            close(float(np.linalg.norm(case["pointwise"][1]["vectors"]["fine"])), scale)
            difference = np.array(case["pointwise"][2]["vectors"]["fine"]) - case["pointwise"][2]["vectors"]["coarse"]
            close(float(np.linalg.norm(difference)), case["peaks"]["total_error"]["maximum"])
            prior_folder = "late_growth_t012" if case["time"] == .12 else "late_growth_t014"
            prior_path = ROOT.parent / prior_folder / "assessment.json"
            prior = json.loads(prior_path.read_text())
            report["prior_assessment_sha256"][str(prior_path.relative_to(RESEARCH))] = sha(prior_path)
            for n, name in ((64, "coarse"), (128, "fine")):
                saved = next(v for v in prior["sampling"] if v["grid"] == n and v["time"] == case["time"])
                close(case["peaks"][name]["maximum"], saved["final_maxima"][1])
                close(case["vorticity_amplification_ratios"][name], saved["amplification_ratios"][1])
            gap = relative(case["peaks"]["coarse"]["maximum"], scale)
            close(gap, case["sampled_vorticity_relative_gap"])
            if case["original_endpoint_vorticity_gate_failed"] != (gap > .10):
                raise ValueError("Original gate result changed")
            for filtered in case["filtered_fine"]:
                coefficients = low_pass(fine, filtered["cutoff"])
                values = direct_curl(coefficients, [filtered["index"]], settings["sampling_grid"])
                close(float(np.linalg.norm(values[0])), filtered["maximum"])
                close(relative(filtered["maximum"], scale), filtered["relative_maximum_gap_to_full"])
                removed = max(0., 1 - spectral_norms(coefficients)["velocity_l2_squared"] / case["norms"]["fine"]["velocity_l2_squared"])
                # A 1 - E_filtered/E_full subtraction loses relative accuracy
                # for tiny removed fractions. Permit only unit-scale roundoff;
                # this is arithmetic verification, not a scientific gate change.
                close(removed, filtered["removed_energy_fraction"],
                      absolute_tolerance=64 * np.finfo(float).eps)
            report["cases"].append({"time": case["time"], "direct_fourier_relative_error": direct_error,
                                    "scalar_peak_gap": gap,
                                    "maximum_sampled_vector_difference_over_fine_peak": case["peaks"]["total_error"]["maximum"] / scale})
            del fields, actual, coarse, fine, lifted, low, shared, extra, field, coefficients
    if not diagnosis["cases"][-1]["original_endpoint_vorticity_gate_failed"]:
        raise ValueError("The original .14 failure must be retained")
    step_count_model = 2749 * 85 / 42
    report["cost_projection"] = {"model": "Prior 128-grid step count times cutoff ratio, assuming comparable velocity bounds; not a completed run or confidence interval",
        "modeled_steps": step_count_model,
        "stepping_hours_from_each_probe_step": [t * step_count_model / 3600 for t in measurements["probes"][1]["rk4_seconds"]],
        "max_dt_step_count": 1120,
        "max_dt_only_scenario_hours": [t * 1120 / 3600 for t in measurements["probes"][1]["rk4_seconds"]],
        "excluded": "Initialization, observations, checkpoint I/O, changes in machine load, and differences in the actual fine trajectory"}
    report["files"] = [{"path": str(p.relative_to(ROOT)), "sha256": sha(p), "bytes": p.stat().st_size}
                       for p in sorted(ROOT.rglob("*")) if p.is_file() and p.name != "audit.json" and
                       "__pycache__" not in p.parts and not p.name.endswith(".tmp")]
    save(args.output, report)
    print(json.dumps({k: v for k, v in report.items() if k != "files"}, indent=2))


if __name__ == "__main__":
    main()
