#!/usr/bin/env python3
"""Assess completed evidence without changing the frozen campaign's thresholds."""
import argparse
import json
from pathlib import Path

from run import driver, save


def history(folder):
    manifest = json.loads((folder / "manifest.json").read_text())
    if manifest["status"] != "completed-trajectory-awaiting-crosschecks":
        raise ValueError("Trajectory incomplete: " + str(folder))
    rows = []
    for run in manifest["runs"]:
        path = folder / ("t%03d.csv" % round(run["endpoint"] * 1000))
        if run["returncode"] or driver.sha(path) != run["evidence_sha256"]:
            raise ValueError("Run failed or evidence changed")
        segment = driver.read_evidence(path)
        if rows and rows[-1] != segment[0]:
            raise ValueError("Restart boundary changed")
        rows.extend(segment[1:] if rows else segment)
    return manifest, rows


def paired(a, b, endpoint):
    ma, ra = a
    mb, rb = b
    if ma["initial_state_sha256"] != mb["initial_state_sha256"]:
        raise ValueError("Initial fields differ")
    ra = [row for row in ra if row["time"] <= endpoint + 1e-12]
    rb = [row for row in rb if row["time"] <= endpoint + 1e-12]
    result = driver.assess_pair(ra, rb, endpoint)
    late_a = [row for row in ma["budgets"] if endpoint / 2 - 1e-12 <= row["time"] <= endpoint + 1e-12]
    late_b = [row for row in mb["budgets"] if endpoint / 2 - 1e-12 <= row["time"] <= endpoint + 1e-12]
    if [row["time"] for row in late_a] != [row["time"] for row in late_b] or not late_a:
        raise ValueError("Critical budget clocks differ or are empty")
    key = "h_half_logarithmic_rate"
    result["late_critical_rates"] = [[row[key] for row in rows] for rows in (late_a, late_b)]
    result["late_critical_rate_times"] = [row["time"] for row in late_a]
    result["horizon_times_maximum_gamma_gap"] = endpoint * max(abs(x[key] - y[key]) for x, y in zip(late_a, late_b))
    if result["horizon_times_maximum_gamma_gap"] > 0.002:
        result["failures"].append("cross-resolution critical growth rate")
        result["numerical_failures"].append("cross-resolution critical growth rate")
    if min(row[key] for rows in (late_a, late_b) for row in rows) <= 0:
        result["failures"].append("nonpositive sampled late critical growth rate")
    if result["sampling_shift"]["vorticity_sample"] > 0.02:
        result["failures"].append("physical vorticity sampling shift")
        result["numerical_failures"].append("physical vorticity sampling shift")
    result["status"] = ("resolution-or-run-failure" if result["numerical_failures"] else
                        "finite-growth-stalled" if result["failures"] else "continue-finite-amplification")
    return result


def timestep(reference, refined, endpoint):
    ma, ra = reference
    mb, rb = refined
    by_time = {round(row["time"], 12): row for row in rb}
    common = [(row, by_time[round(row["time"], 12)]) for row in ra if row["time"] <= endpoint + 1e-12]
    keys = ("h_half_ratio", "l3_dense_ratio", "vorticity_dense_ratio", "enstrophy_ratio", "k_rms_ratio")
    gaps = {key: max(driver.relative(a[key], b[key]) for a, b in common) for key in keys}
    reference_steps, refined_steps = (row["step"] for row in common[-1])
    rate_a = {round(row["time"], 12): row["h_half_logarithmic_rate"] for row in ma["budgets"]}
    rate_gaps = [abs(row["h_half_logarithmic_rate"] - rate_a[round(row["time"], 12)]) for row in mb["budgets"]
                 if endpoint / 2 - 1e-12 <= row["time"] <= endpoint + 1e-12]
    failures = []
    if max(gaps.values()) > 0.001:
        failures.append("timestep diagnostic change")
    if refined_steps / reference_steps < 1.5:
        failures.append("insufficient actual timestep refinement")
    if endpoint * max(rate_gaps) > 0.001:
        failures.append("timestep critical growth-rate change")
    late = [row for row in rb if endpoint - 0.02 - 1e-12 <= row["time"] <= endpoint + 1e-12]
    if any(b["h_half"] <= a["h_half"] for a, b in zip(late, late[1:])):
        failures.append("sampled late H growth reversal under denser time observations")
    return {"status": "passed" if not failures else "failed", "failures": failures,
            "maximum_relative_gaps_at_common_times": gaps,
            "reference_steps": reference_steps, "refined_steps": refined_steps,
            "step_count_factor": refined_steps / reference_steps,
            "horizon_times_maximum_gamma_gap": endpoint * max(rate_gaps),
            "refined_observation_interval": mb["settings"]["observation_interval"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True)
    args = parser.parse_args()
    root = Path(args.root)
    data = {name: history(root / name) for name in ("reference32", "reference64", "reference128", "refined64")
            if (root / name / "manifest.json").exists() and
            json.loads((root / name / "manifest.json").read_text())["status"] == "completed-trajectory-awaiting-crosschecks"}
    result = {"scope": "Empirical finite-amplification comparisons; absent/incomplete legs are not passes", "pairs": {}, "timestep": {}}
    for a, b in ((32, 64), (64, 128)):
        if "reference%d" % a in data and "reference%d" % b in data:
            result["pairs"]["%d/%d" % (a, b)] = {str(t): paired(data["reference%d" % a], data["reference%d" % b], t)
                                                for t in (0.08, 0.10)}
    if "reference64" in data and "refined64" in data:
        result["timestep"] = {str(t): timestep(data["reference64"], data["refined64"], t) for t in (0.08, 0.10)}
    save(root / "assessment.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
