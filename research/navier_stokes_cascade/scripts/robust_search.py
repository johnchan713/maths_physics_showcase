#!/usr/bin/env python3
"""Bounded discovery and held-out validation of finite-amplification candidates.

Only Python's standard library is needed. Every subprocess, failure, coefficient
checksum, and raw trace is retained. These numerical gates are experimental
triage criteria, not sufficient or necessary conditions for a PDE singularity.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess


FAMILIES = ("wave-packets", "vortex-tubes", "orthogonal-bundle")
LIMITS = {
    "cutoff_margin": 0.008,
    "cross_h_final": 0.02,
    "cross_l3_final": 0.02,
    "cross_vorticity_final": 0.10,
    "cross_enstrophy_final": 0.05,
    "cross_scale_final": 0.05,
    "cross_stretch_late_min": 0.10,
    "minimum_h_final": 1.005,
    "minimum_h_late": 1.001,
    "minimum_scale_final": 1.05,
    "minimum_stretch_late": 1.0,
    "perturbation_relative": 0.001,
    "sampling_vorticity_relative": 0.02,
    "stretch_perturbation_relative": 0.02,
    "gain_to_error_factor": 3.0,
    "cross_late_rate_dimensionless": 0.002,
    "perturbation_late_rate_dimensionless": 0.001,
}
COMPARISONS = ("h_final", "l3_final", "vorticity_final", "enstrophy_final",
               "scale_final", "stretch_late_min")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(left, right):
    return abs(left - right) / max(abs(left), abs(right), 1e-15)


def number(row, key):
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {key}")
    return value


def read_rows(path):
    with Path(path).open(newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
    if not rows or any(None in row or None in row.values() for row in rows):
        raise ValueError(f"missing or malformed CSV rows: {path}")
    return rows


def summarize_evidence(rows, horizon, samples, sampling_grid, viscosity, energy):
    """Require a complete common clock and check the actual exported PDE budget."""
    result = {}
    if len(rows) != 2 * (samples + 1):
        raise ValueError("incomplete or duplicate evidence snapshots")
    for resolution in ("coarse", "fine"):
        series = [r for r in rows if r["resolution"] == resolution]
        if len(series) != samples + 1:
            raise ValueError("unbalanced evidence resolutions")
        for i, row in enumerate(series):
            for key in row:
                if key != "resolution":
                    number(row, key)
            if (int(row["snapshot"]) != i or
                    int(row["sampling_grid"]) != sampling_grid or
                    abs(number(row, "time") - horizon * i / samples) > 1e-12):
                raise ValueError("evidence clock or sampling-grid mismatch")
            for key in ("energy", "enstrophy", "h_half", "l3_sample",
                        "vorticity_sample", "k_rms", "viscous_destruction"):
                if number(row, key) <= 0:
                    raise ValueError(f"non-positive {key}")
            destruction = number(row, "viscous_destruction")
            production = number(row, "stretching")
            if (relative(destruction, 2 * viscosity * number(row, "palinstrophy")) > 1e-10 or
                    abs(number(row, "net_enstrophy_rate") - production + destruction) >
                    1e-10 * max(1, abs(production), destruction) or
                    relative(number(row, "production_to_dissipation"),
                             production / destruction) > 1e-10):
                raise ValueError("inconsistent enstrophy budget")
            if number(row, "vorticity_sample") > number(row, "vorticity_fourier_upper_bound") * (1 + 1e-12):
                raise ValueError("sample exceeds Fourier upper bound")
        initial, middle, final = series[0], series[samples // 2], series[-1]
        if relative(number(initial, "energy"), energy) > 1e-10:
            raise ValueError("initial energy mismatch")
        late = series[samples // 2:]
        metrics = {name: number(final, key) / number(initial, key)
                   for name, key in (("h_final", "h_half"), ("l3_final", "l3_sample"),
                                     ("vorticity_final", "vorticity_sample"),
                                     ("enstrophy_final", "enstrophy"),
                                     ("scale_final", "k_rms"))}
        metrics.update({
            "h_late": number(final, "h_half") / number(middle, "h_half"),
            "h_late_min_step": min(number(b, "h_half") / number(a, "h_half")
                                   for a, b in zip(late, late[1:])),
            "stretch_late_min": min(number(r, "production_to_dissipation") for r in late),
            "vorticity_peak_sampled": max(number(r, "vorticity_sample") for r in series) /
                                      number(initial, "vorticity_sample"),
        })
        result[resolution] = metrics
    # The same initial Fourier coefficients must have the same padded samples.
    initial_rows = [r for r in rows if int(r["snapshot"]) == 0]
    for key in ("energy", "enstrophy", "h_half", "l3_sample", "vorticity_sample"):
        if relative(number(initial_rows[0], key), number(initial_rows[1], key)) > 1e-10:
            raise ValueError("coarse/fine initial field mismatch")
    return result


def assess(trace_rows, evidence):
    selected = {}
    for resolution, stage in (("coarse", "selected"), ("fine", "fine-validation")):
        rows = [r for r in trace_rows if r["stage"] == stage and r["resolution"] == resolution]
        if len(rows) != 1:
            raise ValueError("missing or duplicate selected trajectory")
        selected[resolution] = rows[0]
        if rows[0]["search_track"] != "amplification":
            raise ValueError("unexpected optimization track")
        steps = number(rows[0], "steps")
        if steps < 1 or steps != int(steps):
            raise ValueError("invalid accepted-step count")
    cross = {key: relative(evidence["coarse"][key], evidence["fine"][key])
             for key in COMPARISONS}
    conservative = {key: min(evidence[r][key] for r in selected)
                    for key in evidence["coarse"]}
    peak_cutoff = max(number(r, "peak_cutoff_fraction") for r in selected.values())
    numerical_failures = []
    if any(r["valid"] != "true" for r in selected.values()):
        numerical_failures.append("solver safety/invariant gate")
    if peak_cutoff > LIMITS["cutoff_margin"]:
        numerical_failures.append("cutoff margin")
    numerical_failures += [f"cross-resolution {key}" for key, value in cross.items()
                           if value > LIMITS["cross_" + key]]
    physics_failures = []
    for key, limit in (("h_final", "minimum_h_final"), ("h_late", "minimum_h_late"),
                       ("scale_final", "minimum_scale_final")):
        if conservative[key] < LIMITS[limit]:
            physics_failures.append(key)
    if conservative["stretch_late_min"] <= LIMITS["minimum_stretch_late"]:
        physics_failures.append("stretching does not dominate throughout sampled late window")
    if conservative["h_late_min_step"] < 1.0:
        physics_failures.append("late critical-norm reversal")
    return {
        "evidence": evidence, "conservative": conservative, "cross_relative": cross,
        "steps": {grid: int(number(row, "steps")) for grid, row in selected.items()},
        "peak_cutoff": peak_cutoff, "numerical_failures": numerical_failures,
        "physics_failures": physics_failures,
        "profile_route_eligible": all(r["refinement_eligible"] == "true" for r in selected.values()),
        "native_grid_vorticity_disagreement": relative(
            number(selected["coarse"], "peak_vorticity_ratio"),
            number(selected["fine"], "peak_vorticity_ratio")),
    }


def add_late_growth_assessment(result, rows, trace_rows, horizon, start_fraction, samples, temperature):
    """Validate the independent rate clock, not only the aggregate optimizer score.

    The normalized soft minimum is optimistic: min(gamma) <= softmin(gamma).
    Require the actual sampled minimum positive, and compare every common rate.
    This is still a finite sampling gate, not a continuous-time certificate.
    """
    if len(rows) != 2 * samples:
        raise ValueError("incomplete or duplicate late-rate evidence")
    by_grid = {}
    for grid, stage in (("coarse", "selected"), ("fine", "fine-validation")):
        series = [row for row in rows if row["resolution"] == grid]
        if len(series) != samples:
            raise ValueError("unbalanced late-rate evidence")
        for i, row in enumerate(series):
            expected_time = horizon * (start_fraction + (1 - start_fraction) * i / (samples - 1))
            if number(row, "snapshot") != i or abs(number(row, "time") - expected_time) > 1e-12:
                raise ValueError("late-rate clock mismatch")
        rates = [number(row, "critical_log_rate") for row in series]
        minimum = min(rates)
        exponentials = [math.exp(-(rate - minimum) / temperature) for rate in rates]
        normalizer = sum(exponentials)
        weights = [e / normalizer for e in exponentials]
        soft_minimum = minimum - temperature * math.log(normalizer / samples)
        for row, weight in zip(series, weights):
            if abs(number(row, "soft_minimum_weight") - weight) > 1e-12:
                raise ValueError("late-rate adjoint weight mismatch")
        selected = [row for row in trace_rows if row["resolution"] == grid and row["stage"] == stage]
        if len(selected) != 1 or selected[0]["growth_objective"] != "late-rate":
            raise ValueError("late-rate objective/trace mismatch")
        row = selected[0]
        for key, expected in (("horizon", horizon), ("late_window_start", start_fraction),
                              ("late_rate_samples", samples), ("late_rate_temperature", temperature),
                              ("late_rate_minimum", minimum), ("late_rate_maximum", max(rates)),
                              ("late_rate_soft_minimum", soft_minimum)):
            if abs(number(row, key) - expected) > 1e-11 * max(1, abs(expected)):
                raise ValueError("late-rate aggregate/trace mismatch")
        by_grid[grid] = {"minimum": minimum, "maximum": max(rates), "soft_minimum": soft_minimum,
                         "rates": rates}
    disagreement = horizon * max(abs(a - b) for a, b in
                                zip(by_grid["coarse"]["rates"], by_grid["fine"]["rates"]))
    conservative_minimum = min(grid["minimum"] for grid in by_grid.values())
    result["late_growth"] = {"horizon": horizon, "start_fraction": start_fraction, "samples": samples,
                             "temperature": temperature, "evidence": by_grid,
                             "conservative_minimum": conservative_minimum,
                             "cross_dimensionless": disagreement}
    if disagreement > LIMITS["cross_late_rate_dimensionless"]:
        result["numerical_failures"].append("cross-resolution late critical rate")
    if conservative_minimum <= 0:
        result["physics_failures"].append("non-positive sampled late critical rate")


def pareto_vector(result):
    m = result["conservative"]
    vector = (m["h_final"], m["h_late"], m["scale_final"], m["stretch_late_min"],
              -result["peak_cutoff"], -max(result["cross_relative"].values()))
    if "late_growth" in result:
        late = result["late_growth"]
        vector += (late["horizon"] * late["conservative_minimum"],)
    return vector


def dominates(left, right):
    # Do not decide a ranking using differences near the last reported digits.
    tolerances = [1e-4 * max(1.0, abs(a), abs(b)) for a, b in zip(left, right)]
    # Exact weak ordering prevents cycles that epsilon-relaxed weak dominance
    # can introduce. Tolerance only decides whether an improvement is material.
    return (all(a >= b for a, b in zip(left, right)) and
            any(a > b + t for a, b, t in zip(left, right, tolerances)))


def select_finalists(runs, budget):
    usable = [r for r in runs if r.get("result") and not r["result"]["numerical_failures"]]
    frontier = [r for r in usable if not any(
        dominates(pareto_vector(other["result"]), pareto_vector(r["result"]))
        for other in usable if other is not r)]
    # Round-robin objectives plus family diversity reserve room for exploration.
    chosen = []
    axes = (6, 0, 1, 2, 3, 4, 5) if usable and "late_growth" in usable[0]["result"] else (0, 1, 2, 3, 4, 5)
    for axis in axes:
        remaining = [r for r in frontier if r not in chosen]
        if not remaining or len(chosen) >= budget:
            break
        represented = {r["family"] for r in chosen}
        diverse = [r for r in remaining if r["family"] not in represented]
        pool = diverse or remaining
        chosen.append(max(pool, key=lambda r: (pareto_vector(r["result"])[axis], r["id"])))
    return frontier, chosen


def perturbation_failures(reference, repeat, label):
    failures = []
    for resolution in ("coarse", "fine"):
        if "late_growth" in reference:
            if "late_growth" not in repeat:
                failures.append(f"{label}: missing late-rate evidence")
            else:
                late, other = reference["late_growth"], repeat["late_growth"]
                difference = late["horizon"] * abs(late["evidence"][resolution]["minimum"] -
                                                   other["evidence"][resolution]["minimum"])
                if difference > LIMITS["perturbation_late_rate_dimensionless"]:
                    failures.append(f"{label}: {resolution} late critical rate")
        if label == "half-dt" and repeat["steps"][resolution] < 1.9 * reference["steps"][resolution]:
            failures.append(f"half-dt: {resolution} trajectory did not materially refine its steps")
        for key in COMPARISONS + ("h_late",):
            limit = LIMITS["perturbation_relative"]
            if key == "vorticity_final" and label == "sampling":
                limit = LIMITS["sampling_vorticity_relative"]
            if key == "stretch_late_min":
                limit = LIMITS["stretch_perturbation_relative"]
            if relative(reference["evidence"][resolution][key],
                        repeat["evidence"][resolution][key]) > limit:
                failures.append(f"{label}: {resolution} {key}")
    return failures


def gain_error_gate(results):
    values = [r["evidence"][grid]["h_final"] for r in results for grid in ("coarse", "fine")]
    error = max(values) - min(values)
    return {"minimum_gain": min(values) - 1, "observed_spread": error,
            "passed": min(values) - 1 > LIMITS["gain_to_error_factor"] * error}


def run(args):
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)
    optimizer = str(Path(args.optimizer).resolve(strict=True))
    oracle = str(Path(args.oracle).resolve(strict=True)) if args.oracle else None
    manifest = {"schema": 1, "scope": "finite-amplification screening; no PDE proof",
                "config": vars(args), "thresholds": LIMITS, "runs": [], "finalists": [],
                "optimizer_sha256": sha256(optimizer),
                "oracle_sha256": sha256(oracle) if oracle else None,
                "source_sha256": sha256(__file__)}
    research_root = Path(__file__).resolve().parents[1]
    manifest["implementation_sha256"] = {
        str(path.relative_to(research_root)): sha256(path)
        for folder in ("include", "src", "scripts")
        for path in sorted((research_root / folder).rglob("*"))
        if path.is_file() and path.suffix in (".hpp", ".cpp", ".py")}

    def save():
        temporary = output / "manifest.json.tmp"
        temporary.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
        temporary.replace(output / "manifest.json")

    def evaluate(seed, stage, horizon, dt, samples, sample_grid, iterations, source=None,
                 fine_dt=None):
        name = seed["id"] + "-" + stage
        directory = output / name
        directory.mkdir()
        command = [optimizer, "--grid", str(args.grid), "--fine-grid", str(args.fine_grid),
                   "--seed-bandwidth", str(args.bandwidth), "--energy", str(args.energy),
                   "--viscosity", str(args.viscosity), "--final-time", str(horizon),
                   "--dt", str(dt), "--fine-max-dt", str(fine_dt or dt), "--diagnostic-every", "20",
                   "--search-track", "amplification", "--initial-family", seed["family"],
                   "--starts", "1", "--start-offset", str(seed["start"] if stage == "discovery" else 0),
                   "--iterations", str(iterations), "--profile-path-samples", str(samples),
                   "--evidence-grid", str(sample_grid), "--output", str(directory / "trace.csv"),
                   "--evidence-output", str(directory / "evidence.csv"),
                   "--state-output", str(directory / "state.csv")]
        late_samples = args.late_rate_samples if stage != "sampling" else 2 * args.late_rate_samples - 1
        if args.growth_objective == "late-rate":
            command += ["--growth-objective", "late-rate", "--late-window-start", str(args.late_window_start),
                        "--late-rate-samples", str(late_samples),
                        "--late-rate-temperature", str(args.late_rate_temperature),
                        "--late-rate-output", str(directory / "late-rates.csv")]
        if source:
            command += ["--state-input", str(Path(source).resolve())]
        entry = {"id": name, "seed": seed["id"], "family": seed["family"], "stage": stage,
                 "command": command, "input_sha256": sha256(source) if source else None}
        manifest["runs"].append(entry)
        save()
        print(f"Running {name}", flush=True)
        try:
            with (directory / "run.log").open("w") as log:
                completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                           timeout=args.timeout, check=False)
            entry["returncode"] = completed.returncode
            if completed.returncode != 0:
                raise ValueError("optimizer rejected trajectory; see run.log")
            evidence = summarize_evidence(read_rows(directory / "evidence.csv"),
                                          horizon, samples, sample_grid, args.viscosity, args.energy)
            entry["result"] = assess(read_rows(directory / "trace.csv"), evidence)
            if args.growth_objective == "late-rate":
                add_late_growth_assessment(entry["result"], read_rows(directory / "late-rates.csv"),
                                           read_rows(directory / "trace.csv"), horizon,
                                           args.late_window_start, late_samples, args.late_rate_temperature)
                entry["late_rates_sha256"] = sha256(directory / "late-rates.csv")
            entry["state"] = str(directory.relative_to(output) / "state.csv")
            entry["state_sha256"] = sha256(directory / "state.csv")
            entry["evidence_sha256"] = sha256(directory / "evidence.csv")
            entry["trace_sha256"] = sha256(directory / "trace.csv")
            r = entry["result"]
            print(f"  H={r['conservative']['h_final']:.7f} late-H={r['conservative']['h_late']:.7f} "
                  f"cutoff={r['peak_cutoff']:.5g}; numerical={len(r['numerical_failures'])} "
                  f"physics={len(r['physics_failures'])} failures", flush=True)
        except (ValueError, KeyError, OSError, subprocess.TimeoutExpired) as error:
            entry["error"] = str(error)
            print(f"  rejected: {error}", flush=True)
        save()
        return entry

    seeds = [{"id": f"{family}-{start}", "family": family, "start": start}
             for family in args.families.split(",") for start in range(args.starts)]
    if args.checkpoint:
        seeds.append({"id": "checkpoint", "family": "wave-packets", "start": 0,
                      "source": args.checkpoint})
    discoveries = [evaluate(seed, "discovery", args.discovery_time, args.dt, 8,
                            args.fine_grid, args.iterations, seed.get("source")) for seed in seeds]
    frontier, selected = select_finalists(discoveries, args.max_finalists)
    manifest["pareto_frontier"] = [r["id"] for r in frontier]
    manifest["selected_for_holdout"] = [r["id"] for r in selected]
    save()
    for discovery in selected:
        seed = next(s for s in seeds if s["id"] == discovery["seed"])
        source = output / discovery["state"]
        verdict = {"seed": seed["id"], "state": discovery["state"],
                   "state_sha256": discovery["state_sha256"], "failures": [],
                   "status": "screening-only"}
        manifest["finalists"].append(verdict)
        holdout = evaluate(seed, "holdout", args.holdout_time, args.dt, 8, args.fine_grid, 0, source)
        # Cap the fine step at half its observed mean; merely halving an unused
        # maximum can leave an adaptive trajectory completely unchanged.
        fine_dt = min(args.dt / 2, args.holdout_time /
                      (2 * holdout["result"]["steps"]["fine"])) if "result" in holdout else args.dt / 2
        coarse_dt = min(args.dt / 2, args.holdout_time /
                        (2 * holdout["result"]["steps"]["coarse"])) if "result" in holdout else args.dt / 2
        checks = [holdout,
                  evaluate(seed, "half-dt", args.holdout_time, coarse_dt, 8, args.fine_grid, 0, source,
                           fine_dt=fine_dt),
                  evaluate(seed, "sampling", args.holdout_time, args.dt, 16, 2 * args.fine_grid, 0, source)]
        for check in checks:
            if "error" in check:
                verdict["failures"].append(check["id"] + ": " + check["error"])
            else:
                verdict["failures"] += [check["id"] + ": " + f for f in
                                        check["result"]["numerical_failures"] + check["result"]["physics_failures"]]
        if all("result" in c for c in checks):
            results = [c["result"] for c in checks]
            for other, label in zip(results[1:], ("half-dt", "sampling")):
                verdict["failures"] += perturbation_failures(results[0], other, label)
            verdict["gain_error"] = gain_error_gate(results)
            if not verdict["gain_error"]["passed"]:
                verdict["failures"].append("critical gain does not exceed observed numerical spread by 3x")
        if not verdict["failures"]:
            if oracle:
                command = [oracle, "--grid", str(args.fine_grid), "--viscosity", str(args.viscosity),
                           "--energy", str(args.energy), "--final-time", str(args.holdout_time),
                           "--dt", str(args.dt), "--state-input", str(source),
                           "--output", str(output / (seed["id"] + "-fftw.csv"))]
                verdict["oracle_command"] = command
                try:
                    with (output / (seed["id"] + "-fftw.log")).open("w") as log:
                        status = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                                timeout=args.timeout, check=False).returncode
                    verdict["oracle_returncode"] = status
                    if status == 0:
                        verdict["status"] = "validated-finite-amplification-shortlist"
                    else:
                        verdict["failures"].append("independent FFTW trajectory gate")
                except (OSError, subprocess.TimeoutExpired) as error:
                    verdict["failures"].append("independent FFTW: " + str(error))
            else:
                verdict["failures"].append("independent FFTW trajectory not supplied")
        save()
    print(f"Saved {len(manifest['runs'])} runs; "
          f"{sum(f['status'] == 'validated-finite-amplification-shortlist' for f in manifest['finalists'])} "
          "passed all finite-amplification gates. No PDE singularity is established.", flush=True)
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--optimizer", required=True)
    parser.add_argument("--oracle")
    parser.add_argument("--output-dir", required=True, help="new directory; refuses to overwrite a run")
    parser.add_argument("--families", default=",".join(FAMILIES))
    parser.add_argument("--checkpoint", help="optional wave-packet checkpoint with matching bandwidth and energy")
    parser.add_argument("--starts", type=int, default=2)
    parser.add_argument("--iterations", type=int, default=1)
    parser.add_argument("--max-finalists", type=int, default=3)
    parser.add_argument("--grid", type=int, default=16)
    parser.add_argument("--fine-grid", type=int, default=32)
    parser.add_argument("--bandwidth", type=int, default=3)
    parser.add_argument("--energy", type=float, default=10.0)
    parser.add_argument("--viscosity", type=float, default=0.02)
    parser.add_argument("--dt", type=float, default=0.0005)
    parser.add_argument("--discovery-time", type=float, default=0.04)
    parser.add_argument("--holdout-time", type=float, default=0.06)
    parser.add_argument("--timeout", type=float, default=1800.0)
    parser.add_argument("--growth-objective", choices=("endpoint", "late-rate"), default="endpoint")
    parser.add_argument("--late-window-start", type=float, default=0.5)
    parser.add_argument("--late-rate-samples", type=int, default=5)
    parser.add_argument("--late-rate-temperature", type=float, default=0.1)
    args = parser.parse_args()
    if (any(f not in FAMILIES for f in args.families.split(",")) or
            len(set(args.families.split(","))) != len(args.families.split(",")) or
            not 1 <= args.starts <= 64 or not 0 <= args.iterations <= 20 or
            not 0 <= args.max_finalists <= 6 or args.bandwidth < 1 or
            not math.isfinite(args.late_window_start) or not 0 < args.late_window_start < 1 or
            not 2 <= args.late_rate_samples <= 32 or
            not math.isfinite(args.late_rate_temperature) or args.late_rate_temperature <= 0 or
            args.grid < 8 or args.fine_grid <= args.grid or
            args.grid & (args.grid - 1) or args.fine_grid & (args.fine_grid - 1) or
            any(not math.isfinite(v) or v <= 0 for v in
                (args.energy, args.viscosity, args.dt, args.discovery_time, args.holdout_time, args.timeout)) or
            args.holdout_time <= args.discovery_time):
        parser.error("invalid search configuration or holdout must be longer than discovery")
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
