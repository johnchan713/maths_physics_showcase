#!/usr/bin/env python3
"""Frozen-leader resolution check after the failed low-grid held-out screen.

No ascent steps, no threshold relaxation. Resolving a failed coarse run is a
separate follow-up experiment, not a retroactive pass of the original pilot.
"""
import argparse
import json
from pathlib import Path
import subprocess

from compare import ROOT, SEARCH, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--optimizer", required=True)
    parser.add_argument("--oracle", required=True)
    parser.add_argument("--resume", action="store_true", help="reassess completed output and run missing checks")
    args = parser.parse_args()
    optimizer, oracle = (str(Path(p).resolve(strict=True)) for p in (args.optimizer, args.oracle))
    pilot = json.loads((ROOT / "late/manifest.json").read_text())
    leader = next(r for r in pilot["runs"] if r["id"] == "checkpoint-discovery")
    source = ROOT / "late" / leader["state"]
    if sha(source) != leader["state_sha256"]:
        raise ValueError("Leader initial coefficients changed")
    target = ROOT / "leader_32_64"
    target.mkdir(exist_ok=args.resume)
    record = {"scope": "separate frozen-field resolution follow-up; numerical screening only",
              "state_sha256": sha(source), "optimizer_sha256": sha(optimizer), "oracle_sha256": sha(oracle),
              "thresholds": SEARCH.LIMITS, "horizon": 0.06, "runs": [], "failures": []}
    if args.resume:
        record = json.loads((target / "manifest.json").read_text())
        if record["state_sha256"] != sha(source) or record["horizon"] != 0.06:
            raise ValueError("Resume changed the frozen field or horizon")
        record["resume_optimizer_sha256"] = sha(optimizer)

    def save():
        (target / "manifest.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")

    def run(stage, dt, fine_dt, samples, sampling_grid, rate_samples):
        folder = target / stage
        command = [optimizer, "--grid", "32", "--fine-grid", "64", "--seed-bandwidth", "3",
                   "--energy", "10", "--viscosity", "0.02", "--state-input", str(source),
                   "--iterations", "0", "--final-time", "0.06", "--dt", str(dt),
                   "--fine-max-dt", str(fine_dt), "--search-track", "amplification",
                   "--growth-objective", "late-rate", "--late-window-start", "0.5",
                   "--late-rate-temperature", "0.1", "--late-rate-samples", str(rate_samples),
                   "--profile-path-samples", str(samples), "--evidence-grid", str(sampling_grid),
                   "--diagnostic-every", "40", "--output", str(folder / "trace.csv"),
                   "--state-output", str(folder / "state.csv"), "--evidence-output", str(folder / "evidence.csv"),
                   "--late-rate-output", str(folder / "late-rates.csv")]
        existing = [r for r in record["runs"] if r["stage"] == stage]
        if existing:
            entry = existing[0]
            if not args.resume or entry["command"][1:] != command[1:] or entry.get("returncode") != 0:
                raise ValueError("Cannot reuse an incomplete or changed leader command")
            entry["reassessed_existing_output"] = True
        else:
            folder.mkdir()
            entry = {"stage": stage, "command": command, "executable_sha256": sha(optimizer)}
            record["runs"].append(entry)
            save()
            print("Frozen 32/64 leader:", stage, flush=True)
            with (folder / "run.log").open("w") as log:
                entry["returncode"] = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                                      timeout=1800, check=False).returncode
            save()
        if entry["returncode"] != 0:
            entry["failure"] = "solver failure; see log/trace"
            record["failures"].append(stage + ": solver failure")
            save()
            return None
        # CSV metadata describes the writing grid, so moving 16 -> 32 changes
        # its raw file hash. Every physical coefficient and constraint must
        # remain identical; do not confuse metadata with changed initial data.
        def coefficients(path):
            return [{key: value for key, value in row.items()
                     if key not in ("source_grid", "simulation_cutoff")}
                    for row in SEARCH.read_rows(path)]
        if coefficients(folder / "state.csv") != coefficients(source):
            raise ValueError("Frozen replay changed an initial Fourier coefficient or constraint")
        entry["coefficient_identity"] = True
        entry["output_initial_csv_sha256"] = sha(folder / "state.csv")
        rows = SEARCH.read_rows(folder / "trace.csv")
        evidence = SEARCH.summarize_evidence(SEARCH.read_rows(folder / "evidence.csv"),
                                             0.06, samples, sampling_grid, 0.02, 10)
        result = SEARCH.assess(rows, evidence)
        SEARCH.add_late_growth_assessment(result, SEARCH.read_rows(folder / "late-rates.csv"),
                                          rows, 0.06, 0.5, rate_samples, 0.1)
        entry["result"] = result
        entry["files"] = {p.name: sha(p) for p in sorted(folder.iterdir()) if p.is_file()}
        record["failures"] += [stage + ": " + failure for failure in
                                result["numerical_failures"] + result["physics_failures"]]
        save()
        return result

    reference = run("holdout", 0.000125, 0.000125, 8, 64, 5)
    if reference is not None:
        coarse_dt = min(0.0000625, 0.06 / (2 * reference["steps"]["coarse"]))
        fine_dt = min(0.0000625, 0.06 / (2 * reference["steps"]["fine"]))
        refined = run("half-dt", coarse_dt, fine_dt, 8, 64, 5)
        sampled = run("sampling", 0.000125, 0.000125, 16, 128, 9)
        for result, label in ((refined, "half-dt"), (sampled, "sampling")):
            if result is not None:
                record["failures"] += SEARCH.perturbation_failures(reference, result, label)
        if refined is not None and sampled is not None:
            record["gain_error"] = SEARCH.gain_error_gate([reference, refined, sampled])
            if not record["gain_error"]["passed"]:
                record["failures"].append("critical gain does not exceed observed spread by 3x")
    if not record["failures"]:
        command = [oracle, "--grid", "64", "--viscosity", "0.02", "--energy", "10",
                   "--final-time", "0.06", "--dt", "0.000125", "--state-input", str(source),
                   "--diagnostic-every", "40", "--output", str(target / "fftw.csv")]
        record["oracle_command"] = command
        save()
        print("Independent 64-grid trajectory", flush=True)
        with (target / "fftw.log").open("w") as log:
            record["oracle_returncode"] = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                                          timeout=1800, check=False).returncode
        if record["oracle_returncode"] != 0:
            record["failures"].append("independent FFTW trajectory failed")
    record["status"] = "screening-only" if record["failures"] else "validated-finite-amplification-shortlist"
    save()
    print(record["status"], record["failures"], flush=True)


if __name__ == "__main__":
    main()
