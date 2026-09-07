#!/usr/bin/env python3
"""Score endpoint-selected fields under the frozen late objective, without optimizing.

Discovery starts/iteration budgets are matched; holdout selection policies differ.
Retain raw replays so an aggregate score cannot hide a negative sampled rate.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parent
RESEARCH = ROOT.parents[1]
SPEC = importlib.util.spec_from_file_location("robust_search", RESEARCH / "scripts/robust_search.py")
SEARCH = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SEARCH)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def selected_rows(path):
    rows = SEARCH.read_rows(path)
    return rows, {r["resolution"]: r for r in rows if r["stage"] in ("selected", "fine-validation")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--optimizer", required=True)
    args = parser.parse_args()
    optimizer = str(Path(args.optimizer).resolve(strict=True))
    endpoint = json.loads((ROOT / "endpoint/manifest.json").read_text())
    late = json.loads((ROOT / "late/manifest.json").read_text())
    config = late["config"]
    for key in ("grid", "fine_grid", "bandwidth", "energy", "viscosity", "dt", "discovery_time",
                "iterations", "starts", "families", "checkpoint"):
        if config[key] != endpoint["config"][key]:
            raise ValueError("Mismatched discovery protocol: " + key)
    if sha(optimizer) != late["optimizer_sha256"] or sha(optimizer) != endpoint["optimizer_sha256"]:
        raise ValueError("Comparison must use the same optimizer executable")
    target = ROOT / "endpoint_late_replays"
    target.mkdir(exist_ok=False)
    report = {"scope": "four matched discovery comparisons; not a PDE proof or global search",
              "optimizer_sha256": sha(optimizer), "pairs": []}
    discoveries = {r["seed"]: r for r in late["runs"] if r["stage"] == "discovery"}
    for original in (r for r in endpoint["runs"] if r["stage"] == "discovery"):
        matched = discoveries[original["seed"]]
        if "result" not in original or "result" not in matched:
            report["pairs"].append({"seed": original["seed"], "status": "missing-discovery-result"})
            continue
        folder = target / original["seed"]
        folder.mkdir()
        source = ROOT / "endpoint" / original["state"]
        command = [optimizer, "--grid", str(config["grid"]), "--fine-grid", str(config["fine_grid"]),
                   "--seed-bandwidth", str(config["bandwidth"]), "--energy", str(config["energy"]),
                   "--viscosity", str(config["viscosity"]), "--dt", str(config["dt"]),
                   "--fine-max-dt", str(config["dt"]), "--final-time", str(config["discovery_time"]),
                   "--initial-family", original["family"], "--state-input", str(source),
                   "--iterations", "0", "--diagnostic-every", "20", "--profile-path-samples", "8",
                   "--search-track", "amplification", "--growth-objective", "late-rate",
                   "--late-window-start", str(config["late_window_start"]),
                   "--late-rate-samples", str(config["late_rate_samples"]),
                   "--late-rate-temperature", str(config["late_rate_temperature"]),
                   "--output", str(folder / "trace.csv"), "--state-output", str(folder / "state.csv"),
                   "--evidence-output", str(folder / "evidence.csv"),
                   "--late-rate-output", str(folder / "late-rates.csv")]
        entry = {"seed": original["seed"], "command": command, "input_sha256": sha(source)}
        report["pairs"].append(entry)
        print("Scoring matched endpoint state:", original["seed"], flush=True)
        with (folder / "run.log").open("w") as log:
            entry["returncode"] = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                                  timeout=1800, check=False).returncode
        if entry["returncode"] != 0:
            raise RuntimeError("Endpoint scoring replay failed: " + original["seed"])
        if sha(folder / "state.csv") != entry["input_sha256"]:
            raise ValueError("Scoring replay modified frozen initial coefficients")
        # Default late times are a subset of the endpoint clock. The replay
        # therefore must reproduce every exported physical measurement exactly.
        previous_evidence = ROOT / "endpoint" / Path(original["state"]).parent / "evidence.csv"
        if (folder / "evidence.csv").read_bytes() != previous_evidence.read_bytes():
            raise ValueError("Matched scoring replay changed the endpoint trajectory")
        trace, selected = selected_rows(folder / "trace.csv")
        evidence = SEARCH.summarize_evidence(SEARCH.read_rows(folder / "evidence.csv"),
            config["discovery_time"], 8, config["fine_grid"], config["viscosity"], config["energy"])
        result = SEARCH.assess(trace, evidence)
        SEARCH.add_late_growth_assessment(result, SEARCH.read_rows(folder / "late-rates.csv"), trace,
            config["discovery_time"], config["late_window_start"], config["late_rate_samples"],
            config["late_rate_temperature"])
        _, late_selected = selected_rows(ROOT / "late" / Path(matched["state"]).parent / "trace.csv")
        entry.update(endpoint_as_late=result, late_search=matched["result"],
            endpoint_state_sha256=original["state_sha256"], late_state_sha256=matched["state_sha256"],
            matched_endpoint_late_objective=min(float(r["objective"]) for r in selected.values()),
            matched_late_objective=min(float(r["objective"]) for r in late_selected.values()),
            evidence_byte_identical=True)
        entry["conservative_minimum_rate_change"] = (matched["result"]["late_growth"]["conservative_minimum"] -
                                                       result["late_growth"]["conservative_minimum"])
        entry["files"] = {p.name: sha(p) for p in sorted(folder.iterdir()) if p.is_file()}
        (ROOT / "comparison.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    for entry in report["pairs"]:
        if "conservative_minimum_rate_change" in entry:
            print(entry["seed"], "minimum rate change:", entry["conservative_minimum_rate_change"])


if __name__ == "__main__":
    main()
