#!/usr/bin/env python3
"""Verify archived numerical evidence before final publication; fail on drift."""
import json
from pathlib import Path

from compare import ROOT, RESEARCH, SEARCH, sha


def main():
    checked = []

    def check(path, expected):
        if sha(path) != expected:
            raise ValueError("Archived file checksum mismatch: " + str(path))
        checked.append(str(path.relative_to(ROOT)))

    for arm in ("endpoint", "late"):
        manifest = json.loads((ROOT / arm / "manifest.json").read_text())
        for name, expected in manifest["implementation_sha256"].items():
            source = ROOT / "pre_io_source" / name
            if not source.exists():
                source = RESEARCH / name
            if sha(source) != expected:
                raise ValueError("Discovery source is not reproducible: " + name)
        for run in manifest["runs"]:
            if "result" not in run:
                continue  # Failures remain failures, with their original logs.
            directory = ROOT / arm / Path(run["state"]).parent
            for filename, key in (("state.csv", "state_sha256"), ("evidence.csv", "evidence_sha256"),
                                  ("trace.csv", "trace_sha256"), ("late-rates.csv", "late_rates_sha256")):
                if key in run:
                    check(directory / filename, run[key])
        if len(manifest["runs"]) != 10 or len(manifest["finalists"]) != 2:
            raise ValueError("Incomplete bounded pilot")
    comparison = json.loads((ROOT / "comparison.json").read_text())
    if len(comparison["pairs"]) != 4:
        raise ValueError("Incomplete matched comparison")
    for pair in comparison["pairs"]:
        for filename, expected in pair["files"].items():
            check(ROOT / "endpoint_late_replays" / pair["seed"] / filename, expected)
    final = json.loads((ROOT / "validation_final/manifest.json").read_text())
    if len(final["runs"]) != 12 or any(r["returncode"] != 0 for r in final["runs"]):
        raise ValueError("Final local regression pass is incomplete")
    for run in final["runs"]:
        check(ROOT / "validation_final" / (run["name"] + ".log"), run["log_sha256"])
    leader = json.loads((ROOT / "leader_32_64/manifest.json").read_text())
    if len(leader["runs"]) != 3 or "status" not in leader:
        raise ValueError("Leading-field resolution follow-up is unfinished")
    for run in leader["runs"]:
        if "result" in run:
            for filename, expected in run["files"].items():
                check(ROOT / "leader_32_64" / run["stage"] / filename, expected)
    if leader["status"] == "validated-finite-amplification-shortlist":
        if leader["failures"] or leader.get("oracle_returncode") != 0:
            raise ValueError("Validated label is unsupported")
        oracle = SEARCH.read_rows(ROOT / "leader_32_64/fftw.csv")
        if not oracle:
            raise ValueError("Missing independent trajectory evidence")
    files = {str(p.relative_to(ROOT)): sha(p) for p in sorted(ROOT.rglob("*"))
             if p.is_file() and p.name != "audit.json" and "__pycache__" not in p.parts}
    report = {"scope": "checksums and completed numerical records; no continuum error bound",
              "discovery_arms": 2, "discovery_starts": 8, "matched_pairs": 4,
              "local_final_checks": 12, "verified_recorded_hashes": len(checked),
              "leader_status": leader["status"], "files": files}
    (ROOT / "audit.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print("Verified", len(checked), "recorded hashes and", len(files), "archived files")


if __name__ == "__main__":
    main()
