#!/usr/bin/env python3
"""Keep both the passing smoke and the initially discovered underresolved control."""
import argparse
import json
from pathlib import Path
import subprocess

from probe import ROOT, RESEARCH, save, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=ROOT / "probe_controls.json")
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Refusing to replace probe-control evidence")
    initial = RESEARCH / "candidates/wave_k3_late_rate_t004_screening.csv"
    binary = args.binary.resolve(strict=True)
    report = {"status": "passed", "binary_sha256": sha(binary),
              "source_sha256": sha(Path(__file__)), "initial_state_sha256": sha(initial), "controls": []}
    for name, grids, expected in (("resolved-short-control", (32, 64), 0),
                                  ("underresolved-short-control", (16, 32), 1),
                                  ("invalid-grid-pair", (16, 64), 1)):
        command = [str(binary), str(initial)] + [str(n) for n in grids]
        run = subprocess.run(command, capture_output=True, text=True, timeout=30)
        if run.returncode != expected:
            raise ValueError("Unexpected control result: " + name + run.stderr)
        result = json.loads(run.stdout) if run.stdout else None
        if expected == 0 and result["relative_state_l2_gap"] > 1e-10:
            raise ValueError("Passing control disagrees")
        if name == "underresolved-short-control" and (result["relative_state_l2_gap"] <= 1e-10 or
                                                       "state agreement failed" not in run.stderr):
            raise ValueError("Underresolved control failed for the wrong reason")
        if name == "invalid-grid-pair" and "nested pair" not in run.stderr:
            raise ValueError("Invalid-grid control failed for the wrong reason")
        report["controls"].append({"name": name, "command": command, "expected_returncode": expected,
                                   "returncode": run.returncode, "measurements": result, "stderr": run.stderr})
    save(args.output, report)
    print("Three probe controls passed, including expected rejection of the 16/32 case")


if __name__ == "__main__":
    main()
