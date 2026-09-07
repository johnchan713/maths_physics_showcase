#!/usr/bin/env python3
"""Capture real local regression output and hashes; does not simulate passes."""
import argparse
import json
from pathlib import Path
import subprocess

from compare import ROOT, RESEARCH, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", required=True)
    parser.add_argument("--legacy-optimizer", required=True)
    parser.add_argument("--optimizer-name", default="state_optimize")
    parser.add_argument("--output-name", default="validation")
    args = parser.parse_args()
    build = Path(args.build_dir).resolve(strict=True)
    if Path(args.output_name).name != args.output_name:
        raise ValueError("Output name must be a single new directory name")
    target = ROOT / args.output_name
    target.mkdir(exist_ok=False)
    suites = ["test_late_growth", "test_late_growth_ubsan", "test_state_optimizer", "test_galerkin",
              "test_pseudospectral", "test_candidate_score", "test_parameter_optimizer",
              "test_candidate_evidence", "test_fftw", "test_continuation"]
    commands = [(name, [str(build / name)]) for name in suites]
    commands += [
        ("robust_search", ["python3", str(RESEARCH / "tests/test_robust_search.py")]),
        ("late_growth_cli", ["python3", str(RESEARCH / "tests/test_late_growth_cli.py"),
                             "--optimizer", str(build / args.optimizer_name),
                             "--legacy-optimizer", str(Path(args.legacy_optimizer).resolve(strict=True))])]
    record = {"scope": "local regression results; clean CMake and ASan+UBSan are separate CI checks",
              "strict_build_flags": ["-std=c++11", "-O2", "-Wall", "-Wextra", "-Wpedantic", "-Werror"],
              "ubsan_flags": ["-O1", "-g", "-fsanitize=undefined", "-fno-omit-frame-pointer"],
              "runs": []}
    for name, command in commands:
        print("Validating", name, flush=True)
        with (target / (name + ".log")).open("w") as log:
            completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=180)
        record["runs"].append({"name": name, "command": command, "returncode": completed.returncode,
                               "executable_sha256": sha(command[0]) if command[0] != "python3" else None,
                               "log_sha256": sha(target / (name + ".log"))})
        (target / "manifest.json").write_text(json.dumps(record, indent=2) + "\n")
        if completed.returncode != 0:
            raise RuntimeError("Regression failed: " + name)


if __name__ == "__main__":
    main()
