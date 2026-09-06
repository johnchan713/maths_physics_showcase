#!/usr/bin/env python3
"""Reproduce the frozen finalist's longer 32/64 evolution and FFTW check."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--optimizer", required=True)
parser.add_argument("--oracle", required=True)
parser.add_argument("--output-dir", required=True)
args = parser.parse_args()
output = Path(args.output_dir).resolve()
output.mkdir(parents=True, exist_ok=True)
if (output / "commands.json").exists():
    raise SystemExit("Refusing to replace an existing extension run")
research = Path(__file__).resolve().parents[2]
state = research / "candidates/wave_k3_amplification_t004_screening.csv"
optimizer, oracle = str(Path(args.optimizer).resolve()), str(Path(args.oracle).resolve())
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
commands = {
    "state_sha256": sha(state), "optimizer_sha256": sha(optimizer), "oracle_sha256": sha(oracle),
    "paired": [optimizer, "--search-track", "amplification", "--state-input", str(state),
               "--initial-family", "wave-packets", "--seed-bandwidth", "3",
               "--grid", "32", "--fine-grid", "64", "--energy", "10", "--viscosity", "0.02",
               "--dt", "0.000125", "--fine-max-dt", "0.000125", "--final-time", "0.08",
               "--iterations", "0", "--diagnostic-every", "40", "--profile-path-samples", "8",
               "--evidence-grid", "64", "--evidence-output", str(output / "evidence.csv"),
               "--output", str(output / "trace.csv"), "--state-output", str(output / "state.csv")],
    "fftw": [oracle, "--grid", "64", "--energy", "10", "--viscosity", "0.02",
             "--dt", "0.000125", "--final-time", "0.08", "--state-input", str(state),
             "--diagnostic-every", "40", "--output", str(output / "fftw.csv")],
}
for name in ("paired", "fftw"):
    (output / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
    print("Running extension " + name, flush=True)
    with (output / (name + ".log")).open("w") as log:
        result = subprocess.run(commands[name], stdout=log, stderr=subprocess.STDOUT, check=False)
    commands[name + "_returncode"] = result.returncode
    (output / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
    if result.returncode:
        raise SystemExit(result.returncode)
print("Longer 32/64 paired evolution and independent FFTW replay completed.")
