#!/usr/bin/env python3
"""Run the frozen, time- and address-space-bounded cost probe."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import time

ROOT = Path(__file__).resolve().parent
RESEARCH = ROOT.parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads((ROOT / "protocol.json").read_text())
    if (ROOT / "probe.json").exists():
        raise ValueError("Refusing to overwrite a cost-probe record")
    for name, digest in protocol["source_sha256"].items():
        if sha(RESEARCH / name) != digest:
            raise ValueError("Frozen dependency changed: " + name)
    initial = RESEARCH / protocol["initial_state"]
    if sha(initial) != protocol["initial_state_sha256"]:
        raise ValueError("Frozen initial field changed")
    settings = protocol["cost_probe"]
    binary = args.binary.resolve(strict=True)
    command = [str(binary), str(initial)] + [str(n) for n in settings["grids"]]
    report = {"status": "running", "command": command, "binary_sha256": sha(binary),
              "protocol_sha256": sha(ROOT / "protocol.json"),
              "initial_state_sha256": sha(initial), "settings": settings,
              "scope": "Two tiny RK4 steps, not a full finer-grid evolution through .14"}
    save(ROOT / "probe.json", report)

    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (settings["address_space_cap_bytes"],) * 2)

    started = time.monotonic()
    try:
        with (ROOT / "probe.stdout").open("x") as output, (ROOT / "probe.log").open("x") as log:
            completed = subprocess.run(command, stdout=output, stderr=log, preexec_fn=limits,
                                       timeout=settings["timeout_seconds"])
        report["returncode"] = completed.returncode
        report["status"] = "short-probe-passed" if completed.returncode == 0 else "probe-failed"
        if completed.returncode == 0:
            report["measurements"] = json.loads((ROOT / "probe.stdout").read_text())
    except subprocess.TimeoutExpired:
        report["status"] = "probe-timeout"
    finally:
        report["wall_seconds"] = time.monotonic() - started
        report["maximum_child_rss_kib"] = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        report["stdout_sha256"] = sha(ROOT / "probe.stdout")
        report["log_sha256"] = sha(ROOT / "probe.log")
        save(ROOT / "probe.json", report)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
