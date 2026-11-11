#!/usr/bin/env python3
"""Freeze Q169's same-object pressure and comparison audit."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import q169


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENTS = tuple(HERE.parent / name / "evidence.json" for name in q169.PROVIDERS)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    if isinstance(value, Fraction):
        return str(value)
    return value


def generate():
    seen = set()
    counts = {"nodes": 0, "edges": 0, "source_bindings": 0}

    def visit(path):
        path = path.resolve()
        if path in seen:
            return
        seen.add(path)
        counts["nodes"] += 1
        data = json.loads(path.read_text())
        for name, digest in data.get("source_hashes", {}).items():
            assert sha(PROJECT / name) == digest, name
            counts["source_bindings"] += 1
        for name, digest in data.get("dependency_evidence_sha256", {}).items():
            child = PROJECT / name
            assert sha(child) == digest, name
            counts["edges"] += 1
            visit(child)

    for parent in PARENTS:
        visit(parent)
    own = [HERE / name for name in ("README.md", "q169.py", "audit.py", "test_audit.py")]
    provider_readmes = [path.parent / "README.md" for path in PARENTS]
    return {
        "status": q169.record()["status"],
        "parent": "06185eb484f89e9b767f2221e30b4dc093b083f0",
        "closure": counts,
        "dependency_evidence_sha256": {
            str(path.relative_to(PROJECT)): sha(path) for path in PARENTS
        },
        "source_hashes": {
            str(path.relative_to(PROJECT)): sha(path)
            for path in own + provider_readmes
        },
        "providers": encode(q169.provider_contracts()),
        "pressure": encode(q169.pressure_contract()),
        "commutator": encode(q169.commutator_contract()),
        "radius": encode(q169.radius_contract()),
        "remaining_gate": encode(q169.remaining_gate()),
        "frontier": encode(q169.frontier()),
        "symbolic": encode(q169.record()),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--verify-record", type=Path)
    args = parser.parse_args()
    record = generate()
    if args.verify_record:
        assert record == json.loads(args.verify_record.read_text())
    with args.output.open("x") as output:
        json.dump(record, output, indent=2, sort_keys=True)
        output.write("\n")
    print(record["status"], record["closure"])
