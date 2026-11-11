#!/usr/bin/env python3
"""Q164 dependency replay; diagnostics do not replace continuum estimates."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import q164

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENTS = tuple(HERE.parent / name / "evidence.json" for name in q164.PROVIDERS)


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
    own = [HERE / name for name in ("README.md", "q164.py", "audit.py", "test_audit.py")]
    provider_readmes = [path.parent / "README.md" for path in PARENTS]
    return {
        "status": q164.record()["status"],
        "parent": "6d0dacc5b2ec662857cb40eafdc750c22bcc3c84",
        "closure": counts,
        "dependency_evidence_sha256": {
            str(path.relative_to(PROJECT)): sha(path) for path in PARENTS
        },
        "source_hashes": {
            str(path.relative_to(PROJECT)): sha(path)
            for path in own + provider_readmes
        },
        "operator": encode(q164.operator_contract()),
        "source": encode(q164.source_contract()),
        "contraction": encode(q164.contraction_contract()),
        "equations": encode(q164.equation_contract()),
        "pressure": encode(q164.pressure_contract()),
        "providers": encode(q164.provider_contracts()),
        "remaining_gate": encode(q164.remaining_gate()),
        "frontier": encode(q164.frontier()),
        "symbolic": encode(q164.record()),
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
