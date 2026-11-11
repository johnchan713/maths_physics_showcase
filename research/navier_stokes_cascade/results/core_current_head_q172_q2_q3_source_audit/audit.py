#!/usr/bin/env python3
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import q172

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENTS = tuple(HERE.parent / name / "evidence.json" for name in q172.PROVIDERS)

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def encode(value):
    if isinstance(value, dict): return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [encode(v) for v in value]
    if isinstance(value, Fraction): return str(value)
    return value

def generate():
    seen = set(); counts = {"nodes": 0, "edges": 0, "source_bindings": 0}
    def visit(path):
        path = path.resolve()
        if path in seen: return
        seen.add(path); counts["nodes"] += 1
        data = json.loads(path.read_text())
        for name, digest in data.get("source_hashes", {}).items():
            assert sha(PROJECT / name) == digest, name
            counts["source_bindings"] += 1
        for name, digest in data.get("dependency_evidence_sha256", {}).items():
            child = PROJECT / name
            assert sha(child) == digest, name
            counts["edges"] += 1; visit(child)
    for parent in PARENTS: visit(parent)
    own = [HERE / name for name in ("README.md", "q172.py", "audit.py", "test_audit.py")]
    readmes = [path.parent / "README.md" for path in PARENTS]
    return {
        "status": q172.record()["status"],
        "parent": "e4bfe868abf9d8c1276cb3780a35d0282e3580ea",
        "closure": counts,
        "dependency_evidence_sha256": {str(p.relative_to(PROJECT)): sha(p) for p in PARENTS},
        "source_hashes": {str(p.relative_to(PROJECT)): sha(p) for p in own + readmes},
        "review": encode(q172.review()),
        "symbolic": encode(q172.record()),
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--verify-record", type=Path)
    args = parser.parse_args(); record = generate()
    if args.verify_record: assert record == json.loads(args.verify_record.read_text())
    with args.output.open("x") as output:
        json.dump(record, output, indent=2, sort_keys=True); output.write("\n")
    print(record["status"], record["closure"])
