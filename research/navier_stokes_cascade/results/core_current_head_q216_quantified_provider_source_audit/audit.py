#!/usr/bin/env python3
import argparse
import hashlib
import json
from pathlib import Path
import q216

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENTS = tuple(HERE.parent / name / "evidence.json" for name in q216.PROVIDERS)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
    own = [HERE / name for name in ("README.md", "q216.py", "audit.py", "test_audit.py")]
    readmes = [parent.parent / "README.md" for parent in PARENTS]
    sources = own + readmes
    return {
        "status": q216.record()["status"],
        "parent": "c0912fa025e08f6cac42f53598efda5e53ea6a8c",
        "closure": counts,
        "dependency_evidence_sha256": {
            str(path.relative_to(PROJECT)): sha(path) for path in PARENTS
        },
        "source_hashes": {
            str(path.relative_to(PROJECT)): sha(path) for path in sources
        },
        "review": q216.review(),
        "symbolic": q216.record(),
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
