#!/usr/bin/env python3
import argparse
import hashlib
import json
from pathlib import Path
import q213

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENTS = tuple(HERE.parent / name / "evidence.json" for name in q213.PROVIDERS)
EXTRA_SOURCES = (
    HERE.parent / "core_leading_profile_certificate/README.md",
    HERE.parent / "core_leading_background_seminorm_reconstruction/README.md",
    HERE.parent / "core_same_band_gap_reconstruction/README.md",
    HERE.parent / "core_full_path_ode_reconstruction/README.md",
    HERE.parent / "core_fixed_profile_projector_reconstruction/README.md",
    HERE.parent / "core_compact_slow_inverse_reconstruction/README.md",
)


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
    own = [HERE / name for name in ("README.md", "q213.py", "audit.py", "test_audit.py")]
    readmes = [parent.parent / "README.md" for parent in PARENTS]
    sources = own + readmes + list(EXTRA_SOURCES)
    return {
        "status": q213.record()["status"],
        "parent": "9d84632658fde16893666aa59372e3986bd55fa9",
        "closure": counts,
        "dependency_evidence_sha256": {
            str(path.relative_to(PROJECT)): sha(path) for path in PARENTS
        },
        "source_hashes": {
            str(path.relative_to(PROJECT)): sha(path) for path in sources
        },
        "review": q213.review(),
        "symbolic": q213.record(),
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
