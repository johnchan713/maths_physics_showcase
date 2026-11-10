"""Q148 frozen provenance replay; diagnostics are not continuum proofs."""
import argparse
import hashlib
import json
from pathlib import Path
import q148

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENTS = (
    HERE.parent / "core_current_head_q147_repair_pressure_join" / "evidence.json",
    HERE.parent / "core_coefficient_support_recheck" / "evidence.json",
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
    paths = [HERE / name for name in ("README.md", "q148.py", "audit.py", "test_audit.py")]
    providers = (
        "core_current_head_q147_repair_pressure_join",
        "core_coefficient_support_recheck",
        "coefficient_support_foundation",
        "coefficient_support_audit",
        "core_current_head_q79_ambient_value_bridge",
        "core_current_head_q80_support_saturated_tuple",
        "core_current_head_q120_closed_borel_cutoff",
        "core_current_head_q122_positive_order_repair",
        "physical_pulse_transfer",
    )
    paths += [HERE.parent / name / "README.md" for name in providers]
    return {
        "status": q148.record()["status"],
        "parent": "45efb625b15fd4039b21086895a4c2fcecf1986f",
        "closure": counts,
        "dependency_evidence_sha256": {
            str(path.relative_to(PROJECT)): sha(path) for path in PARENTS
        },
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in paths},
        "symbolic": q148.record(),
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
