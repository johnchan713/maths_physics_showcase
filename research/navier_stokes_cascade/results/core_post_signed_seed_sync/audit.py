import argparse
import hashlib
import json
from pathlib import Path
from controls import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_first_sourced_physical_tail_closure",
    "core_post_signed_residual_join",
    "core_post_signed_mean_closure",
    "core_physical_tail_join",
    "finite_stage_induction",
    "core_finite_support_closure",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    frozen, inherited = {}, {}
    for name in DEPENDENCIES:
        path = HERE.parent/name/"evidence.json"
        record = json.loads(path.read_text())
        frozen[str(path.relative_to(PROJECT))] = sha(path)
        for rel, digest in record["source_hashes"].items():
            inherited[name+":"+rel] = sha(PROJECT/rel) == digest
    result = json.loads(json.dumps(review(), default=str))
    sources = [HERE/name for name in (
        "README.md", "ledger.py", "controls.py", "audit.py",
        "test_audit.py", "requirements.txt", ".gitignore")]
    sources += [
        HERE.parent/"core_first_sourced_physical_tail_closure/closure.py",
        HERE.parent/"core_post_signed_residual_join/bounds.py",
        HERE.parent/"core_post_signed_mean_closure/budgets.py",
        HERE.parent/"finite_stage_induction/induction.py",
    ]
    return dict(
        status="conditional-post-signed-induction-seed-compact-sync",
        review=result,
        gates=dict(exact_controls=all(result["checks"].values()),
                   six_frozen_records=len(frozen) == 6,
                   historical_sources=all(inherited.values())),
        dependency_evidence_sha256=frozen,
        inherited_source_checks=inherited,
        source_hashes={str(path.relative_to(PROJECT)): sha(path) for path in sources},
        flags=dict(
            post_signed_seed_compact_graph_fully_expanded=True,
            corrected_seed_compatible_with_half_line_rate_recurrence=True,
            later_stage_compact_graphs_synchronized=False,
            later_actual_source_inequalities_certified=False,
            primitive_contracts_independently_certified=False,
            complete_residual_decay_proved=False,
            blowup_certified=False,
        ),
    )


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Post-signed seed synchronization, provenance or scope changed")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify-record", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("Use a fresh output; frozen evidence is never overwritten")
    record = generate()
    validate(record)
    if args.verify_record:
        validate(json.loads(args.verify_record.read_text()))
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True)+"\n")
    print("post-signed seed compact synchronization:",
          len(record["review"]["checks"]), "exact controls;",
          len(record["inherited_source_checks"]), "historical source checks")


if __name__ == "__main__":
    main()
