import argparse, hashlib, json
from pathlib import Path
from controls import review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_initial_a0_physical_tail_closure",
    "first_sourced_foundation",
    "core_physical_tail_join",
    "regular_source_foundation",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    dependency_hashes, inherited = {}, {}
    for name in DEPENDENCIES:
        path = HERE.parent/name/"evidence.json"
        record = json.loads(path.read_text())
        dependency_hashes[str(path.relative_to(PROJECT))] = sha(path)
        for rel, digest in record["source_hashes"].items():
            inherited[name+":"+rel] = sha(PROJECT/rel) == digest
    result = json.loads(json.dumps(review(), default=str))
    paths = [HERE/name for name in (
        "README.md", "stage_sync.py", "controls.py", "audit.py",
        "test_audit.py", "requirements.txt", ".gitignore")]
    return dict(
        status="conditional-first-real-physical-stage-synchronization-obstruction",
        review=result,
        gates=dict(exact_controls=all(result["checks"].values()),
                   four_frozen_records=len(dependency_hashes) == 4,
                   historical_sources=all(inherited.values())),
        dependency_evidence_sha256=dependency_hashes,
        inherited_source_checks=inherited,
        source_hashes={str(path.relative_to(PROJECT)): sha(path) for path in paths},
        flags=dict(
            stage_mixing_excluded=True,
            initial_compact_tail_closure_preserved=True,
            post_sourced_regular_gain_preserved=True,
            post_sourced_changed_compact_DAG_closed=False,
            complete_physical_residual_decay_proved=False,
            blowup_certified=False,
        ),
    )


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Physical-stage synchronization record, provenance or scope changed")


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
    print("initial physical stage synchronization:",
          len(record["review"]["checks"]), "exact controls;",
          len(record["inherited_source_checks"]), "historical source checks")


if __name__ == "__main__":
    main()
