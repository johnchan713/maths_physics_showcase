#!/usr/bin/env python3
"""Freeze the corrected same-field localization and force join."""
import argparse
import hashlib
import json
from pathlib import Path
from controls import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_corrected_endpoint_join",
    "core_complete_potential_diagonal_join",
    "force_foundation",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dependency_review():
    records, sources = {}, {}
    for name in DEPENDENCIES:
        path = HERE.parent / name / "evidence.json"
        record = json.loads(path.read_text())
        records[str(path.relative_to(PROJECT))] = sha(path)
        for rel, digest in record.get("source_hashes", {}).items():
            target = PROJECT / rel
            sources[name + ":" + rel] = target.exists() and sha(target) == digest
    return records, sources


def generate():
    frozen, inherited = dependency_review()
    result = review()
    own = [HERE / name for name in (
        "README.md", "localization.py", "controls.py", "audit.py",
        "test_audit.py", "requirements.txt", ".gitignore",
    )]
    gates = dict(
        exact_controls=all(result["checks"].values()),
        three_frozen_dependency_records=len(frozen) == 3,
        dependency_sources_unchanged=len(inherited) == 39 and all(inherited.values()),
        conditional_scope=(
            result["scope"]["compact_smooth_force_and_zero_initial_data_follow_conditionally"]
            and not result["scope"]["primitive_corner_and_punctured_bounds_independently_proved"]
            and not result["scope"]["global_competitor_comparison_or_blowup_proved"]
        ),
    )
    return dict(
        status="conditional-corrected-same-field-localization-force-and-energy-join",
        review=result,
        gates=gates,
        dependency_evidence_sha256=frozen,
        inherited_source_checks=inherited,
        source_hashes={str(path.relative_to(PROJECT)): sha(path) for path in own},
    )


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Corrected force join, provenance or scope changed")


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify-record", type=Path)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit("Use a fresh output; frozen evidence is never overwritten")
    record = generate()
    validate(record)
    if args.verify_record:
        validate(json.loads(args.verify_record.read_text()))
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print("corrected same-field force join:", len(record["review"]["checks"]),
          "exact sections;", len(record["inherited_source_checks"]),
          "dependency source checks")


if __name__ == "__main__":
    main()
