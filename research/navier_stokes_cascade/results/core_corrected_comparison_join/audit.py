#!/usr/bin/env python3
"""Freeze corrected-field pressure recovery and competitor comparison."""
import argparse
import hashlib
import json
from pathlib import Path
from controls import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_corrected_force_join",
    "comparison_foundation",
    "global_solution_comparison",
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
        "README.md", "comparison.py", "controls.py", "audit.py",
        "test_audit.py", "requirements.txt", ".gitignore",
    )]
    gates = dict(
        exact_controls=all(result["checks"].values()),
        three_frozen_dependency_records=len(frozen) == 3,
        dependency_sources_unchanged=len(inherited) == 35 and all(inherited.values()),
        conditional_scope=(
            result["scope"]["global_smooth_bounded_energy_competitor_exclusion_follows_conditionally_from_inherited_core_growth"]
            and not result["scope"]["primitive_field_force_and_positive_core_contracts_independently_proved"]
            and not result["scope"]["unconditional_admissible_force_or_blowup_proved"]
        ),
    )
    return dict(
        status="conditional-corrected-same-field-pressure-comparison-join",
        review=result, gates=gates,
        dependency_evidence_sha256=frozen,
        inherited_source_checks=inherited,
        source_hashes={str(path.relative_to(PROJECT)): sha(path) for path in own},
    )


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Corrected comparison join, provenance or scope changed")


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
    print("corrected comparison join:", len(record["review"]["checks"]),
          "exact sections;", len(record["inherited_source_checks"]),
          "dependency source checks")


if __name__ == "__main__":
    main()
