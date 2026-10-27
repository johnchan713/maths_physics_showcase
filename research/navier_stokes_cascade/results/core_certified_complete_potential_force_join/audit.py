#!/usr/bin/env python3
"""Freeze the certified complete-potential localization and force join."""
import argparse
import hashlib
import json
from pathlib import Path
from controls import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_certified_complete_potential_endpoint_join",
    "core_certified_complete_potential_diagonal_join",
    "core_corrected_force_join",
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
        "test_audit.py", "protocol.json", "requirements.txt", ".gitignore",
    )]
    gates = dict(
        exact_controls=all(result["checks"].values()),
        four_frozen_dependency_records=len(frozen) == 4,
        dependency_sources_unchanged=bool(inherited) and all(inherited.values()),
        conditional_scope=(
            result["scope"]["compact_smooth_force_and_zero_initial_data_follow_conditionally"]
            and not result["scope"]["primitive_corner_and_punctured_bounds_independently_proved"]
            and not result["scope"]["global_competitor_comparison_or_blowup_proved"]
        ),
    )
    return dict(
        status="certified-domain-conditional-localization-force-and-energy-join",
        review=result,
        gates=gates,
        dependency_evidence_sha256=frozen,
        inherited_source_checks=inherited,
        source_hashes={str(path.relative_to(PROJECT)): sha(path) for path in own},
        remaining_obligations=[
            "recover the fixed-gauge pressure and close expanding-ball comparison for every admissible smooth finite-energy competitor with the same zero datum and force",
            "discharge the retained primitive background, root, inverse, support, pressure and physical-transfer estimates before any unconditional force or blow-up claim",
        ],
    )


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Certified complete-potential force join, provenance or scope changed")


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
    print("certified complete-potential force join:", len(record["review"]["checks"]),
          "exact sections;", len(record["inherited_source_checks"]),
          "dependency source checks")


if __name__ == "__main__":
    main()
