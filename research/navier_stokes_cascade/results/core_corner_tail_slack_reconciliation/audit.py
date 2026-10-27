#!/usr/bin/env python3
"""Freeze the closed-corner cutoff-tail slack correction."""
import argparse
import hashlib
import json
from pathlib import Path

from slack import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_corner_flatness_obstruction",
    "core_certified_complete_potential_diagonal_join",
    "core_complete_potential_diagonal_join",
    "canonical_local_realization",
    "finite_stage_induction",
    "core_certified_complete_potential_endpoint_join",
    "core_certified_complete_potential_force_join",
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
        "README.md", "slack.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    scope = result["scope"]
    gates = {
        "eight_exact_controls": len(result["checks"]) == 8 and all(result["checks"].values()),
        "seven_frozen_dependencies": len(frozen) == len(DEPENDENCIES),
        "all_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "correction_and_scope_preserved": (
            not scope["full_exponent_no_log_tail_supremum_proved"]
            and not scope["full_exponent_no_log_tail_supremum_required"]
            and scope["half_gain_closed_tail_transfer_proved_conditionally"]
            and not scope["actual_primitive_closed_coefficients_proved"]
            and not scope["compact_force_comparison_or_blowup_unconditional"]
        ),
    }
    return {
        "status": "closed-corner-cutoff-tail-slack-reconciled",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": result,
        "gates": gates,
        "dependency_evidence_sha256": frozen,
        "inherited_source_checks": inherited,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "prove the actual fixed-stage residual coefficient bound uniformly through the closed angular poles on each compact X interval",
            "prove the actual increment representative coefficient and logarithmic-power bound on the same closed chart before applying the now-closed cutoff-tail transfer",
            "then replay corner localization, compact forcing and comparison on the same field without strengthening the cutoff-tail quantifier",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Corner-tail slack correction, provenance or scope changed")


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
    print("corner-tail slack reconciliation:", len(record["review"]["checks"]),
          "exact controls;", len(record["inherited_source_checks"]),
          "dependency source checks")


if __name__ == "__main__":
    main()
