#!/usr/bin/env python3
"""Freeze the same-field protected-core reconciliation."""

import argparse
import hashlib
import json
from pathlib import Path

from core_path import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_growth_activation",
    "core_all_order_foundation",
    "core_leading_background_borel_join",
    "core_complete_potential_diagonal_join",
    "core_certified_complete_potential_diagonal_join",
    "core_scale_compatible_corner_force_replay",
    "core_certified_complete_potential_force_join",
    "core_certified_complete_potential_comparison_join",
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
        "README.md", "core_path.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    conclusions, scope = result["analytic_conclusions"], result["scope"]
    gates = {
        "six_exact_controls": len(result["checks"]) == 6 and all(result["checks"].values()),
        "eight_frozen_dependencies": len(frozen) == len(DEPENDENCIES),
        "all_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "same_field_lower_bound_joined": (
            conclusions["leading_cartesian_swirl_coefficient_is_E0_at_z_zero"]
            and conclusions["rule_defined_positive_order_background_is_o_one_after_qA_normalization"]
            and conclusions["all_annular_diagonal_increments_vanish_exactly_on_the_path"]
            and conclusions["same_complete_diagonal_has_positive_core_lower_bound"]
        ),
        "scientific_scope_preserved": (
            not scope["stress_growth_used_as_velocity_lower_bound"]
            and not scope["annular_vanishing_incorrectly_applied_to_background_coefficients"]
            and not scope["unconditional_provider_chain_or_blowup_claimed"]
        ),
    }
    return {
        "status": "same-complete-diagonal-positive-core-bound-reconciled",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": result,
        "gates": gates,
        "dependency_evidence_sha256": frozen,
        "inherited_source_checks": inherited,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "audit the consolidated field-force-comparison dependency graph for one common selected profile, scale sequences, gauges and localization",
            "enumerate and discharge every remaining conditional continuum premise before an admissible-force or blowup claim",
            "retain the fixed-time then radius then protected-core order of limits in the final comparison replay",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Protected-core reconciliation, provenance, or scope changed")


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
    print("same-field protected core:", len(record["review"]["checks"]),
          "exact controls;", len(record["inherited_source_checks"]),
          "dependency source checks")


if __name__ == "__main__":
    main()
