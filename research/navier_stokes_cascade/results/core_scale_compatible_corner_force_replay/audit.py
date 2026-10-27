#!/usr/bin/env python3
"""Freeze the scale-compatible corner/force replay."""

import argparse
import hashlib
import json
from pathlib import Path

from corner_force import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_primitive_margin_normalization_reconciliation",
    "core_closed_chart_coefficient_join",
    "core_corner_tail_slack_reconciliation",
    "core_corner_flatness_obstruction",
    "core_certified_complete_potential_diagonal_join",
    "core_certified_complete_potential_endpoint_join",
    "core_certified_complete_potential_force_join",
    "core_same_field_pressure_reconciliation",
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
        "README.md", "corner_force.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    scope = result["scope"]
    conclusions = result["analytic_conclusions"]
    gates = {
        "six_exact_controls": len(result["checks"]) == 6 and all(result["checks"].values()),
        "eight_frozen_dependencies": len(frozen) == len(DEPENDENCIES),
        "all_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "corner_and_force_replay_joined": (
            conclusions["all_fixed_residual_jets_are_cartesian_corner_flat_relative_to_the_certified_provider_chain"]
            and conclusions["same_field_compact_force_and_preterminal_energy_follow_relative_to_punctured_endpoint_bounds"]
        ),
        "remaining_core_scope_preserved": (
            not scope["same_field_positive_core_lower_bound_independently_proved_here"]
            and not scope["unconditional_blowup_or_competitor_exclusion_proved"]
        ),
    }
    return {
        "status": "scale-compatible-corner-flatness-and-force-replayed",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": result,
        "gates": gates,
        "dependency_evidence_sha256": frozen,
        "inherited_source_checks": inherited,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "prove the positive physical swirl coefficient on the same complete-potential diagonal core path",
            "prove exact vanishing of every correction on that path or the summed one-half leading bound",
            "only then replay the fixed-time expanding-ball limit before considering competitor exclusion",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Corner/force replay, provenance, or scope changed")


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
    print("corner/force replay:", len(record["review"]["checks"]),
          "exact controls;", len(record["inherited_source_checks"]),
          "dependency source checks")


if __name__ == "__main__":
    main()
