#!/usr/bin/env python3
"""Freeze the current-head Q5 weighted inverse and log-source join."""

import argparse
import hashlib
import json
from pathlib import Path

from weighted_inverse import DEPENDENCIES, review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    result = review()
    dependency_paths = [HERE.parent / name / "evidence.json" for name in DEPENDENCIES]
    own = [HERE / name for name in (
        "README.md", "weighted_inverse.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    scope = result["scope"]
    return {
        "status": "current-head-Q5-weighted-inverse-and-log-source-domain-bound",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": result,
        "gates": {
            "ten_exact_graph_and_interface_controls": len(result["checks"]) == 10 and all(result["checks"].values()),
            "ten_direct_provider_records": len(dependency_paths) == 10 and all(path.exists() for path in dependency_paths),
            "scientific_scope_preserved": scope["Q1_Q4_object_is_a_direct_dependency"] and scope["Q5_weighted_inverse_and_log_source_contract_bound_on_same_later_tail"] and scope["finite_stage_log_source_induction_instantiated"] and not scope["raw_pure_power_supremum_promoted_as_exact_current_source_domain"] and not scope["old_zero_margin_equality_erased"] and not scope["physical_characteristic_shift_removed"] and not scope["actual_A1_axial_velocity_removed"] and not scope["six_group_complete_residual_replayed"] and not scope["hashes_or_tests_treated_as_continuum_proof"] and not scope["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"],
        },
        "dependency_evidence_sha256": {
            str(path.relative_to(PROJECT)): sha(path) for path in dependency_paths
        },
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "instantiate the six-group residual replay on this exact Q1-Q5 object with separate compact-source, regular-residual and cutoff-tail logarithmic accounts",
            "retain complete repeated covariance, fresh axial velocity and whole moving-pressure means at every finite stage",
            "include unit viscosity, all potential and pressure cutoff derivatives and both nonlinear cross orientations plus the new square",
            "only after the residual closes may corner/endpoint, localization, compact force, energy, comparison, exclusion or blow-up claims be reconsidered",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Q5 graph, weighted inverse/source interfaces, provenance or scope changed")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--verify-record", type=Path)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit("Use a fresh output; frozen evidence is never overwritten")
    record = generate()
    validate(record)
    if args.verify_record:
        validate(json.loads(args.verify_record.read_text()))
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    graph = record["review"]["graph"]
    print("current-head Q5 weighted-inverse join:", graph["nodes"], "nodes;", graph["edges"], "edges;", graph["source_hash_checks"], "source checks")


if __name__ == "__main__":
    main()
