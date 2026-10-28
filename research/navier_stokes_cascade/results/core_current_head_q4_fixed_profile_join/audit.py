#!/usr/bin/env python3
"""Freeze the current-head Q4 fixed-profile compensator/projector join."""

import argparse
import hashlib
import json
from pathlib import Path

from fixed_profiles import DEPENDENCIES, review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    result = review()
    dependency_paths = [HERE.parent / name / "evidence.json" for name in DEPENDENCIES]
    own = [HERE / name for name in (
        "README.md", "fixed_profiles.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    scope = result["scope"]
    return {
        "status": "current-head-Q4-fixed-profile-compensator-and-same-choice-projectors-bound",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": result,
        "gates": {
            "ten_exact_graph_and_interface_controls": len(result["checks"]) == 10 and all(result["checks"].values()),
            "ten_direct_provider_records": len(dependency_paths) == 10 and all(path.exists() for path in dependency_paths),
            "scientific_scope_preserved": scope["Q1_Q3_field_is_a_direct_dependency"] and scope["Q4_fixed_profile_compensator_projectors_bound_on_same_later_tail"] and scope["same_choice_projector_dependent_objects_recomputed"] and not scope["normalization_only_used_to_control_arbitrary_profile_jets"] and not scope["moving_support_enclosure_identified_with_actual_support"] and not scope["historical_unspecified_projectors_identified_with_new_choice"] and not scope["uniformity_over_profiles_orders_or_infinite_cycles_claimed"] and not scope["compact_slow_weighted_inverse_contract_closed"] and not scope["hashes_or_tests_treated_as_continuum_proof"] and not scope["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"],
        },
        "dependency_evidence_sha256": {
            str(path.relative_to(PROJECT)): sha(path) for path in dependency_paths
        },
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "join the shifted compact maps, zero-Haar temporal inverse and fixed slow response to these exact same-choice coefficients on their full weighted domains and whole-path supports",
            "prove the actual generated-source weighted suprema at every finite request rather than inferring them from operator boundedness",
            "then join complete covariance and fresh axial/full-pressure means, retaining both Pi1/Pi2 tails and all reserved supports",
            "only after those joins may the full viscosity/cutoff/interaction residual, force, comparison or blow-up claims be reconsidered",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Q4 graph, profile/projector interfaces, provenance or scope changed")


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
    print("current-head Q4 fixed-profile join:", graph["nodes"], "nodes;", graph["edges"], "edges;", graph["source_hash_checks"], "source checks")


if __name__ == "__main__":
    main()
