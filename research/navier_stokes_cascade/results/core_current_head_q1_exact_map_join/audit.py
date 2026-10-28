#!/usr/bin/env python3
"""Freeze the current-head Q1 exact defining-map join."""

import argparse
import hashlib
import json
from pathlib import Path

from exact_map import DEPENDENCIES, provider_contracts, review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    result = review()
    dependency_paths = [HERE.parent / name / "evidence.json" for name in DEPENDENCIES]
    own = [HERE / name for name in (
        "README.md", "exact_map.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    scope = result["scope"]
    return {
        "status": "current-head-Q1-exact-core-roots-and-leading-pressure-bound-on-one-profile",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": result,
        "gates": {
            "seven_exact_graph_and_interface_controls": len(result["checks"]) == 7 and all(result["checks"].values()),
            "four_Q1_contract_groups": len(result["contracts"]) == len(provider_contracts()) == 4,
            "fourteen_direct_provider_records": len(dependency_paths) == 14 and all(path.exists() for path in dependency_paths),
            "scientific_scope_preserved": scope["Q1_provider_chain_bound_on_one_current_head_object"] and scope["Q1_interface_uniqueness_and_pressure_compatibility_rechecked"] and not scope["every_provider_theorem_independently_reproved_here"] and not scope["Q2_regional_continuum_infima_discharged"] and not scope["Q3_full_frame_primitive_domain_discharged"] and not scope["hashes_or_tests_treated_as_continuum_proof"] and not scope["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"],
        },
        "dependency_evidence_sha256": {
            str(path.relative_to(PROJECT)): sha(path) for path in dependency_paths
        },
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "audit Q2 on this exact Q1 profile by proving the positive core, determinant, orientation, cone and factored-edge continuum infima on the actual finite regional cover",
            "then audit Q3 original full-frame primitive bounds, complete path and phase jets and the accepted physical domain on the same common tail",
            "retain every later pressure component; the leading pressure identity does not collapse the full corrected pressure",
            "only after Q2 and Q3 may the already replayed actual-pulse, residual, force and comparison implications be reconsidered for unconditional promotion",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Q1 graph, interface controls, provenance or scope changed")


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
    print("current-head Q1 exact-map join:", graph["nodes"], "nodes;", graph["edges"], "edges;", graph["source_hash_checks"], "source checks")


if __name__ == "__main__":
    main()
