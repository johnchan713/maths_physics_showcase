#!/usr/bin/env python3
"""Freeze the current-head Q3 full-frame primitive/path-domain join."""

import argparse
import hashlib
import json
from pathlib import Path

from full_frame import DEPENDENCIES, review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    result = review()
    dependency_paths = [HERE.parent / name / "evidence.json" for name in DEPENDENCIES]
    own = [HERE / name for name in (
        "README.md", "full_frame.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    scope = result["scope"]
    return {
        "status": "current-head-Q3-full-frame-whole-path-and-phase-domain-bound",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": result,
        "gates": {
            "nine_exact_graph_and_interface_controls": len(result["checks"]) == 9 and all(result["checks"].values()),
            "seventeen_direct_provider_records": len(dependency_paths) == 17 and all(path.exists() for path in dependency_paths),
            "scientific_scope_preserved": scope["Q1_Q2_object_is_a_direct_dependency"] and scope["Q3_full_frame_primitive_path_and_phase_domain_bound_on_same_later_tail"] and not scope["old_numeric_Cphi_sufficient_gate_used"] and not scope["unspecified_historical_first_band_certified"] and not scope["every_provider_theorem_independently_reproved_here"] and not scope["fixed_projector_compensator_and_weighted_inverse_contracts_closed"] and not scope["hashes_or_tests_treated_as_continuum_proof"] and not scope["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"],
        },
        "dependency_evidence_sha256": {
            str(path.relative_to(PROJECT)): sha(path) for path in dependency_paths
        },
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "bind the fixed compensator/projector profiles, actual supports and every fixed jet to this Q1-Q3 field and accepted later tail",
            "join the compact and slow inverse maps with their weighted domains and whole-path support on the same field",
            "then join generated harmonic sources, complete covariance and fresh axial/full-pressure means, retaining Pi1/Pi2 and all reserved supports",
            "only after those contracts may the full viscosity/cutoff/interaction residual, force, comparison or blow-up claims be reconsidered",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Q3 graph, frame/path interfaces, provenance or scope changed")


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
    print("current-head Q3 full-frame join:", graph["nodes"], "nodes;", graph["edges"], "edges;", graph["source_hash_checks"], "source checks")


if __name__ == "__main__":
    main()
