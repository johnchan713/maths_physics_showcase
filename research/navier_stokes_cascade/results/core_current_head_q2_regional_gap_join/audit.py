#!/usr/bin/env python3
"""Freeze the current-head Q2 regional continuum and factored-edge join."""

import argparse
import hashlib
import json
from pathlib import Path

from regional_gaps import DEPENDENCIES, review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    result = review()
    dependency_paths = [HERE.parent / name / "evidence.json" for name in DEPENDENCIES]
    own = [HERE / name for name in (
        "README.md", "regional_gaps.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    scope = result["scope"]
    return {
        "status": "current-head-Q2-regional-continuum-factored-edges-and-same-band-gap-bound",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": result,
        "gates": {
            "nine_exact_graph_and_interface_controls": len(result["checks"]) == 9 and all(result["checks"].values()),
            "thirteen_direct_provider_records": len(dependency_paths) == 13 and all(path.exists() for path in dependency_paths),
            "scientific_scope_preserved": scope["Q1_exact_object_is_a_direct_dependency"] and scope["Q2_regional_continuum_and_factored_edge_chain_bound_on_that_object"] and scope["Q2_positive_same_band_radius_rechecked"] and not scope["every_provider_theorem_independently_reproved_here"] and not scope["positive_absolute_stress_infimum_at_flat_edges_claimed"] and not scope["historical_uncertified_existing_tilt_gate_used"] and not scope["Q3_full_frame_primitive_domain_discharged"] and not scope["hashes_tests_or_interval_diagnostics_treated_as_continuum_proof"] and not scope["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"],
        },
        "dependency_evidence_sha256": {
            str(path.relative_to(PROJECT)): sha(path) for path in dependency_paths
        },
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "audit Q3 on this Q1-Q2 object: the original full-frame primitive, moving-normal path ODE, complete phase jets and accepted physical domain on the same later tail",
            "then audit fixed compensator/projector supports and compact/slow inverse weighted domains before promoting generated sources",
            "retain complete covariance means including axial velocity, curl remainders, Pi1/Pi2 and every positive-order, wave, mean, localization and forced pressure component",
            "only after those obligations may complete viscosity/cutoff/interaction residual, force, exclusion or blow-up claims be reconsidered",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Q2 graph, continuum interfaces, provenance or scope changed")


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
    print("current-head Q2 regional-gap join:", graph["nodes"], "nodes;", graph["edges"], "edges;", graph["source_hash_checks"], "source checks")


if __name__ == "__main__":
    main()
