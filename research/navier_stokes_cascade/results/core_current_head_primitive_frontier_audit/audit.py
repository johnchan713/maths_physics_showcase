#!/usr/bin/env python3
"""Freeze the current-head lower primitive-frontier audit."""

import argparse
import hashlib
import json
from pathlib import Path

from current_frontier import BOUNDARY_WITNESSES, START, primitive_frontier, review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    result = review()
    dependency = HERE.parent / START / "evidence.json"
    own = [HERE / name for name in (
        "README.md", "current_frontier.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    scope = result["scope"]
    return {
        "status": "current-head-five-contract-frontier-corrected-to-lower-core-root-regional-frame-antecedents",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": result,
        "gates": {
            "seven_graph_frontier_controls": len(result["checks"]) == 7 and all(result["checks"].values()),
            "one_current_head_dependency": dependency.exists(),
            "three_live_primitive_groups_enumerated": len(result["live_frontier"]) == len(primitive_frontier()) == 3,
            "three_boundary_witnesses_reached": len(result["graph"]["boundary_links"]) == len(BOUNDARY_WITNESSES) and all(result["graph"]["boundary_links"].values()),
            "scientific_scope_preserved": scope["former_five_contracts_reconstructed_on_current_objects"] and scope["downstream_replay_complete_relative_to_lower_frontier"] and not scope["Q1_Q2_Q3_independently_reproved_here"] and not scope["hashes_or_tests_treated_as_continuum_proof"] and not scope["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"],
        },
        "dependency_evidence_sha256": {str(dependency.relative_to(PROJECT)): sha(dependency)},
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "prove Q1 from the pressure-family fixed-point map and every exact finite root, including completed backward-pressure identification on the same profile",
            "prove Q2 as continuum infima and factored-edge coefficients on the actual regional cover, without grid sampling or a false positive absolute-stress floor",
            "prove Q3 for the original full frame, complete path and phase jets and accepted physical domain, then evaluate the finite selector rows only after their positive margins and finite constants exist",
            "only then reconsider unconditional promotion of the already reconstructed source, residual, localization, force, energy and comparison implications",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Current-head graph, primitive frontier, provenance or scope changed")


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
    print("current-head primitive frontier audit:", graph["nodes"], "nodes;", graph["edges"], "edges;", graph["source_hash_checks"], "source checks")


if __name__ == "__main__":
    main()
