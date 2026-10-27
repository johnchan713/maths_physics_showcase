#!/usr/bin/env python3
"""Freeze the complete dependency and quantifier audit."""

import argparse
import hashlib
import json
from pathlib import Path

from dependency_graph import ROOTS, review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    result = review()
    dependencies = {
        str((HERE.parent / name / "evidence.json").relative_to(PROJECT)):
        sha(HERE.parent / name / "evidence.json")
        for name in ROOTS
    }
    own = [HERE / name for name in (
        "README.md", "dependency_graph.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    return {
        "status": "complete-common-object-quantifier-audited-with-provider-frontier",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": result,
        "gates": {
            "seven_exact_controls": len(result["checks"]) == 7 and all(result["checks"].values()),
            "thirteen_direct_frozen_dependencies": len(dependencies) == len(ROOTS),
            "common_object_and_quantifier_join": result["scope"]["common_object_and_quantifier_join_checked"],
            "scientific_scope_preserved": (
                not result["scope"]["tests_or_hashes_treated_as_continuum_proof"]
                and not result["scope"]["unconditional_pulse_residual_force_or_blowup_claimed"]
            ),
        },
        "dependency_evidence_sha256": dependencies,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "independently reconstruct the master provider condition for one fixed all-order profile and every fixed finite-stage representative",
            "then replay the complete residual and localized force from those reconstructed estimates rather than from conditional provider flags",
            "only after that replay use the already valid fixed-T then radius then protected-core limit to assess unconditional competitor exclusion",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Complete dependency graph, quantifiers, provenance, or scope changed")


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
    with args.output.open("x") as output:
        output.write(json.dumps(record, indent=2, sort_keys=True) + "\n")
    graph = record["review"]["graph"]
    print("complete dependency audit:", graph["nodes"], "nodes;", graph["edges"], "edges;", graph["source_hash_checks"], "source checks")


if __name__ == "__main__":
    main()
