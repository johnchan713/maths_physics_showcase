#!/usr/bin/env python3
"""Freeze the direct compact/slow provider reconstruction."""

import argparse
import hashlib
import json
from pathlib import Path

from provider_reconstruction import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_complete_dependency_quantifier_audit",
    "core_fixed_profile_projector_join",
    "core_compact_slow_inverse_domain_join",
    "core_certified_finite_stage_domain_join",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    result = review()
    dependencies = {
        str((HERE.parent / name / "evidence.json").relative_to(PROJECT)):
        sha(HERE.parent / name / "evidence.json")
        for name in DEPENDENCIES
    }
    own = [HERE / name for name in (
        "README.md", "provider_reconstruction.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    return {
        "status": "direct-compact-slow-operators-proved-with-actual-source-supremum-frontier",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": result,
        "gates": {
            "nine_exact_controls": len(result["checks"]) == 9 and all(result["checks"].values()),
            "four_frozen_dependencies": len(dependencies) == 4,
            "direct_operator_reconstruction": result["scope"]["compact_and_slow_operator_bounds_reconstructed_directly"],
            "scientific_scope_preserved": (
                not result["scope"]["actual_generated_source_tail_supremum_proved"]
                and not result["scope"]["unconditional_complete_residual_force_comparison_or_blowup_claimed"]
            ),
        },
        "dependency_evidence_sha256": dependencies,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "prove the normalized weighted source supremum uniformly over the accepted band tail for every fixed stage and requested finite jet",
            "derive those bounds from the actual matching-stage covariance, axial velocity, moving pressure and support formulas rather than a provider flag",
            "then insert the direct compact and slow estimates into the six-group complete residual and replay localization, force and comparison",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Direct provider reconstruction, provenance, or scope changed")


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
    print("direct compact/slow provider reconstruction:", len(record["review"]["checks"]), "controls; actual source supremum remains open")


if __name__ == "__main__":
    main()
