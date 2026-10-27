#!/usr/bin/env python3
"""Freeze log-reserved pressure recovery and expanding-ball comparison."""

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from pressure_comparison import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_log_reserved_preterminal_energy_replay",
    "core_log_reserved_localization_force_replay",
    "core_corrected_comparison_join",
    "comparison_foundation",
    "core_same_field_positive_core_reconciliation",
    "core_complete_dependency_quantifier_audit",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    if isinstance(value, Fraction):
        return str(value)
    return value


def generate():
    result = review()
    dependencies = {
        str((HERE.parent / name / "evidence.json").relative_to(PROJECT)):
        sha(HERE.parent / name / "evidence.json")
        for name in DEPENDENCIES
    }
    own = [HERE / name for name in (
        "README.md", "pressure_comparison.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    scope = result["scope"]
    return {
        "status": "log-reserved-pressure-and-expanding-ball-comparison-replayed-conditionally",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": encode(result),
        "gates": {
            "ten_exact_controls": len(result["checks"]) == 10 and all(result["checks"].values()),
            "six_frozen_dependencies": len(dependencies) == 6,
            "conditional_comparison_joined": (
                scope["fixed_gauge_pressure_gradient_recovered_conditionally"]
                and scope["same_field_expanding_ball_comparison_replayed_conditionally"]
                and scope["global_smooth_finite_energy_competitor_exclusion_follows_conditionally"]
            ),
            "scientific_scope_preserved": (
                not scope["finite_energy_hypothesis_removed"]
                and not scope["mixed_pressure_moment_isotropy_assumed"]
                and not scope["future_velocity_or_pressure_constructed"]
                and not scope["all_frozen_continuum_contracts_independently_reproved_here"]
                and not scope["unconditional_admissible_force_or_blowup_proved"]
            ),
        },
        "dependency_evidence_sha256": dependencies,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "audit the current-head dependency graph to enumerate every still-conditional continuum contract",
            "independently discharge any remaining primitive or operator contract before promoting admissible-force or exclusion statements",
            "do not claim an unconditional Navier-Stokes blow-up theorem from executable or numerical evidence",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Pressure/comparison replay, provenance, or scope changed")


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
    print("pressure/comparison replay:", len(record["review"]["checks"]),
          "controls; unconditional promotion remains barred")


if __name__ == "__main__":
    main()
