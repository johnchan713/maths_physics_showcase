#!/usr/bin/env python3
"""Freeze the log-reserved preterminal energy replay."""

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from preterminal_energy import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_log_reserved_localization_force_replay",
    "core_log_reserved_endpoint_compatibility_join",
    "core_corrected_force_join",
    "core_certified_complete_potential_force_join",
    "force_foundation",
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
        "README.md", "preterminal_energy.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    scope = result["scope"]
    return {
        "status": "log-reserved-preterminal-energy-replayed-conditionally",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": encode(result),
        "gates": {
            "ten_exact_controls": len(result["checks"]) == 10 and all(result["checks"].values()),
            "five_frozen_dependencies": len(dependencies) == 5,
            "preterminal_energy_joined": (
                scope["same_field_preterminal_energy_identity_replayed_conditionally"]
                and scope["separate_finite_energy_and_total_dissipation_bounds_obtained_conditionally"]
            ),
            "scientific_scope_preserved": (
                not scope["terminal_velocity_or_pressure_regular_at_origin_proved"]
                and not scope["future_velocity_or_pressure_constructed"]
                and not scope["fixed_gauge_pressure_recovery_replayed_here"]
                and not scope["expanding_ball_comparison_or_competitor_exclusion_proved"]
                and not scope["unconditional_pulse_residual_force_or_blowup_claimed"]
            ),
        },
        "dependency_evidence_sha256": dependencies,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "replay fixed-gauge pressure recovery with the compact force-divergence term retained",
            "prove the expanding-ball comparison for every admissible smooth finite-energy competitor with the same datum and force",
            "only then revisit competitor exclusion; no unconditional blow-up claim is available",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Energy replay, provenance, or scope changed")


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
    print("preterminal energy replay:", len(record["review"]["checks"]),
          "controls; pressure recovery and comparison remain open")


if __name__ == "__main__":
    main()
