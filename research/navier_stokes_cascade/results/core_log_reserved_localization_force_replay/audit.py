#!/usr/bin/env python3
"""Freeze the log-reserved localization and compact-force replay."""

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from localization_force import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_log_reserved_endpoint_compatibility_join",
    "core_log_reserved_cartesian_corner_replay",
    "core_log_reserved_complete_residual_replay",
    "core_corrected_force_join",
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
        "README.md", "localization_force.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    scope = result["scope"]
    return {
        "status": "log-reserved-localization-and-compact-force-replayed-conditionally",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": encode(result),
        "gates": {
            "ten_exact_controls": len(result["checks"]) == 10 and all(result["checks"].values()),
            "five_frozen_dependencies": len(dependencies) == 5,
            "localization_and_force_joined": (
                scope["same_field_axisymmetric_localization_replayed_conditionally"]
                and scope["one_smooth_compact_force_extension_constructed_conditionally"]
            ),
            "scientific_scope_preserved": (
                not scope["terminal_velocity_or_pressure_regular_at_origin_proved"]
                and not scope["future_velocity_or_pressure_constructed"]
                and not scope["preterminal_energy_replayed_here"]
                and not scope["pressure_recovery_comparison_or_blowup_proved"]
                and not scope["unconditional_pulse_residual_or_force_claimed"]
            ),
        },
        "dependency_evidence_sha256": dependencies,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "replay the fixed-time preterminal energy identity for this same localized field and compact force",
            "replay fixed-gauge pressure recovery including force divergence",
            "replay expanding-ball comparison before any competitor-exclusion or blow-up claim",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Localization/force replay, provenance, or scope changed")


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
    print("localization/force replay:", len(record["review"]["checks"]),
          "controls; energy and comparison remain open")


if __name__ == "__main__":
    main()
