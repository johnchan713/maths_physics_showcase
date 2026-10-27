#!/usr/bin/env python3
"""Freeze the corrected closed-time endpoint compatibility join."""

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from endpoint_compatibility import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_log_reserved_cartesian_corner_replay",
    "core_log_reserved_complete_residual_replay",
    "core_corrected_endpoint_join",
    "endpoint_foundation",
    "current_endpoint_jets",
    "core_corrected_force_join",
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
        "README.md", "endpoint_compatibility.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    return {
        "status": "log-reserved-corner-and-punctured-endpoint-inputs-joined-conditionally",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": encode(result),
        "gates": {
            "ten_exact_controls": len(result["checks"]) == 10 and all(result["checks"].values()),
            "six_frozen_dependencies": len(dependencies) == 6,
            "closed_time_inputs": result["scope"]["closed_time_localization_inputs_joined_conditionally"],
            "scientific_scope_preserved": (
                not result["scope"]["terminal_velocity_or_pressure_regular_at_origin_proved"]
                and not result["scope"]["localization_and_compact_force_replayed_here"]
                and not result["scope"]["unconditional_pulse_residual_force_comparison_or_blowup_claimed"]
            ),
        },
        "dependency_evidence_sha256": dependencies,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "replay the exact axisymmetric localization formula with all twelve residual terms and unit viscosity",
            "construct and glue the compact smooth future force from the resulting terminal force jets",
            "then replay preterminal energy, fixed-gauge pressure recovery and expanding-ball comparison",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Endpoint compatibility, provenance, or scope changed")


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
    print("endpoint compatibility join:", len(record["review"]["checks"]),
          "controls; localization and compact force remain open")


if __name__ == "__main__":
    main()
