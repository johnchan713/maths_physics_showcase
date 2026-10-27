#!/usr/bin/env python3
"""Freeze the log-reserved Cartesian-corner residual replay."""

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from corner_replay import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_log_reserved_complete_residual_replay",
    "core_log_renormalized_actual_source_induction",
    "core_corner_flatness_obstruction",
    "core_corner_tail_slack_reconciliation",
    "core_primitive_margin_normalization_reconciliation",
    "core_closed_chart_coefficient_join",
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
        "README.md", "corner_replay.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    return {
        "status": "log-reserved-six-group-residual-transferred-to-cartesian-corner-conditionally",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": encode(result),
        "gates": {
            "ten_exact_controls": len(result["checks"]) == 10 and all(result["checks"].values()),
            "six_frozen_dependencies": len(dependencies) == 6,
            "cartesian_corner_flatness": result["scope"][
                "same_field_all_fixed_residual_jets_cartesian_corner_flat_conditionally"
            ],
            "scientific_scope_preserved": (
                not result["scope"]["punctured_field_endpoint_jets_reproved_here"]
                and not result["scope"]["localization_force_energy_comparison_replayed_here"]
                and not result["scope"]["unconditional_pulse_residual_force_or_blowup_claimed"]
            ),
        },
        "dependency_evidence_sha256": dependencies,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "join the inherited punctured one-sided field and pressure jets to this corrected Cartesian residual flatness",
            "replay localization and compact forcing with the extra time derivative, viscosity and every cutoff derivative",
            "then replay energy, fixed-gauge pressure recovery and expanding-ball comparison before any blow-up claim",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Cartesian-corner replay, provenance, or scope changed")


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
    print("log-reserved Cartesian corner:", len(record["review"]["checks"]),
          "controls; localization and force replay remain open")


if __name__ == "__main__":
    main()
