#!/usr/bin/env python3
"""Freeze the log-reserved complete residual replay."""

import argparse
import hashlib
import json
from pathlib import Path

from replay import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_log_renormalized_actual_source_induction",
    "core_post_signed_residual_join",
    "core_finite_stage_compact_recurrence",
    "core_corner_tail_slack_reconciliation",
    "core_complete_potential_diagonal_join",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    if hasattr(value, "numerator") and hasattr(value, "denominator"):
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
        "README.md", "replay.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    return {
        "status": "log-reserved-complete-six-group-local-residual-replayed-conditionally",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": encode(result),
        "gates": {
            "ten_exact_controls": len(result["checks"]) == 10 and all(result["checks"].values()),
            "five_frozen_dependencies": len(dependencies) == 5,
            "six_group_local_residual": result["scope"][
                "log_reserved_six_group_complete_local_residual_flat_conditionally"
            ],
            "scientific_scope_preserved": (
                not result["scope"]["stage_uniform_constants_or_fixed_q_contraction_used"]
                and not result["scope"]["unconditional_pulse_force_comparison_or_blowup_claimed"]
                and not result["scope"]["localization_force_energy_and_comparison_replayed_after_source_correction"]
            ),
        },
        "dependency_evidence_sha256": dependencies,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "replay the punctured closed-corner endpoint and Cartesian flatness transfer on the corrected same field",
            "replay localization and compact forcing including viscosity and cutoff derivatives",
            "replay energy, fixed-gauge pressure recovery and expanding-ball comparison before any blow-up claim",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Log-reserved residual replay, provenance, or scope changed")


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
    print("log-reserved residual replay:", len(record["review"]["checks"]),
          "controls; localization/force/comparison replay remains open")


if __name__ == "__main__":
    main()
