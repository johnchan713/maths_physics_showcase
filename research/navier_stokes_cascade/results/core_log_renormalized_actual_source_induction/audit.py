#!/usr/bin/env python3
"""Freeze the log-renormalized actual-source induction."""

import argparse
import hashlib
import json
from pathlib import Path

from source_induction import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_direct_compact_slow_provider_reconstruction",
    "core_primitive_margin_normalization_reconciliation",
    "core_prefast_weighted_source_join",
    "core_next_sourced_mean_tail_join",
    "core_post_signed_mean_closure",
    "core_finite_stage_compact_recurrence",
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
        "README.md", "source_induction.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    return {
        "status": "log-renormalized-actual-source-tail-induction-proved-with-residual-replay-frontier",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": result,
        "gates": {
            "ten_exact_controls": len(result["checks"]) == 10 and all(result["checks"].values()),
            "six_frozen_dependencies": len(dependencies) == 6,
            "actual_source_induction": result["scope"]["log_renormalized_actual_source_supremum_proved_for_every_fixed_stage"],
            "scientific_scope_preserved": (
                not result["scope"]["pure_power_source_supremum_without_spending_slack_proved"]
                and not result["scope"]["six_group_complete_residual_replayed_with_log_reserve"]
                and not result["scope"]["unconditional_force_comparison_or_blowup_claimed"]
            ),
        },
        "dependency_evidence_sha256": dependencies,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "insert the corrected logarithmic source domain into the direct compact and slow estimates",
            "allocate explicit reserve in all six residual groups, including viscosity, cutoffs and both interaction orientations",
            "then replay localization, compact force, energy and fixed-gauge comparison on the same field",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Log-source induction, provenance, or scope changed")


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
    print("log-renormalized source induction:", len(record["review"]["checks"]),
          "controls; six-group residual replay remains open")


if __name__ == "__main__":
    main()
