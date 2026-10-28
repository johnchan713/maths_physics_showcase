#!/usr/bin/env python3
"""Freeze the current-chain log-reserved six-group residual reconciliation."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from residual import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "current-chain-log-reserved-six-group-local-residual-reconciled-conditionally"
DEPENDENCIES = (
    "core_current_log_source_reconciliation",
    "core_log_reserved_complete_residual_replay",
    "core_compact_slow_inverse_reconstruction",
    "core_log_renormalized_actual_source_induction",
    "core_finite_stage_compact_recurrence",
    "core_complete_potential_diagonal_join",
    "core_corner_tail_slack_reconciliation",
    "core_post_signed_residual_join",
)
FLAGS = {
    "current_source_chain_and_residual_use_same_objects": True,
    "extra_A0_A1_integration_derivative_cost_carried_into_residual_stage": True,
    "compact_regular_and_cutoff_log_reserves_are_disjoint": True,
    "time_viscosity_pressure_both_crosses_and_new_square_retained": True,
    "complete_covariance_axial_velocity_and_pressure_memory_retained": True,
    "same_locally_finite_field_supports_gauges_and_cutoff_sequence_used": True,
    "stage_uniform_constants_or_fixed_q_contraction_proved": False,
    "corner_endpoint_force_energy_or_comparison_replayed_here": False,
    "unconditional_pulse_residual_force_exclusion_or_blowup_proved": False,
}
REMAINING = [
    "instantiate the corrected Cartesian-corner and punctured-endpoint residual transfer on this exact current field",
    "then revalidate all twelve localization terms, compact future force, preterminal energy and fixed-gauge expanding-ball comparison",
    "retain the unresolved frozen exact core/root and regional/frame hypotheses before any unconditional promotion",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dependency_review():
    records, sources = {}, {}
    for name in DEPENDENCIES:
        path = HERE.parent / name / "evidence.json"
        record = json.loads(path.read_text())
        records[str(path.relative_to(PROJECT))] = sha(path)
        for relative, digest in record.get("source_hashes", {}).items():
            source = PROJECT / relative
            sources[name + ":" + relative] = source.exists() and sha(source) == digest
    return records, sources


def generate():
    records, inherited = dependency_review()
    result = review()
    own = [HERE / name for name in (
        "README.md", "residual.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    providers = [
        HERE.parent / "core_current_log_source_reconciliation/reconcile.py",
        HERE.parent / "core_log_reserved_complete_residual_replay/replay.py",
        HERE.parent / "core_complete_potential_diagonal_join/diagonal.py",
    ]
    gates = {
        "eight_exact_current_residual_controls": len(result["checks"]) == 8 and all(result["checks"].values()),
        "eight_frozen_dependencies": len(records) == len(DEPENDENCIES),
        "all_declared_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "scientific_scope_retained": (
            result["scope"]["six_group_log_reserved_local_residual_instantiated_on_current_chain"]
            and not result["scope"]["stage_uniform_constants_or_fixed_q_contraction_used"]
            and not result["scope"]["corner_endpoint_force_energy_comparison_replayed_here"]
            and not result["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"]
        ),
    }
    return {
        "status": STATUS,
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "flags": deepcopy(FLAGS),
        "review": result,
        "gates": gates,
        "dependency_evidence_sha256": records,
        "inherited_source_checks": inherited,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own + providers},
        "remaining_obligations": list(REMAINING),
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Current residual reconciliation, provenance or scientific scope changed")


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
    print(
        "current residual reconciliation:", len(record["review"]["checks"]),
        "exact controls;", len(record["inherited_source_checks"]),
        "dependency source checks",
    )


if __name__ == "__main__":
    main()
