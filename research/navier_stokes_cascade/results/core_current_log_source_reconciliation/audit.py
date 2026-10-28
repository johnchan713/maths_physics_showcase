#!/usr/bin/env python3
"""Freeze current-chain instantiation of the log-renormalized source induction."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from reconcile import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "current-five-contract-log-renormalized-actual-source-induction-reconciled-conditionally"
DEPENDENCIES = (
    "core_current_head_conditional_frontier",
    "core_leading_background_seminorm_reconstruction",
    "core_leading_profile_certificate",
    "core_same_band_gap_reconstruction",
    "core_full_path_ode_reconstruction",
    "core_fixed_profile_projector_reconstruction",
    "core_compact_slow_inverse_reconstruction",
    "core_log_renormalized_actual_source_induction",
    "core_finite_stage_compact_recurrence",
)
FLAGS = {
    "one_current_profile_tail_path_projector_and_inverse_chain_used": True,
    "fixed_current_primitive_leaves_enter_with_log_degree_zero": True,
    "actual_seed_and_sourced_generic_formula_graphs_retained": True,
    "extra_compact_integration_costs_two_source_derivatives_explicitly": True,
    "all_fixed_enlarged_requests_are_accepted": True,
    "six_covariance_axial_velocity_and_complete_pressure_memory_retained": True,
    "pure_power_zero_margin_reintroduced_as_current_frontier": False,
    "stage_uniform_constants_or_log_degrees_proved": False,
    "complete_residual_or_downstream_replays_reproved_here": False,
    "unconditional_pulse_residual_force_exclusion_or_blowup_proved": False,
}
REMAINING = [
    "instantiate the separately reserved six-group residual replay on these exact current objects, keeping distinct compact-source, regular-residual and cutoff-tail reserve accounts",
    "then revalidate the corrected corner and endpoint transfer, localization, compact force, preterminal energy and fixed-gauge comparison on that same current field",
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
            source_path = PROJECT / relative
            sources[name + ":" + relative] = source_path.exists() and sha(source_path) == digest
    return records, sources


def generate():
    records, inherited = dependency_review()
    result = review()
    own = [HERE / name for name in (
        "README.md", "reconcile.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    providers = [
        HERE.parent / "core_log_renormalized_actual_source_induction/source_induction.py",
        HERE.parent / "core_compact_slow_inverse_reconstruction/inverse.py",
        HERE.parent / "core_finite_stage_compact_recurrence/recurrence.py",
    ]
    gates = {
        "twelve_exact_instantiation_controls": len(result["checks"]) == 12 and all(result["checks"].values()),
        "nine_frozen_dependencies": len(records) == len(DEPENDENCIES),
        "all_declared_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "scientific_scope_retained": (
            result["scope"]["log_source_induction_instantiated_on_current_five_contract_chain"]
            and not result["scope"]["stage_uniform_constants_or_log_degrees_proved"]
            and not result["scope"]["complete_residual_or_downstream_force_comparison_replayed_here"]
            and not result["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"]
            and not FLAGS["pure_power_zero_margin_reintroduced_as_current_frontier"]
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
        raise ValueError("Current log-source reconciliation, provenance or scientific scope changed")


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
        "current log-source reconciliation:", len(record["review"]["checks"]),
        "exact controls;", len(record["inherited_source_checks"]),
        "dependency source checks",
    )


if __name__ == "__main__":
    main()
