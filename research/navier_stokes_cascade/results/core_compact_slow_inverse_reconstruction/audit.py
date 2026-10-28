#!/usr/bin/env python3
"""Freeze the compact/slow inverse reconstruction on the current same-choice tail."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from inverse import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "compact-slow-weighted-inverse-contract-reconstructed-on-current-same-choice-tail"
DEPENDENCIES = (
    "core_fixed_profile_projector_reconstruction",
    "core_compact_slow_inverse_domain_join",
    "core_rho_projector_signed_recompute",
    "core_initial_haar_row_budget",
    "core_native_haar_primitive_join",
    "core_symbolic_weighted_maps",
    "core_fast_inverse_budget",
    "core_slow_inverse_budget",
    "core_real_refresh_obligation_join",
    "core_initial_final_row_obstruction",
    "core_finite_support_closure",
    "core_direct_compact_slow_provider_reconstruction",
    "core_log_renormalized_actual_source_induction",
)
FLAGS = {
    "current_fixed_profile_and_same_choice_projector_tail_used": True,
    "shifted_compact_maps_joined_on_complete_weighted_domain": True,
    "coefficient_differentiation_introduces_no_shift_parameter_loss": True,
    "zero_haar_temporal_inverse_and_actual_A1_velocity_retained": True,
    "fixed_slow_response_joined_on_complete_reserved_potential_hull": True,
    "whole_moving_pressure_product_and_both_pressure_tails_retained": True,
    "all_fixed_output_request_graphs_are_finite": True,
    "whole_path_and_zero_extension_support_contract_joined": True,
    "old_pure_power_zero_margin_identity_retained": True,
    "old_pure_power_zero_margin_is_current_global_frontier": False,
    "later_log_renormalized_actual_source_route_identified": True,
    "later_log_source_induction_or_residual_replay_reproved_here": False,
    "uniformity_over_derivative_order_stage_width_or_tilt_proved": False,
    "unconditional_pulse_force_residual_improvement_or_blowup_proved": False,
}
REMAINING = [
    "reconcile the reconstructed primitive contracts on one current object with the frozen log-renormalized actual-source induction and its extra compact-integration reserve",
    "replay complete covariance, fresh axial velocity, full moving-pressure means, viscosity, cutoffs and all interactions on those exact reconstructed objects",
    "only after that same-object replay, revalidate localization, compact force, energy and fixed-gauge comparison before any exclusion or blow-up claim",
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
    protocol = json.loads((HERE / "protocol.json").read_text())
    own = [HERE / name for name in (
        "README.md", "inverse.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    providers = [
        HERE.parent / "core_symbolic_weighted_maps/README.md",
        HERE.parent / "core_fast_inverse_budget/README.md",
        HERE.parent / "core_slow_inverse_budget/README.md",
        HERE.parent / "core_initial_final_row_obstruction/README.md",
        HERE.parent / "core_direct_compact_slow_provider_reconstruction/README.md",
        HERE.parent / "core_log_renormalized_actual_source_induction/README.md",
    ]
    gates = {
        "all_exact_domain_request_and_route_controls": bool(result["checks"]) and all(result["checks"].values()),
        "thirteen_frozen_dependencies": len(records) == len(DEPENDENCIES),
        "all_declared_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "scientific_scope_retained": (
            FLAGS["current_fixed_profile_and_same_choice_projector_tail_used"]
            and FLAGS["shifted_compact_maps_joined_on_complete_weighted_domain"]
            and FLAGS["old_pure_power_zero_margin_identity_retained"]
            and not FLAGS["old_pure_power_zero_margin_is_current_global_frontier"]
            and FLAGS["later_log_renormalized_actual_source_route_identified"]
            and not FLAGS["later_log_source_induction_or_residual_replay_reproved_here"]
            and not FLAGS["unconditional_pulse_force_residual_improvement_or_blowup_proved"]
        ),
    }
    return {
        "status": STATUS,
        "protocol": protocol,
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
        raise ValueError("Compact/slow inverse reconstruction, provenance, route classification or scope changed")


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
        "compact/slow inverse reconstruction:", len(record["review"]["checks"]),
        "exact controls;", len(record["inherited_source_checks"]),
        "dependency source checks",
    )


if __name__ == "__main__":
    main()
