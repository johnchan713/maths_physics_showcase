#!/usr/bin/env python3
"""Freeze the fixed-profile/projector reconstruction on the current tail."""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from projectors import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "fixed-profile-same-choice-projector-contract-reconstructed-on-current-full-path-tail"
DEPENDENCIES = (
    "core_full_path_ode_reconstruction",
    "core_fixed_profile_projector_join",
    "core_fixed_profile_tame",
    "core_moving_pressure_compensator_join",
    "core_reserved_geometry_recheck",
    "core_log_pullback_projector_join",
    "core_rho_squared_projector_realization",
    "core_rho_projector_signed_recompute",
    "core_finite_support_closure",
)
FLAGS = {
    "current_full_path_tail_used": True,
    "one_fixed_smooth_pressure_profile_with_profile_dependent_constants_used": True,
    "normalization_only_projector_jet_inference_rejected": True,
    "explicit_new_same_profile_rho_squared_projectors_selected": True,
    "both_exact_projector_moment_normalizations_proved": True,
    "actual_support_strictly_contained_in_moving_enclosure": True,
    "moving_enclosure_identified_with_actual_support": False,
    "fixed_order_band_uniform_jets_joined": True,
    "reserved_weight_floor_used_only_on_actual_moving_support": True,
    "same_choice_source_signed_target_covariance_and_pressure_recomputation_required": True,
    "historical_unspecified_projectors_or_targets_identified_with_new_choice": False,
    "uniformity_over_all_profiles_orders_or_cycles_proved": False,
    "compact_or_slow_weighted_inverse_domains_closed_here": False,
    "actual_high_jet_current_fields_or_complete_residual_certified": False,
    "unconditional_pulse_force_residual_improvement_or_blowup_proved": False,
}
REMAINING = [
    "join the compact and slow inverse maps to these same-choice coefficients on their full weighted domains and whole-path supports",
    "supply the actual linked current-field high-jet inputs and strict epsilon margins at every requested finite stage without derivative-order uniformity",
    "join fresh axial/full-pressure means, generated interactions, viscosity and every cutoff term in the complete residual",
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
        "README.md", "projectors.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    providers = [
        HERE.parent / "core_fixed_profile_tame/proof_checks.py",
        HERE.parent / "core_moving_pressure_compensator_join/bounds.py",
        HERE.parent / "core_rho_squared_projector_realization/bounds.py",
        HERE.parent / "core_rho_projector_signed_recompute/bounds.py",
        HERE.parent / "core_reserved_geometry_recheck/proof_checks.py",
    ]
    gates = {
        "all_exact_quantifier_support_and_identity_controls": bool(result["checks"]) and all(result["checks"].values()),
        "nine_frozen_dependencies": len(records) == len(DEPENDENCIES),
        "all_declared_dependency_sources_match": bool(inherited) and all(inherited.values()),
        "scientific_scope_retained": (
            FLAGS["normalization_only_projector_jet_inference_rejected"]
            and FLAGS["current_full_path_tail_used"]
            and FLAGS["actual_support_strictly_contained_in_moving_enclosure"]
            and not FLAGS["moving_enclosure_identified_with_actual_support"]
            and FLAGS["same_choice_source_signed_target_covariance_and_pressure_recomputation_required"]
            and not FLAGS["historical_unspecified_projectors_or_targets_identified_with_new_choice"]
            and not FLAGS["compact_or_slow_weighted_inverse_domains_closed_here"]
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
        raise ValueError("Fixed-profile/projector contract, provenance or scope changed")


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
        "fixed-profile/projector reconstruction:", len(record["review"]["checks"]),
        "exact controls;", len(record["inherited_source_checks"]),
        "dependency source checks",
    )


if __name__ == "__main__":
    main()
