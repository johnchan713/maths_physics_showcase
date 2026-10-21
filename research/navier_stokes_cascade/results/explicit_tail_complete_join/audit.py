#!/usr/bin/env python3
"""Freeze the explicit-tail to canonical residual and fresh-force join."""

import argparse
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from join_checks import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = "330cdbd8a200e6349bdaf48d652ac91105eac2d5"
STATUS = "conditional-explicit-tail-complete-canonical-residual-and-fresh-force-join"
DEPENDENCIES = (
    "explicit_tilt_pulse_tail",
    "canonical_value_domain",
    "canonical_primary_foundation",
    "canonical_initial_mean_join",
    "source_envelope_edge_audit",
    "current_target_cycle_bridge",
    "current_residual_ledger_audit",
    "canonical_local_realization",
    "canonical_endpoint_join",
    "certified_tail_reconstruction",
    "canonical_comparison_join",
)
FLAGS = {
    "one_later_tail_dominates_explicit_primary_canonical_and_edge_value_gates": True,
    "inherited_and_canonical_primary_fields_or_roots_identified": False,
    "canonical_primary_uses_own_reference_frame_entrance_roots_and_current_means": True,
    "all_six_covariances_axial_velocity_complete_pressure_and_Pi1_Pi2_retained": True,
    "current_signed_fast_slow_targets_recomputed_after_every_update": True,
    "generated_regular_sources_have_reviewed_full_path_edge_and_entrance_inputs_conditionally": True,
    "complete_finite_residual_keeps_viscosity_cutoffs_both_crosses_squares_and_tails": True,
    "fresh_realization_localization_Borel_force_energy_and_same_force_comparison_joined": True,
    "historical_first_band_realization_force_and_schedules_reidentified": False,
    "fixed_stage_estimate_promoted_to_uncontrolled_stage_diagonal": False,
    "retained_core_background_defining_map_and_inverse_premises_all_discharged": False,
    "analytical_join_distinguished_from_exact_finite_diagnostics": True,
    "unconditional_pulse_complete_residual_admissible_blowup_or_blowup_proved": False,
    "formal_certificate": False,
}
REMAINING = [
    "independently audit and discharge the retained exact core background and global defining-map existence inverse gap weighted-edge and finite-contraction hypotheses",
    "do not identify the distinct canonical field and fresh force with the uncertified historical first-band realization and force",
    "do not promote conditional continuum lemmas or exact finite tests to an unconditional Navier-Stokes blow-up theorem",
]


def protocol():
    return {
        "parent_commit": PARENT,
        "tail": "ell_sharp=max(ell_explicit,ell_canonical,ell_edge); q_sharp=min(qstar/64,2^(-ell_sharp-2)); one derivative-order- and stage-independent later tail",
        "objects": "inherited and canonical primaries remain distinct; canonical Br entrance Hc ac and current means are rebuilt, never copied by O(1/L) closeness",
        "covariance": "all six complete-curl covariance entries, axial velocity, compact pressure, whole dZ(c_rho P), Pi1/Pi2 and actual A0/A1 remainders retained",
        "sources": "current generated regular expressions retain f/P full-path jets, Gaussian envelope, both flat edges, zero entrance collars and complete pressure-normal inverse",
        "residual": "exact angular complement, cutoff curl, temporal additive ledger, viscosity, both old-new orientations, squares and current mean-wave interactions; finite gain 1/2-kappa",
        "force": "fresh whole-interval realization diagonal, certified localization, all twelve force groups, fresh Borel scales, energy and same-force comparison",
        "scope": "conditional compatibility join; historical first band/field/force remain uncertified; retained continuum foundation bundle remains; no unconditional theorem",
    }


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    if type(value) is F:
        return f"{value.numerator}/{value.denominator}"
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    return value


def records():
    return {
        name: json.loads((HERE.parent / name / "evidence.json").read_text())
        for name in DEPENDENCIES
    }


def dependency_review():
    rows = records()
    checks = {
        "explicit_parent_supplies_new_tilt_and_tail": rows["explicit_tilt_pulse_tail"]["flags"]["selected_cutoff_mixed_gate_has_finite_later_tail"],
        "canonical_value_domain_keeps_one_domain_for_all_jets": rows["canonical_value_domain"]["flags"]["one_value_domain_for_every_fixed_jet"],
        "canonical_primary_is_not_inherited_primary": not rows["canonical_primary_foundation"]["flags"]["canonical_and_inherited_physical_primaries_or_roots_exactly_equal"],
        "canonical_initial_mean_keeps_all_six_covariances": rows["canonical_initial_mean_join"]["flags"]["canonical_real_complete_curl_all_six_covariance_inputs_joined"],
        "source_grammar_has_full_path_inputs": rows["source_envelope_edge_audit"]["flags"]["generated_regular_source_grammar_has_full_path_and_entrance_inputs_relative_to_current_target_estimates"],
        "current_targets_use_exact_state_differences": rows["current_target_cycle_bridge"]["flags"]["current_target_estimate_derived_by_exact_state_difference_not_raw_field_bound"],
        "complete_finite_residual_is_conditional": rows["current_residual_ledger_audit"]["flags"]["conditional_complete_finite_residual_estimate_relative_to_listed_inputs"],
        "canonical_realization_keeps_complete_residual": rows["canonical_local_realization"]["flags"]["all_fixed_Cartesian_jets_of_complete_residual_flat_relative_to_retained_lemmas"],
        "endpoint_keeps_full_pressure_memory": rows["canonical_endpoint_join"]["flags"]["all_pressure_memory_axial_velocity_and_cutoff_curl_terms_retained"],
        "distinct_tail_has_fresh_force": rows["certified_tail_reconstruction"]["flags"]["fresh_compact_force_energy_and_same_force_comparison_joined_conditionally"],
        "comparison_uses_same_canonical_field_and_force": rows["canonical_comparison_join"]["flags"]["same_canonical_complete_field_force_zero_datum_and_pressure_used"],
        "all_dependencies_retain_remaining_obligations": all(
            rows[name]["remaining_obligations"] for name in DEPENDENCIES
        ),
    }
    return {"checks": checks}


def dependencies():
    return {
        str((HERE.parent / name / "evidence.json").relative_to(PROJECT)): sha(
            HERE.parent / name / "evidence.json"
        )
        for name in DEPENDENCIES
    }


def sources():
    paths = [HERE / name for name in (
        "README.md", "join_checks.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    paths += [
        HERE.parent / name / filename
        for name in DEPENDENCIES
        for filename in ("README.md", "evidence.json")
    ]
    return {str(path.relative_to(PROJECT)): sha(path) for path in paths}


def generate():
    if json.loads((HERE / "protocol.json").read_text()) != protocol():
        raise ValueError("Protocol changed")
    exact = encode(review())
    routed = dependency_review()
    gates = {
        "exact_join_checks": all(exact["checks"].values()),
        "dependency_routes": all(routed["checks"].values()),
        "eleven_frozen_dependency_records_bound": len(dependencies()) == 11,
        "all_scope_flags_are_strict_booleans": all(type(value) is bool for value in FLAGS.values()),
    }
    return {
        "status": STATUS,
        "protocol": protocol(),
        "flags": deepcopy(FLAGS),
        "reviews": {"exact_join": exact, "dependency_routing": routed},
        "gates": gates,
        "dependency_evidence_sha256": dependencies(),
        "source_hashes": sources(),
        "remaining_obligations": REMAINING[:],
    }


def validate(record):
    expected = generate()
    if json.dumps(expected, sort_keys=True) != json.dumps(record, sort_keys=True):
        raise ValueError("Explicit-tail complete-join record changed")
    if not all(expected["gates"].values()):
        raise ValueError("Explicit-tail complete-join gate failed")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify-record", type=Path)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit("Use a fresh output; frozen evidence is never overwritten")
    record = generate()
    validate(record)
    if args.verify_record:
        validate(json.loads(args.verify_record.read_text()))
    with args.output.open("x") as stream:
        stream.write(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(
        f"{STATUS}: {len(record['gates'])} gates; explicit later tail joined "
        "conditionally to canonical covariance, complete residual and fresh force."
    )


if __name__ == "__main__":
    main()
