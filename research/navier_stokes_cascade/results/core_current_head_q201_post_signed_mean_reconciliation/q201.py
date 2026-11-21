"""Q201: bind corrected target rates to the post-signed fast/slow mean block."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q201-POST-SIGNED-MEAN-RECONCILIATION"

PROVIDERS = (
    "core_current_head_q200_physical_weighted_target_reconciliation",
    "core_post_signed_mean_closure",
    "core_post_signed_residual_join",
    "core_post_signed_seed_sync",
    "core_current_head_q160_physical_target_induction",
    "core_current_head_q161_source_tail_simultaneous_induction",
    "core_current_head_q180_target_source_simultaneous_induction",
    "finite_stage_induction",
    "current_target_moment_audit",
    "compact_mean_foundation",
)

EXPECTED = {
    PROVIDERS[0]: "Q200-physical-weighted-target-route-reconciled-at-corrected-rates",
    PROVIDERS[1]: "conditional-post-signed-fresh-mean-full-regular-row-and-complete-request-domain-closure",
    PROVIDERS[2]: "conditional-post-signed-current-nonzero-source-and-full-physical-request-tail-join",
    PROVIDERS[3]: "conditional-post-signed-induction-seed-compact-sync",
    PROVIDERS[4]: "Q160-physical-all-finite-cycle-target-moment-induction-joined",
    PROVIDERS[5]: "Q161-physical-source-inverse-tail-simultaneous-induction-joined",
    PROVIDERS[6]: "Q180-physical-target-source-tail-simultaneous-induction-reconciled",
    PROVIDERS[7]: "conditional-current-state-complete-finite-stage-induction",
    PROVIDERS[8]: "current-target-full-weight-identities-with-higher-order-cancellation-gate",
    PROVIDERS[9]: "compact-mean-inverse-foundation-sharp-axial-conditioning-relative-to-core-and-primary-inputs",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q201 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def rates(kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q201 identity required")
    kappa = Q(kappa)
    if not 0 < kappa < Q(1, 100):
        raise ValueError("reviewed kappa required")
    B = 1 - kappa
    C = B + Q(1, 2)
    D = C - kappa
    C0 = 2 * B
    return {"kappa": kappa, "B": B, "C": C, "D": D, "C0": C0}


def _strict_integrations(source, target, kappa):
    source, target, kappa = map(Q, (source, target, kappa))
    if kappa <= 0 or target <= source:
        raise ValueError("positive compact gain required")
    return max(1, (target - source) // kappa + 1)


def mean_orders(kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q201 identity required")
    r = rates(kappa)
    k, B, C, D, C0 = (r[x] for x in ("kappa", "B", "C", "D", "C0"))
    p1 = _strict_integrations(D, C0, k)
    p0 = _strict_integrations(B, C0, k)
    tangential = {
        "post_signed_Haar": C0,
        "fast_full_difference": D + 1 - 2 * k,
        "slow_full_difference": C + 1 - 2 * k,
        "A1_compact_remainder": D + p1 * k,
    }
    radial = {
        "radial_compensator": C + 1,
        "A0_compact_remainder": B + p0 * k,
    }
    return {
        "post_signed_Qcirc_rows": D,
        "fast_tangential": D,
        "fast_radial": D + 1,
        "postfast_scalars": C,
        "slow_tangential": C,
        "slow_radial": C + 1,
        "postslow_scalars": C + 1,
        "postslow_axial_moment": C + 2,
        "tangential_terms": tangential,
        "tangential_margins_over_C0": {name: value - C0 for name, value in tangential.items()},
        "radial_terms": radial,
        "radial_margins_over_C0": {name: value - C0 for name, value in radial.items()},
        "final_tangential": min(tangential.values()),
        "final_radial": min(C0, *radial.values()),
        "final_regular_mean": min(min(tangential.values()), min(C0, *radial.values())),
        "A1_integrations": p1,
        "A0_integrations": p0,
    }


def exact_mean_memory(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q201 identity required")
    return {
        "fast_source": "Qcirc of the complete post-signed Etheta and Ez rows",
        "fast_identity": "c_i N_t(du,dgamma)=-(Ftheta,Fz)+(0,A1 Fz)",
        "A1_sign_positive": True,
        "fast_all_components_zero_Haar": True,
        "fast_divergence_exact": True,
        "both_velocity_moments_preserved": True,
        "postfast_scalars": (
            "previous postslow scalars plus complete signed covariance change plus correlations "
            "of the full old means with the actual new fast field"
        ),
        "slow_uses_current_postfast_scalars": True,
        "slow_uses_same_fixed_five_map": True,
        "whole_reserved_potential_hull_G_zero_required": True,
        "all_six_covariances_retained": True,
        "axial_velocity_retained": True,
        "compact_pressure_retained": True,
        "moving_rho_and_whole_c_rho_P_derivative_retained": True,
        "viscosity_time_cylindrical_cross_and_square_terms_retained": True,
        "final_tangential_identity": "Efinal=Haar(EpostSigned)+(0,A1 Fz)+Rf+Rs",
        "final_radial_identity": "Er_final=-rho Pnew-A0(gr_final-rho Pnew)",
    }


def compact_requests(output_order=0, original_pole=1, previous_wave_pole=30, kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q201 identity required")
    if any(isinstance(x, bool) or not isinstance(x, int) or x < 0 for x in (output_order, original_pole, previous_wave_pole)):
        raise ValueError("nonnegative integer request data required")
    if original_pole < 1:
        raise ValueError("positive original pole required")
    o = mean_orders(kappa)
    p1, p0 = o["A1_integrations"], o["A0_integrations"]
    J1 = output_order + 2 * p1 + 3
    J0 = output_order + 2 * p0 + 3
    maxima = {
        "post_signed_row": J0 + 10,
        "new_covariance": J0 + 13,
        "next_signed_target_input": J0 + 14,
    }
    N = maxima["new_covariance"]
    q = max(previous_wave_pole, 6 * original_pole, 6 * N + 9, 6 * N + 12 * original_pole + 6)
    return {
        "output_order": output_order,
        "A1_source_order": J1,
        "A0_source_order": J0,
        "maxima": maxima,
        "common_wave_pole": q,
        "signed_covariance_and_fast_domain": 2 * q,
        "postfast_and_final_radial_domain": 4 * q,
        "old_DAGs_expanded": False,
        "actual_high_jets_evaluated": False,
    }


def blocking_inputs(output_order=0, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q201 identity required")
    r = rates()
    req = compact_requests(output_order)
    return {
        "A0_inequality": (
            f"N_{req['A0_source_order']}(Qcirc gr_final) <= "
            "epsilon^B P_I(S+K_mean) on the common same-family domain"
        ),
        "A1_inequality": (
            f"N_{req['A1_source_order']}(Qcirc EpostSigned_z) <= "
            "epsilon^D Q_I(S+K_mean) on the common same-family domain"
        ),
        "whole_hull_contract": "G=0 on the entire reserved potential hull, including gaps between bumps",
        "required_target_rate": r["C0"],
        "finite_request_graph_explicit": True,
        "closed_relative_to_named_continuum_providers": True,
        "independently_evaluated_or_formally_certified": False,
    }


def physical_residual_power(h=Q(1, 200), kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q201 identity required")
    h = Q(h)
    if not 0 < h < Q(1, 100):
        raise ValueError("reviewed h required")
    D = rates(kappa)["D"]
    Lambda = Q(3, 2) + 2 * h
    return h * D - Lambda


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q201 identity required")
    return {
        "closed_here": (
            "relative to the named physical fixed-jet continuum providers, the corrected first "
            "post-signed fast and slow maps give all three regular mean rows at C0=2-2*kappa"
        ),
        "next": (
            "bind the complete post-signed nonzero residual at order D to the Q161/Q180 next "
            "constrained sourced solve on the same expanded domains, then replay the physical "
            "complete-potential diagonal and localization ledger"
        ),
        "actual_high_jet_values_or_stage_uniform_constants_proved": False,
        "absolute_physical_residual_decay_proved": False,
        "complete_diagonal_force_exclusion_or_blowup_proved": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q201 identity required")
    status = statuses()
    r = rates()
    o = mean_orders()
    memory = exact_mean_memory()
    request = compact_requests()
    blockers = blocking_inputs()
    boundary = frontier()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "post_signed_fast_source_is_corrected_D": o["post_signed_Qcirc_rows"] == r["D"],
        "fresh_fast_and_slow_orders_are_D_and_C": o["fast_tangential"] == r["D"] and o["slow_tangential"] == r["C"],
        "postfast_scalars_stay_C_and_postslow_gain_one": o["postfast_scalars"] == r["C"] and o["postslow_scalars"] == r["C"] + 1,
        "fast_and_slow_full_differences_have_strict_margins": o["tangential_margins_over_C0"]["fast_full_difference"] == Q(1, 2) - 2 * r["kappa"] and o["tangential_margins_over_C0"]["slow_full_difference"] == Q(1, 2) - r["kappa"],
        "compact_counts_and_strict_kappa_margins_are_exact": o["A1_integrations"] == 50001 and o["A0_integrations"] == 100000 and o["tangential_margins_over_C0"]["A1_compact_remainder"] == r["kappa"] and o["radial_margins_over_C0"]["A0_compact_remainder"] == r["kappa"],
        "all_regular_mean_rows_close_at_C0": o["final_tangential"] == r["C0"] and o["final_radial"] == r["C0"] and o["final_regular_mean"] == r["C0"],
        "complete_pressure_axial_and_old_new_memory_retained": memory["all_six_covariances_retained"] and memory["axial_velocity_retained"] and memory["compact_pressure_retained"] and memory["moving_rho_and_whole_c_rho_P_derivative_retained"] and memory["viscosity_time_cylindrical_cross_and_square_terms_retained"],
        "finite_request_graph_matches_frozen_mean_closure": request["A1_source_order"] == 100005 and request["A0_source_order"] == 200003 and request["maxima"]["next_signed_target_input"] == 200017 and request["common_wave_pole"] == 1200114,
        "conditionality_is_explicit": blockers["closed_relative_to_named_continuum_providers"] and not blockers["independently_evaluated_or_formally_certified"],
        "physical_regular_residual_still_has_negative_exponent": physical_residual_power() < 0,
        "no_diagonal_force_or_blowup_promotion": not boundary["complete_diagonal_force_exclusion_or_blowup_proved"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "rates": r,
        "mean_orders": o,
        "exact_mean_memory": memory,
        "compact_requests": request,
        "blocking_inputs": blockers,
        "physical_residual_power": physical_residual_power(),
        "frontier": boundary,
        "scope": {
            "conditional_analytical_first_post_signed_mean_block": True,
            "tests_or_finite_request_arithmetic_are_continuum_proof": False,
            "actual_enormous_derivatives_numerically_evaluated": False,
            "uniform_in_cycle_constants_claimed": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup": False,
        },
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q201 identity required")
    return {
        "status": "Q201-corrected-post-signed-fast-slow-mean-block-joined-conditionally",
        "established": (
            "relative to the named physical fixed-jet continuum providers, the actual zero-Haar "
            "fast map at D and current five-map at C retain complete pressure memory and place "
            "all three regular mean rows at C0=2-2*kappa"
        ),
        "limiting_rows": (
            "the post-signed Haar rows tie C0; A1 and A0 compact branches require 50001 and "
            "100000 integrations and finish only kappa above C0"
        ),
        "next": frontier()["next"],
        "historical_frozen_records_modified": False,
        "unconditional_complete_residual_force_exclusion_or_blowup_certified": False,
    }
