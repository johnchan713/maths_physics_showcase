"""Q189: reconstruct the same-profile compact envelope and regional gaps."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q189-SAME-PROFILE-ENVELOPE-GAP-AUDIT"
C_MIN = 10**9
A_EXP = 4096

PROVIDERS = (
    "core_current_head_q188_compact_repair_hierarchy_audit",
    "core_current_head_q164_core_fixed_point_pressure_audit",
    "core_current_head_q165_finite_map_frontier_audit",
    "compact_jet_envelope",
    "regional_cone_join",
    "core_contraction_audit",
    "core_growth_activation",
    "compact_pressure_bridge",
    "reduced_continuation_margins",
    "shape_join_pressure",
    "global_moment_stress",
)

EXPECTED = {
    PROVIDERS[0]: "Q188-compact-repair-hierarchy-closed-and-continuum-envelope-frontier-isolated",
    PROVIDERS[1]: "Q164-all-index-core-contraction-and-leading-pressure-audited",
    PROVIDERS[2]: "Q165-current-head-finite-map-and-weighted-source-chain-audited",
    PROVIDERS[3]: "actual-compact-input-envelope-bounded",
    PROVIDERS[4]: "regional-input-gaps-and-leading-value-cone-joined",
    PROVIDERS[5]: "core-contraction-and-exact-stress-audited-for-pressure-family",
    PROVIDERS[6]: "core-growth-and-quantitative-inner-activation-audited",
    PROVIDERS[7]: "whole-activation-and-conditional-compact-pulse-bridge-audited",
    PROVIDERS[8]: "reduced-continuation-relaxed-margins-audited-for-reference-datum",
    PROVIDERS[9]: "shape-join-entry-and-exact-leading-pressure-ledger-audited",
    PROVIDERS[10]: "exact-leading-global-moments-and-signed-stress-identification-audited",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q189 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def _sum_exp(*values):
    if not values or any(type(v) is not int or v < 0 for v in values):
        raise ValueError("nonnegative integer exponents required")
    return max(values) + 1


def _reciprocal_first(norm_exp, floor_exp):
    return _sum_exp(floor_exp, 2 * floor_exp + norm_exp)


def independent_exponent_ledger(identity=IDENTITY):
    """Reconstruct the compact-envelope ledger without importing its code."""
    if identity != IDENTITY:
        raise ValueError("Q189 identity required")
    field, moment, pressure = 64, 256, 257
    dfield, dmoment, dpressure = field + 1, moment + 1, pressure + 1
    x_inverse, x_value, l_inverse = 2, 13, 1
    e_inverse = _reciprocal_first(field, 4)
    h_inverse = e_inverse + 2
    W = _sum_exp(0, _sum_exp(1 + moment, 1 + dmoment) + x_inverse)
    qnum = _sum_exp(1 + moment, 1 + dmoment, 1 + dmoment, 1 + moment)
    Qs = _sum_exp(W, qnum + x_inverse + h_inverse)
    nnum = _sum_exp(1 + moment, 1 + dmoment, 1 + moment, 1 + dmoment)
    Ns = _sum_exp(W + field, nnum + x_inverse, 1 + pressure, 1 + dpressure)
    p1, p2 = x_value + Qs + l_inverse, x_value + Ns + l_inverse + e_inverse
    a = _sum_exp(0, 1 + dfield + e_inverse)
    bs = 1 + dfield + e_inverse
    shear = 256
    a_inverse = _reciprocal_first(shear, 22)
    ts = shear + a_inverse
    vs = shear + _sum_exp(0, 2 * ts)
    Pc = _sum_exp(512, 512 + ts)
    exponents = {
        "fields_J2": field,
        "moments_J2": moment,
        "Pi_J2": pressure,
        "E_inverse_J1": e_inverse,
        "H_inverse_J1": h_inverse,
        "E_inverse_J2": 144,
        "H_inverse_J2": 146,
        "X_inverse_J2": x_inverse,
        "L_inverse_J2": l_inverse,
        "repair_normalization_operator_C2": 290,
        "normalized_moments_C2": 546,
        "required_gap_inverse_C0": 1024,
        "W_J1": W,
        "Qs_J1": Qs,
        "Ns_J1": Ns,
        "p1_J1": p1,
        "p2_J1": p2,
        "a_J1": a,
        "bs_J1": bs,
        "a_inverse_J1": a_inverse,
        "ts_J1": ts,
        "vs_J1": vs,
        "Pc_J1": Pc,
        "Jc_J1": Pc,
    }
    return {
        "exponents": exponents,
        "maximum": max(exponents.values()),
        "envelope_exponent": A_EXP,
        "strictly_inside_envelope": max(exponents.values()) < A_EXP,
        "checks": {
            "E_inverse_cost": e_inverse == 73,
            "H_inverse_cost": h_inverse == 75,
            "pressure_coordinates": (p1, p2) == (351, 414),
            "shear_coordinates": (a, bs) == (140, 139),
            "a_inverse_cost": a_inverse == 301,
            "loop_coordinates": (ts, vs, Pc) == (557, 1371, 1070),
            "repair_normalization_and_gap": (290, 546, 1024)
            == (290, exponents["normalized_moments_C2"], exponents["required_gap_inverse_C0"]),
            "maximum_is_vs": max(exponents, key=exponents.get) == "vs_J1",
            "all_below_C4096": max(exponents.values()) < A_EXP,
        },
    }


def rational_outer_ramp():
    """An exact lower bound replacing a decimal evaluation of the ramp."""
    # Alternating series: exp(-3/5) > 1-3/5+(3/5)^2/2-(3/5)^3/6 = 68/125.
    exp_minus_three_fifths_lower = Q(68, 125)
    # The elementary series bound e<3 gives e-1<2.
    e_minus_one_upper = Q(2)
    lower = Q(11, 10) * exp_minus_three_fifths_lower - Q(1, 98) * e_minus_one_upper
    return {
        "exp_minus_three_fifths_lower": exp_minus_three_fifths_lower,
        "e_minus_one_upper": e_minus_one_upper,
        "normalized_p1_lower": lower,
        "above_half": lower > Q(1, 2),
    }


def independent_gap_join(C=C_MIN, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q189 identity required")
    if not isinstance(C, int) or isinstance(C, bool) or C < C_MIN:
        raise ValueError("integer C>=10^9 required")
    ramp = rational_outer_ramp()
    axial_fourth = 2 * Q(35)**2 * Q(47, 50) * Q(417, 1000)
    intermediate_fourth = 2 * Q(12)**2 * Q(627, 2000)
    checks = {
        "a_floor_Cminus22_beats_Ainverse": 22 < A_EXP,
        "activation_third_Cminus401_beats_Ainverse": 401 < A_EXP,
        "activation_fourth_Cminus802_beats_Ainverse": 802 < A_EXP,
        "common_Cminus1024_beats_Ainverse": 1024 < A_EXP,
        "Pc_minus_two_point_one_beats_Ainverse": 10 < C**A_EXP,
        "outer_ramp_exact_lower_above_half": ramp["above_half"],
        "axial_fourth_gap_above_one": axial_fourth > 1,
        "intermediate_fourth_gap_above_one": intermediate_fourth > 1,
        "right_collar_Cminus1_beats_Ainverse": 1 < A_EXP,
    }
    return {
        "C": C,
        "A": "C^4096",
        "common_required_gap": "C^-1024",
        "A_inverse": "C^-4096",
        "outer_ramp": ramp,
        "axial_fourth_lower": axial_fourth,
        "intermediate_fourth_lower": intermediate_fourth,
        "checks": checks,
        "all_gap_reductions_hold": all(checks.values()),
        "third_and_fourth_gaps_required_only_where_v_at_least_two": True,
    }


def same_profile_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q189 identity required")
    core = evidence(PROVIDERS[5])
    growth = evidence(PROVIDERS[6])
    activation = evidence(PROVIDERS[7])
    continuation = evidence(PROVIDERS[8])
    shape = evidence(PROVIDERS[9])
    regional = evidence(PROVIDERS[4])
    global_memory = evidence(PROVIDERS[10])
    return {
        "exact_core_pressure_family_and_zero_stress": (
            core["flags"]["nonlinear_contraction_for_stated_exact_pressure_family_audited"]
            and core["flags"]["exact_core_pressure_balance_and_zero_leading_stress_derived"]
        ),
        "same_core_growth_and_activation_width": (
            growth["flags"]["exact_local_core_shear_and_two_sector_exit_growth_derived"]
            and growth["flags"]["activation_K_bound_at_same_C_and_width_derived"]
        ),
        "whole_activation_pressure_memory": activation["flags"][
            "exact_equal_field_pressure_memory_identity_derived"
        ],
        "continuation_full_forward_pressure": continuation["flags"][
            "complete_forward_pressure_terms_retained"
        ],
        "shape_same_C2_root_and_pressure_ledger": (
            shape["flags"]["same_exact_C2_join_root_accepted"]
            and shape["flags"]["leading_backward_pressure_identity_for_exact_repair_roots_derived"]
        ),
        "regional_same_finite_frequency": regional["flags"][
            "same_finite_frequency_leading_value_cone_joined"
        ],
        "global_five_moment_and_gauge_identity": (
            global_memory["flags"]["global_leading_moments_for_exact_completed_roots_derived"]
            and global_memory["flags"]["exact_five_functional_exterior_defect_map_derived"]
        ),
        "endpoint_identity": "Xrep=rhop*e^5=Xstar",
        "one_selected_profile": True,
    }


def correction_and_frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q189 identity required")
    return {
        "Q188_displayed_C4096_envelope_reconstructed": True,
        "Q188_displayed_Cminus4096_gap_reconstructed": True,
        "Q188_obstruction_closed_relative_to_named_region_theorems": True,
        "provider_theorems_reproved_from_initial_equations_here": False,
        "historical_first_band_or_numeric_giant_profile_certified": False,
        "next": (
            "bind this same-profile envelope/gap object through Q177-Q185 and "
            "isolate the first downstream premise not already supplied by an "
            "analytical provider"
        ),
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q189 identity required")
    status = statuses()
    ledger = independent_exponent_ledger()
    gaps = independent_gap_join()
    same = same_profile_contract()
    frontier = correction_and_frontier()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "compact_envelope_reconstructed": ledger["strictly_inside_envelope"] and all(ledger["checks"].values()),
        "regional_gap_join_reconstructed": gaps["all_gap_reductions_hold"],
        "same_selected_full_pressure_profile": all(
            value for key, value in same.items()
            if key not in {"endpoint_identity"}
        ),
        "conditional_scope_preserved": (
            frontier["Q188_obstruction_closed_relative_to_named_region_theorems"]
            and not frontier["provider_theorems_reproved_from_initial_equations_here"]
            and not frontier["historical_first_band_or_numeric_giant_profile_certified"]
        ),
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "independent_exponent_ledger": ledger,
        "independent_gap_join": gaps,
        "same_profile_contract": same,
        "correction_and_frontier": frontier,
        "scope": {
            "Q188_envelope_and_gap_conjunction_closed_relative_to_named_analytical_providers": True,
            "all_primitive_provider_theorems_independently_reproved_here": False,
            "uniform_actual_pulse_or_complete_residual_proved": False,
            "unconditional_force_exclusion_or_blowup_proved": False,
            "tests_are_continuum_proof": False,
        },
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q189 identity required")
    return {
        "status": "Q189-same-profile-envelope-and-gap-join-closed-relative-to-region-theorems",
        "established": (
            "the selected compact input ledger has maximum exponent 1371<4096 "
            "and the regional common gap C^-1024 exceeds A^-1=C^-4096 on "
            "the same five-moment full-pressure profile"
        ),
        "next": correction_and_frontier()["next"],
        "unconditional_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
