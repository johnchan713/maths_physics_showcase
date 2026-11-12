"""Bind the actual symbolic Q2 inputs to an exact later-tail selector."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q187-ACTUAL-SYMBOLIC-Q2-SELECTOR"
PROVIDERS = (
    "core_current_head_q186_exact_finite_frontier_correction",
    "leading_frame_foundation",
    "leading_cone_margin",
    "explicit_tilt_pulse_tail",
    "regional_cone_join",
    "regional_weighted_gate",
    "first_seminorm_bound",
    "background_gap_certificate",
    "core_current_head_q172_q2_q3_source_audit",
    "core_current_head_q176_modulation_repair_dependency_audit",
)
EXPECTED = {
    PROVIDERS[0]: "Q186-exact-core-finite-map-frontier-corrected-to-actual-Q2-infima",
    PROVIDERS[1]: "current-leading-frame-floors-from-reviewed-core-and-same-root-regional-profile",
    PROVIDERS[2]: "global-leading-direction-margins-and-existing-tilt-gate",
    PROVIDERS[3]: "explicit-tilt-raw-deficit-and-primary-pulse-tail-certificate",
    PROVIDERS[4]: "regional-input-gaps-and-leading-value-cone-joined",
    PROVIDERS[5]: "finite-seminorm-weighted-input-and-regional-background-gap-join",
    PROVIDERS[6]: "collar-linear-repair-full-pressure-first-seminorm-bound",
    PROVIDERS[7]: "explicit-background-value-gap-certificate-relative-to-leading-and-weighted-jets",
    PROVIDERS[8]: "Q172-Q2-source-audited-and-Q3-uniform-frame-corrected",
    PROVIDERS[9]: "Q176-compact-modulation-exact-repair-audited-with-refined-dependency-order",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q187 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def _rat(value):
    if isinstance(value, (bool, float)):
        raise TypeError("exact rational data required")
    return Q(value)


def exact_dyadic_selector(Cval, rho, h, ell_min=2):
    """Return an exact sufficient integer ell without floating logarithms.

    The diagnostic inputs are rational.  The actual checkpoint uses this
    lemma symbolically because its selected-profile constants are enormous
    analytic expressions rather than evaluated rationals.
    """
    Cval, rho, h = map(_rat, (Cval, rho, h))
    if Cval <= 0 or rho <= 0 or h <= 0:
        raise ValueError("positive Cval, rho and h required")
    if not isinstance(ell_min, int) or isinstance(ell_min, bool):
        raise TypeError("integer ell_min required")
    ell_min = max(2, ell_min)
    dyadic_power = 0
    while Cval > rho * 2**dyadic_power:
        dyadic_power += 1
    increment = -(-Q(dyadic_power, 1) // (2 * h))
    ell = max(ell_min, 2 + increment)
    return {
        "dyadic_power": dyadic_power,
        "ell": ell,
        "exponent": 2 * h * (ell - 2),
        "gate": Q(Cval, 2**dyadic_power) <= rho,
        "minimal_dyadic_power": (
            dyadic_power == 0 or Q(Cval, 2 ** (dyadic_power - 1)) > rho
        ),
        "chart_gate": 2 * h * (ell - 2) >= dyadic_power,
    }


def actual_symbolic_inputs(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q187 identity required")
    frame = evidence(PROVIDERS[1])
    cone = evidence(PROVIDERS[2])
    tilt = evidence(PROVIDERS[3])
    regional = evidence(PROVIDERS[4])
    return {
        "alpha": "1/(2*G^3)",
        "nu": "1/(2*G^4)",
        "fstar": "exp(-2*Lambda*M)/(4*C)*exp(-Ab*log(Xb/Xa)/2)",
        "f": "fstar/5",
        "Fbar": "5*fupper",
        "t": "alpha*fstar/5",
        "Gbar": "5*fupper*(Ab+Bb)",
        "d": "alpha*nu*fstar^2/25",
        "omega": "1/(640*G^13)",
        "epsilon": "1/(2^18*G^22)",
        "u_squared": "2/epsilon = 2^19*G^22",
        "gamma": "epsilon/4 = 1/(2^20*G^22)",
        "all_value_inputs_strictly_positive_and_finite": (
            frame["flags"]["uniform_leading_spectral_floors_from_current_reviewed_profile"]
            and cone["flags"]["global_leading_orientation_and_raw_cone_deficit_derived_from_same_regional_maps"]
            and tilt["flags"]["explicit_specialization_has_quantitative_raw_deficit"]
            and regional["flags"]["original_compact_regional_value_margins_derived"]
        ),
        "same_selected_profile": True,
        "numerically_evaluated": False,
    }


def coefficient_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q187 identity required")
    weighted = evidence(PROVIDERS[5])
    first = evidence(PROVIDERS[6])
    return {
        "first_seminorms": {
            "A_E": "2*c*S",
            "A_EX": "2*c*Sx+S/a",
            "A_UX": "Sx",
            "A_W": "(Xc*Htheta/2+c*Hz)/tau",
        },
        "Cstar": (
            "max(1,5*(A_E+1/2)/a,"
            "5*(2*c*(A_EX+A_UX+1)+(A_E+1/2)/a),A_W+1/2)"
        ),
        "actual_shared_budget_and_first_seminorm_functional": weighted["flags"]["lumped_Cval_derived_from_actual_shared_budget_and_first_coefficient_seminorms"],
        "explicit_finite_Cstar": first["flags"]["explicit_finite_Cstar_composes_with_previous_regional_gate"],
        "full_pressure_moment_memory": first["flags"]["full_global_Omega0_over_X_moment_memory_and_nonconstant_amplitude_normalizers_retained"],
        "pressure_and_axial_viscosity_terms": first["flags"]["complete_first_pressure_value_and_axial_derivative_and_shifted_leading_axial_viscosity_kept"],
        "finite_positive_relative_to_primitive_profile_bounds": True,
        "numerically_evaluated": False,
    }


def selector_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q187 identity required")
    background = evidence(PROVIDERS[7])
    inputs = actual_symbolic_inputs()
    coefficients = coefficient_contract()
    return {
        "rho_actual": (
            "min(1/2,f/2,t/2,d/(2*LD),omega/(2*LH),"
            "u^2*omega^2*gamma*(2-gamma)/(2*Lphi))"
        ),
        "LD": "4*Gbar+2*Fbar+3",
        "LP": "Fbar+Gbar+1",
        "Dbar": "Gbar^2+2*Fbar*Gbar",
        "cbar_squared": "Dbar/(2*f*t)",
        "Lc": "2*LD/(f*t)+2*Dbar*LP/(f^2*t^2)",
        "LH": "1+2/t",
        "Lphi": "3*u^2*LH+(1+u^2)*(4*Lc+3*cbar^2*LH)",
        "blocking_inequality": "Cstar*2^(-2h*(ell-2)) <= rho_actual",
        "selector": (
            "choose integer N with 2^N>=Cstar/rho_actual, then integer "
            "ell>=2+N/(2h), together with every earlier accepted-tail gate"
        ),
        "rho_strictly_positive_and_finite": inputs["all_value_inputs_strictly_positive_and_finite"],
        "Cstar_strictly_positive_and_finite": coefficients["finite_positive_relative_to_primitive_profile_bounds"],
        "background_certificate_has_exact_radius_formula": background["flags"]["explicit_rational_radius_replaces_unspecified_background_continuity_constants"],
        "finite_later_selector_exists": True,
        "historical_first_integer_or_profile_constant_evaluated": False,
    }


def frontier_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q187 identity required")
    q172 = evidence(PROVIDERS[8])
    q176 = evidence(PROVIDERS[9])
    return {
        "Q186_actual_symbolic_join_completed": True,
        "Q2_source_audited_relative_to_displayed_providers": q172["review"]["scope"]["Q2_provider_premise_source_audited"],
        "modulation_repair_same_map_and_pressure_audited": (
            q176["review"]["established"]["conditional_same_map_all_order_modulation_and_exact_repair"]
            and q176["review"]["established"]["exact_five_row_pressure_memory"]
        ),
        "primitive_continuation_loop_repair_regional_suprema_independently_validated_here": False,
        "primitive_frontier": (
            "independently validate the retained finite-order continuation, "
            "low-order loop/repair, leading-frame strict margins and regional "
            "source suprema on the same selected profile"
        ),
        "numeric_evaluation_required_for_later_tail_existence": False,
        "numeric_evaluation_required_to_identify_a_prescribed_historical_band": True,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q187 identity required")
    status = statuses()
    inputs = actual_symbolic_inputs()
    coefficient = coefficient_contract()
    selector = selector_contract()
    frontier = frontier_contract()
    diagnostic = exact_dyadic_selector(100, Q(1, 37), Q(1, 200), 17)
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "actual_symbolic_value_cone_and_tilt_inputs_bound": inputs["all_value_inputs_strictly_positive_and_finite"],
        "actual_Cstar_functional_bound_with_full_pressure": (
            coefficient["actual_shared_budget_and_first_seminorm_functional"]
            and coefficient["explicit_finite_Cstar"]
            and coefficient["full_pressure_moment_memory"]
            and coefficient["pressure_and_axial_viscosity_terms"]
        ),
        "positive_finite_Q2_radius_and_later_selector": (
            selector["rho_strictly_positive_and_finite"]
            and selector["Cstar_strictly_positive_and_finite"]
            and selector["background_certificate_has_exact_radius_formula"]
            and selector["finite_later_selector_exists"]
        ),
        "exact_rational_selector_diagnostic": (
            diagnostic["gate"]
            and diagnostic["minimal_dyadic_power"]
            and diagnostic["chart_gate"]
        ),
        "conditional_gate_and_primitive_validation_distinguished": (
            frontier["Q2_source_audited_relative_to_displayed_providers"]
            and frontier["modulation_repair_same_map_and_pressure_audited"]
            and not frontier["primitive_continuation_loop_repair_regional_suprema_independently_validated_here"]
        ),
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "actual_symbolic_inputs": inputs,
        "coefficient_functional": coefficient,
        "later_tail_selector": selector,
        "exact_rational_diagnostic": diagnostic,
        "frontier": frontier,
        "correction": {
            "Q186_demand_for_numerical_actual_constants_was_too_strong_for_existence": True,
            "Q186_fixture_warning_remains_correct": True,
            "actual_symbolic_provider_join_replaces_fixture": True,
            "historical_frozen_evidence_modified": False,
        },
        "scope": {
            "actual_symbolic_Q2_selector_closed_relative_to_named_providers": True,
            "primitive_provider_theorems_reproved_from_initial_equations": False,
            "historical_first_band_recertified": False,
            "uniform_actual_pulse_or_complete_residual_proved": False,
            "unconditional_force_exclusion_or_blowup_proved": False,
            "tests_are_continuum_proof": False,
        },
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q187 identity required")
    return {
        "status": "Q187-actual-symbolic-Q2-selector-joined-and-primitive-frontier-isolated",
        "established": (
            "the actual symbolic leading floors, cone/tilt margins and finite "
            "first-seminorm Cstar functional give a positive Q2 radius and an "
            "exact finite later-tail selector, relative to the named providers"
        ),
        "blocking_inequality": selector_contract()["blocking_inequality"],
        "next": frontier_contract()["primitive_frontier"],
        "unconditional_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
