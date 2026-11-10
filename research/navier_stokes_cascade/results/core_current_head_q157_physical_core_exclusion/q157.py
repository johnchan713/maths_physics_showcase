"""Bind the positive core path to the Q153--Q156 physical field."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q157-same-physical-diagonal-core-path-and-conditional-exclusion"
PROVIDERS = (
    "core_same_field_positive_core_reconciliation",
    "core_current_head_q149_primary_pulse_join",
    "core_certified_complete_potential_diagonal_join",
    "core_current_head_q153_physical_source_tail_diagonal_join",
    "core_current_head_q155_physical_force_energy_join",
    "core_current_head_q156_physical_pressure_comparison",
)


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def fixed_parameters(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q157 identity required")
    h = Q(1, 200)
    return {
        "h": h,
        "A": Q(1, 2) + h,
        "X_a": Q(4, 101),
        "X_in": Q(1, 101),
        "two_h_reciprocal": 100,
    }


def path_coordinates(tau, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q157 identity required")
    tau = Q(tau)
    p = fixed_parameters()
    if tau <= 0:
        raise ValueError("Positive tau required")
    return {
        "q": tau,
        "eta": Q(0),
        "X": p["X_in"],
        "r_squared": 2 * p["X_in"] * tau,
        "z": Q(0),
        "t": 1 - tau,
    }


def core_threshold(e0, background_constant, q0, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q157 identity required")
    e0, background_constant, q0 = map(Q, (e0, background_constant, q0))
    if min(e0, background_constant, q0) <= 0:
        raise ValueError("Positive core data required")
    ratio = min(Q(1), e0 / (2 * background_constant))
    return min(q0, ratio ** fixed_parameters()["two_h_reciprocal"])


def combined_threshold(e0, background_constant, q0, tau0,
                       radius_squared, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q157 identity required")
    tau0, radius_squared = map(Q, (tau0, radius_squared))
    if min(tau0, radius_squared) <= 0:
        raise ValueError("Positive localization data required")
    p = fixed_parameters()
    q_core = core_threshold(e0, background_constant, q0)
    plateau = min(tau0 / 2, radius_squared / (8 * p["X_in"]))
    return {
        "q_core": q_core,
        "plateau": plateau,
        "tau_star": min(q_core, plateau),
    }


def same_physical_binding(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q157 identity required")
    p = fixed_parameters()
    return {
        "background": (
            "Q149 uses the same finalized leading profile and rule-defined "
            "Borel background as the positive-core provider"
        ),
        "diagonal": (
            "Q153 applies the certified complete-potential diagonal to the "
            "inherited physical family with unchanged closed annular supports"
        ),
        "support_hole": "every stage increment vanishes for X<X_a",
        "path_inside_hole": p["X_in"] < p["X_a"],
        "background_remainder": (
            "q^A u_theta=e0+R_bg; |R_bg|<=C_bg q^(2h)"
        ),
        "localization": (
            "Q155 has c=1 and grad c=0 on the late core, hence u=u0 there"
        ),
        "canonical_physical_primary_identification_used": False,
        "pressure_gauge_changed": False,
    }


def lower_bound_statement(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q157 identity required")
    p = fixed_parameters()
    return {
        "leading_coefficient": (
            "e0=E0(X_in,0)=sqrt(2 X_in)*varphi_star*Phi/C>0"
        ),
        "path": "x_tau=(sqrt(2 X_in tau),0,0), t=1-tau",
        "exponent": p["A"],
        "bound": (
            "|u_theta(x_tau,1-tau)|>=(e0/2) tau^(-A) "
            "for 0<tau<=tau_star"
        ),
        "tau_star": (
            "min(q0,(e0/(2 C_bg))^100,tau0/2,r_L^2/(8 X_in))"
        ),
        "conditional": True,
    }


def exclusion_statement(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q157 identity required")
    return {
        "comparison_input": (
            "Q156: v=u on R^3 x [0,T] for every T<1 for each locally smooth "
            "same-force zero-datum competitor with fixed-interval L2 bound"
        ),
        "terminal_input": (
            "Q157 path tends to (0,1) and |u_theta| tends to infinity"
        ),
        "conclusion": (
            "no such competitor can also be smooth on a neighborhood of (0,1)"
        ),
        "global_smooth_finite_energy_competitor_excluded_conditionally": True,
        "unconditional_blowup_claimed": False,
        "future_velocity_or_pressure_constructed": False,
    }


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q157 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    expected = {
        PROVIDERS[0]: "same-complete-diagonal-positive-core-bound-reconciled",
        PROVIDERS[1]: "Q148-schedule-frontier-corrected-and-uniform-primary-pulse-family-joined",
        PROVIDERS[2]: "certified-domain-conditional-complete-potential-diagonal-join",
        PROVIDERS[3]: "Q153-physical-source-tail-and-complete-potential-diagonal-bound-conditionally",
        PROVIDERS[4]: "Q155-physical-localization-compact-force-and-energy-bound-conditionally",
        PROVIDERS[5]: "Q156-physical-fixed-gauge-pressure-and-expanding-ball-comparison-conditionally",
    }
    if any(data[name].get("status") != status for name, status in expected.items()):
        raise ValueError("Q157 provider status mismatch")
    core = data[PROVIDERS[0]]["review"]
    conclusions = core["analytic_conclusions"]
    needed = (
        "leading_cartesian_swirl_coefficient_is_E0_at_z_zero",
        "rule_defined_positive_order_background_is_o_one_after_qA_normalization",
        "all_annular_diagonal_increments_vanish_exactly_on_the_path",
        "same_complete_diagonal_has_positive_core_lower_bound",
        "localized_field_preserves_the_bound_on_its_unit_plateau",
    )
    if not all(conclusions.get(key) for key in needed):
        raise ValueError("Positive-core provider contract missing")
    q149 = data[PROVIDERS[1]]["symbolic"]
    if "same finalized leading/background profile" not in q149["established"]:
        raise ValueError("Q149 same background binding missing")
    if not data[PROVIDERS[3]]["symbolic"][
            "complete_potential_diagonal_applicable_conditionally"]:
        raise ValueError("Q153 physical diagonal binding missing")
    if not data[PROVIDERS[5]]["symbolic"][
            "preterminal_comparison_certified_conditionally"]:
        raise ValueError("Q156 comparison missing")
    return {
        "same_finalized_background": True,
        "physical_diagonal_has_unchanged_annular_support_hole": True,
        "Q155_localization_plateau": True,
        "Q156_preterminal_comparison": True,
        "all_provider_conclusions_remain_conditional": True,
    }


def record():
    return {
        "status": "Q157-same-physical-core-path-and-competitor-exclusion-conditionally",
        "correction": (
            "Q156 correctly refused to import the older path without a family binding; "
            "Q149's exact shared Borel background and Q153's unchanged annular support "
            "hole now provide that binding without identifying canonical and physical pulses"
        ),
        "established": (
            "the exact Q153-Q155 physical field satisfies the protected-core lower bound "
            "with A=101/200 and an explicit positive combined threshold; with Q156 this "
            "conditionally excludes a same-force zero-datum smooth finite-energy competitor"
        ),
        "scope": (
            "analytical composition of named conditional continuum providers; exact tests "
            "check arithmetic and provenance but do not prove those providers"
        ),
        "next": (
            "audit the current physical Q147-Q157 dependency and quantifier closure, list "
            "every still-conditional continuum primitive, and prove or precisely obstruct "
            "the first unresolved one before any unconditional theorem claim"
        ),
        "same_physical_core_and_conditional_exclusion_certified": True,
        "unconditional_admissible_force_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
    }
