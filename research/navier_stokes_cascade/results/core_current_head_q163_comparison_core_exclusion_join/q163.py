"""Join Q162 to fixed-gauge comparison and the protected core, with legal limits."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
IDENTITY = "Q163-same-localized-field-comparison-core-exclusion-join"
PROVIDERS = (
    "core_current_head_q162_diagonal_endpoint_force_join",
    "core_current_head_q156_physical_pressure_comparison",
    "core_current_head_q157_physical_core_exclusion",
    "comparison_foundation",
    "core_same_field_positive_core_reconciliation",
    "core_current_head_q1_exact_map_join",
)


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


q156 = load("q163_q156", "core_current_head_q156_physical_pressure_comparison/q156.py")
q157 = load("q163_q157", "core_current_head_q157_physical_core_exclusion/q157.py")


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def comparison_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q163 identity required")
    pressure = q156.pressure_ledger()
    radius = q156.radius_ledger()
    comparison = q156.comparison_ledger()
    return {
        "pressure_difference": pressure["difference"],
        "pressure_recovery": pressure["recovery"],
        "common_force_cancels_once": pressure["force_cancellations"] == 1,
        "finite_energy_essential": pressure["finite_energy_essential"],
        "competitor_pressure_decay_assumed": pressure["competitor_pressure_decay_assumed"],
        "compact_mixed_flux": radius["compact_mixed_flux"],
        "worst_radius_decay": radius["worst_decay"],
        "cutoff_power": radius["cutoff"],
        "comparison": comparison,
    }


def core_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q163 identity required")
    parameters = q157.fixed_parameters()
    binding = q157.same_physical_binding()
    threshold = q157.combined_threshold(1, 2, Q(1, 8), Q(1, 16), Q(1, 4))
    return {
        "parameters": parameters,
        "binding": binding,
        "lower_bound": q157.lower_bound_statement(),
        "threshold_example": threshold,
        "threshold_positive": threshold["tau_star"] > 0,
        "path_inside_annular_support_hole": binding["path_inside_hole"],
        "localization_plateau_preserved": "c=1 and grad c=0" in binding["localization"],
    }


def limit_order(tau=None, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q163 identity required")
    threshold = core_contract()["threshold_example"]["tau_star"]
    tau = threshold / 2 if tau is None else Q(tau)
    if not 0 < tau <= threshold:
        raise ValueError("tau must lie on the protected path")
    t = 1 - tau
    T = 1 - tau / 2
    return {
        "tau": tau,
        "t": t,
        "T_tau": T,
        "strict_time_order": 0 < t < T < 1,
        "first": "fix T_tau<1 and recover the pressure gradient",
        "second": "send the comparison radius R to infinity at this fixed T_tau",
        "third": "deduce v=u at t=1-tau, then let tau decrease to zero",
        "uniform_comparison_constant_as_T_to_one_required": False,
        "reason": (
            "the R-limit gives exact equality separately for each tau; no bound on "
            "C_T or L_T uniform as T approaches one is used"
        ),
    }


def contradiction_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q163 identity required")
    order = limit_order()
    A = core_contract()["parameters"]["A"]
    return {
        "comparison_applies_at_path_time": order["strict_time_order"],
        "path": "x_tau=(sqrt(2 X_in tau),0,0) tends to zero",
        "growth": "|u_theta(x_tau,1-tau)| >= (e0/2) tau^(-A)",
        "positive_exponent": A > 0,
        "global_smooth_competitor_property": (
            "bounded velocity on one compact space-time neighborhood of (0,1)"
        ),
        "conditional_conclusion": (
            "no globally smooth same-force zero-datum finite-energy competitor exists"
        ),
        "future_velocity_or_pressure_constructed": False,
        "unconditional_blowup_claimed": False,
    }


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q163 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    expected = {
        PROVIDERS[0]: "Q162-physical-diagonal-endpoint-localized-force-chain-joined",
        PROVIDERS[1]: "Q156-physical-fixed-gauge-pressure-and-expanding-ball-comparison-conditionally",
        PROVIDERS[2]: "Q157-same-physical-core-path-and-competitor-exclusion-conditionally",
        PROVIDERS[3]: "current-source-compact-pressure-comparison-relative-to-retained-core-and-finite-analytic-lemmas",
        PROVIDERS[4]: "same-complete-diagonal-positive-core-bound-reconciled",
        PROVIDERS[5]: "current-head-Q1-exact-core-roots-and-leading-pressure-bound-on-one-profile",
    }
    for name, status in expected.items():
        if data[name].get("status") != status:
            raise ValueError("Q163 provider status mismatch: " + name)
    q162 = data[PROVIDERS[0]]["symbolic"]
    q156 = data[PROVIDERS[1]]["symbolic"]
    q157 = data[PROVIDERS[2]]["symbolic"]
    comparison_flags = data[PROVIDERS[3]]["flags"]
    core = data[PROVIDERS[4]]["review"]["analytic_conclusions"]
    q1 = data[PROVIDERS[5]]["review"]
    required = (
        q162["P6_complete_diagonal_endpoint_localized_force_contract_discharged"],
        q156["preterminal_comparison_certified_conditionally"],
        q157["same_physical_core_and_conditional_exclusion_certified"],
        comparison_flags["pressure_gradient_recovered_without_assumed_competitor_pressure_decay"],
        comparison_flags["full_expanding_ball_energy_remainder_improved_from_four_thirds_to_two"],
        core["same_complete_diagonal_has_positive_core_lower_bound"],
        core["localized_field_preserves_the_bound_on_its_unit_plateau"],
        q1["checks"]["all_ten_Q1_provider_statuses_match"],
        q1["checks"]["nonzero_pressure_homogeneous_mode_is_endpoint_unbounded"],
    )
    if not all(required):
        raise ValueError("Q163 comparison/core provider contract missing")
    return {
        "same_Q162_field_force_pressure_gauge_and_localization": True,
        "fixed_T_then_radius_then_core_limit_order": True,
        "P7_conditional_consequence_closed": True,
        "historical_conditional_records_preserved": True,
    }


def remaining_foundation_gate(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q163 identity required")
    comparison = evidence("comparison_foundation")
    q1 = evidence("core_current_head_q1_exact_map_join")
    return {
        "comparison_foundation_all_dependencies_discharged": comparison["flags"][
            "all_retained_analytic_dependencies_discharged"
        ],
        "Q1_every_provider_theorem_independently_reproved": q1["review"]["scope"][
            "every_provider_theorem_independently_reproved_here"
        ],
        "first_audit": (
            "independently rederive the selected all-index core-map ball invariance and "
            "strict contraction, recover the original equations, and identify its pressure "
            "with the completed leading backward pressure on the Q1-Q162 object"
        ),
        "required_inequalities": (
            "T(B_r) subset B_r and ||T(Phi)-T(Psi)||_X <= lambda "
            "||Phi-Psi||_X with one explicit lambda<1, followed by "
            "Pi_axis+C_p,infinity=0 in the same fixed gauge"
        ),
        "dependency_hashes_or_passing_tests_suffice": False,
        "unconditional_theorem_gate_closed": False,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q163 identity required")
    return (
        {"id": "P1", "status": "closed-relative-to-named-providers"},
        {"id": "P2", "status": "closed-relative-to-named-providers"},
        {"id": "P3", "status": "closed-relative-to-named-providers"},
        {"id": "P4", "status": "closed-relative-to-named-providers"},
        {"id": "P5", "status": "closed-relative-to-named-providers"},
        {"id": "P6", "status": "closed-relative-to-named-providers"},
        {"id": "P7", "status": "closed-as-conditional-consequence-by-Q163"},
        {"id": "U1", "status": "next-independent-foundation-audit"},
    )


def record():
    return {
        "status": "Q163-same-object-comparison-and-core-exclusion-joined-conditionally",
        "correction": (
            "the terminal comparison needs no constants uniform as T approaches one: for "
            "each protected tau choose T_tau=1-tau/2, take R to infinity at fixed T_tau, "
            "and only then send tau to zero"
        ),
        "established": (
            "Q162's exact localized field, complete fixed-gauge pressure and compact force "
            "satisfy the Q156 comparison hypotheses; Q157's path lies on the same unit "
            "plateau and forces conditional exclusion of a globally smooth same-force "
            "zero-datum finite-energy competitor"
        ),
        "scope": (
            "same-object implication relative to the named analytical continuum providers; "
            "no future field is constructed and dependency hashes or tests are not a proof "
            "of the retained provider theorems"
        ),
        "next": remaining_foundation_gate()["first_audit"],
        "P7_comparison_and_core_exclusion_replayed_conditionally": True,
        "unconditional_admissible_force_breakdown_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
    }
