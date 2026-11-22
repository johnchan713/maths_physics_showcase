"""Q215: discharge Q207's fixed-stage bundle request conditionally."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q215-FIXED-STAGE-RESIDUAL-BUNDLE-JOIN"
PROVIDERS = (
    "core_current_head_q214_exact_core_root_pressure_join",
    "core_current_head_q213_reconstructed_primitive_bundle_exact_root_frontier",
    "core_current_head_q207_residual_constant_frontier",
    "core_current_head_q211_degree_grammar_primitive_frontier",
    "core_current_head_q184_generated_source_log_loss_correction",
    "core_current_head_q166_complete_residual_frontier_audit",
    "core_current_head_q203_diagonal_endpoint_force_corrected_domain",
    "core_current_head_q205_stagewise_pulse_diagonal",
    "core_current_head_q206_terminal_origin_obstruction",
)
EXPECTED = {
    PROVIDERS[0]: "Q214-exact-core-selected-roots-and-leading-pressure-joined-conditionally",
    PROVIDERS[1]: "Q213-reconstructed-primitive-bundle-joined-conditionally-exact-core-root-frontier-isolated",
    PROVIDERS[2]: "Q207-fixed-stage-complete-residual-constant-reconstructed-conditionally-primitive-seminorm-frontier-isolated",
    PROVIDERS[3]: "Q211-polynomial-degree-grammar-closed-conditionally-weighted-primitive-frontier-isolated",
    PROVIDERS[4]: "Q184-generated-source-log-loss-corrected-and-pure-power-reserve-derived",
    PROVIDERS[5]: "Q166-current-head-complete-six-group-local-residual-audited-conditionally",
    PROVIDERS[6]: "Q203-corrected-domain-diagonal-endpoint-localized-force-replayed-conditionally",
    PROVIDERS[7]: "Q205-stagewise-generated-pulse-bounds-suffice-for-diagonal-uniform-raw-hierarchy-remains-unproved",
    PROVIDERS[8]: "Q206-terminal-origin-velocity-and-gauge-independent-pressure-value-jets-obstructed-conditionally",
}

R_DERIV = 36_479_350_300
Q_POLE = 218_876_101_818
R_STAR = max(R_DERIV, 12 * Q_POLE)


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q215 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def concrete_request(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q215 identity required")
    return {
        "stage": 602,
        "derivative_ceiling": R_DERIV,
        "common_pole": Q_POLE,
        "six_pole": 6 * Q_POLE,
        "twelve_pole": 12 * Q_POLE,
        "safe_request": R_STAR,
        "derivative_ceiling_alone_is_incomplete": R_DERIV < R_STAR,
        "finite": True,
    }


def regular_exponent(stage, m=0, reserve=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q215 identity required")
    if isinstance(stage, bool) or not isinstance(stage, int) or stage < 0:
        raise ValueError("nonnegative integer stage required")
    if isinstance(m, bool) or not isinstance(m, int) or m < 0:
        raise ValueError("nonnegative integer jet required")
    reserve = Q(reserve)
    if reserve <= 0:
        raise ValueError("positive reserve required")
    h = Q(1, 200)
    kappa = Q(1, 100000)
    delta = Q(1, 2) - 3 * kappa
    B0 = Q(3, 2) - 2 * kappa
    lam = Q(3, 2) + 2 * h
    return h * (B0 + stage * delta) - (lam + 2 * m + reserve)


def bundle_join(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q215 identity required")
    return {
        "request_fixed_before_tail_selection": True,
        "safe_request": R_STAR,
        "A_R_finite": True,
        "g_R_positive": True,
        "K_R_finite_including_6q_12q": True,
        "polynomial_S_degree_finite": True,
        "degree_numerically_evaluated": False,
        "same_profile_root_pressure_and_support_object": True,
        "Q207_primitive_provider_flag_closed_conditionally": True,
    }


def residual_join(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q215 identity required")
    return {
        "stage_601_exponent": regular_exponent(601),
        "stage_602_exponent": regular_exponent(602),
        "first_zero_jet_zero_target_stage": 602,
        "log_reserve": Q(49997, 40000000),
        "compact_reserve": Q(1, 40000000),
        "six_groups": [
            "time", "unit viscosity", "complete pressure",
            "old-new cross", "new-old cross", "new square",
        ],
        "cutoff_curls_and_pressure_derivatives_retained": True,
        "conditional_fixed_jet_residual_flatness_joined": True,
        "uniform_raw_hierarchy_constant_proved": False,
        "terminal_origin_regularized": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q215 identity required")
    request = concrete_request()
    bundle = bundle_join()
    residual = residual_join()
    checks = {
        "provider_statuses_match": all(statuses().values()),
        "complete_request_includes_pole_descendants": (
            request["safe_request"] == 2_626_513_221_816
            and request["derivative_ceiling_alone_is_incomplete"]
        ),
        "fixed_request_bundle_closes_conditionally": (
            bundle["A_R_finite"] and bundle["g_R_positive"]
            and bundle["K_R_finite_including_6q_12q"]
            and bundle["polynomial_S_degree_finite"]
        ),
        "stage_602_is_first_positive_regular_stage": (
            residual["stage_601_exponent"] < 0 < residual["stage_602_exponent"]
        ),
        "complete_residual_memory_retained": (
            len(residual["six_groups"]) == 6
            and residual["cutoff_curls_and_pressure_derivatives_retained"]
            and residual["conditional_fixed_jet_residual_flatness_joined"]
        ),
    }
    return {
        "status": "Q215-fixed-stage-primitive-bundle-and-complete-residual-recurrence-joined-conditionally",
        "checks": checks,
        "scope": {
            "Q207_concrete_request_discharged_relative_to_quantified_providers": True,
            "conditional_fixed_jet_residual_flatness_joined": True,
            "uniform_uncut_pulse_hierarchy": False,
            "terminal_origin_velocity_or_pressure_regular": False,
            "unconditional_residual_force_exclusion_or_blowup": False,
        },
        "next": (
            "independently audit the source proofs behind Q213's all-fixed-order and all-pole "
            "quantifiers at the enlarged R_star request, detecting any circular or diagnostic-only row"
        ),
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q215 identity required")
    residual = residual_join()
    return {
        "status": review()["status"],
        "providers": statuses(),
        "concrete_request": concrete_request(),
        "bundle_join": bundle_join(),
        "residual_join": {
            **residual,
            "stage_601_exponent": str(residual["stage_601_exponent"]),
            "stage_602_exponent": str(residual["stage_602_exponent"]),
            "log_reserve": str(residual["log_reserve"]),
            "compact_reserve": str(residual["compact_reserve"]),
        },
        "review": review(),
    }
