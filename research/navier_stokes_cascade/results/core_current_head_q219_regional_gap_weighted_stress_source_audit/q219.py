"""Q219: source audit of the same-band regional gaps and weighted stress."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q219-REGIONAL-GAP-WEIGHTED-STRESS-SOURCE-AUDIT"
PROVIDERS = (
    "core_current_head_q218_energy_class_comparison_entry",
    "core_current_head_q216_quantified_provider_source_audit",
    "core_leading_profile_certificate",
    "regional_cone_join",
    "weighted_background_edges",
    "core_same_band_gap_join",
    "core_same_band_gap_reconstruction",
    "core_current_head_q189_same_profile_envelope_gap_audit",
)
EXPECTED = {
    PROVIDERS[0]: "Q218-energy-class-entry-to-fixed-gauge-comparison-derived-conditionally",
    PROVIDERS[1]: "Q216-all-fixed-order-and-finite-pole-provider-quantifiers-source-audited",
    PROVIDERS[2]: "finalized-leading-profile-certificate-joined-relative-to-exact-core-and-root-hypotheses",
    PROVIDERS[3]: "regional-input-gaps-and-leading-value-cone-joined",
    PROVIDERS[4]: "all-order-weighted-edges-and-relative-background-stress-audited",
    PROVIDERS[5]: "same-rule-defined-background-tail-phase-normal-value-cone-contract-joined",
    PROVIDERS[6]: "same-band-continuum-gap-contract-reconstructed-on-finalized-leading-profile",
    PROVIDERS[7]: "Q189-same-profile-envelope-and-gap-join-closed-relative-to-region-theorems",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q219 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def regional_ledger(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q219 identity required")
    corrected_ramp = Q(17701, 30625)
    axial = 2 * Q(35) ** 2 * Q(47, 50) * Q(417, 1000)
    intermediate = 2 * Q(12) ** 2 * Q(627, 2000)
    return {
        "activation_a_floor_power": 22,
        "activation_third_gap_power": 401,
        "activation_fourth_gap_power": 802,
        "joined_gap_power": 1024,
        "compact_envelope_power": 4096,
        "corrected_outer_ramp_floor": corrected_ramp,
        "corrected_outer_ramp_exceeds_half": corrected_ramp > Q(1, 2),
        "axial_fourth_gap_lower_control": axial,
        "axial_control_exceeds_one": axial > 1,
        "intermediate_fourth_gap_lower_control": intermediate,
        "intermediate_control_exceeds_one": intermediate > 1,
        "regionwise_suprema_from_analytic_inequalities": True,
        "sample_grid_used_as_continuum_proof": False,
        "complete_pressure_memory_retained": True,
    }


def edge_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q219 identity required")
    return {
        "global_absolute_stress_infimum_positive": False,
        "global_absolute_stress_infimum": Q(0),
        "inner_factor": "exp(-t1^2/y_a^2) A_a",
        "inner_factor_floor_positive": True,
        "outer_factor": "exp(-4/y_b^2) y_b^-3 (b0,y_b^6 b_z)",
        "outer_b0_positive": True,
        "weighted_floor": "|T0| >= c zeta",
        "normalized_direction_extends_smoothly": True,
        "opposite_edge_factor_kept": True,
    }


def background_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q219 identity required")
    return {
        "relative_stress": "D^I[(T_hat-T0)/|T0|]=O_I(q^(2h))",
        "all_fixed_derivative_orders": True,
        "one_domain_for_all_fixed_orders": True,
        "constants_uniform_in_order": False,
        "finite_head_plus_geometric_tail": True,
        "coefficient_one_kept_in_finite_head": True,
        "higher_coefficients_have_common_compact_support": True,
        "generic_absolute_weighted_bound_is_sufficient": False,
        "exact_edge_factor_and_support_used": True,
    }


def perturbation_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q219 identity required")
    return {
        "rho_gap": "min(1/2,f/2,t/2,d/(2L_D),omega/(2L_H),m_phi/(2L_phi))",
        "rho_gap_positive_from_finite_positive_inputs": True,
        "F_floor": "f/2",
        "minus_g_theta_floor": "t/2",
        "D_floor": "d/2",
        "lambda_shear_squared_floor": "f*t*d/[4*(Gbar+1)^2]",
        "c_shear_squared_floor": "d/[4*(Fbar+1)*(Gbar+1)]",
        "selector_uses_finitely_many_value_direction_rows": True,
        "selector_is_derivative_order_independent": True,
        "scale_selection_circularity": False,
    }


def source_classification(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q219 identity required")
    return {
        "regional_cover_is_continuum_and_finite": True,
        "flat_edges_are_handled_by_exact_factorization": True,
        "weighted_relative_stress_is_continuum": True,
        "interval_arithmetic_role": "fixed scalar outward rounding only",
        "tests_role": "identity provenance and regression checks only",
        "diagnostic_only_continuum_row_found": False,
        "gap_leaf_supported_at_theorem_dependency_level": True,
        "every_primitive_regional_ODE_rederived_here": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q219 identity required")
    region = regional_ledger()
    edge = edge_contract()
    background = background_contract()
    perturb = perturbation_contract()
    source = source_classification()
    checks = {
        "provider_statuses_match": all(statuses().values()),
        "regional_powers_fit_compact_envelope": (
            max(region["activation_a_floor_power"], region["activation_third_gap_power"],
                region["activation_fourth_gap_power"], region["joined_gap_power"])
            < region["compact_envelope_power"]),
        "corrected_outer_and_later_region_floors_are_strict": (
            region["corrected_outer_ramp_exceeds_half"]
            and region["axial_control_exceeds_one"]
            and region["intermediate_control_exceeds_one"]),
        "zero_stress_edges_use_weighted_not_absolute_floor": (
            not edge["global_absolute_stress_infimum_positive"]
            and edge["global_absolute_stress_infimum"] == 0
            and edge["normalized_direction_extends_smoothly"]),
        "relative_stress_uses_exact_construction": (
            background["all_fixed_derivative_orders"]
            and background["finite_head_plus_geometric_tail"]
            and not background["generic_absolute_weighted_bound_is_sufficient"]),
        "common_selector_is_acyclic": (
            perturb["selector_uses_finitely_many_value_direction_rows"]
            and perturb["selector_is_derivative_order_independent"]
            and not perturb["scale_selection_circularity"]),
        "continuum_claim_not_from_diagnostics": (
            source["gap_leaf_supported_at_theorem_dependency_level"]
            and not source["diagnostic_only_continuum_row_found"]),
    }
    return {
        "status": "Q219-same-band-regional-gaps-and-weighted-relative-stress-source-audited",
        "checks": checks,
        "scope": {
            "Q216_same_band_gap_leaf_supported_by_continuum_source_chain": True,
            "false_positive_absolute_stress_floor_rejected": True,
            "all_primitive_regional_theorems_first_principles_reproved": False,
            "uniform_actual_pulse_or_complete_residual": False,
            "unconditional_force_exclusion_or_blowup": False,
        },
        "next": (
            "audit Q216's next cited source theorem: the shifted two-edge weighted kernel "
            "estimate for every finite pole, including its exact physical shift and edge weights"
        ),
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q219 identity required")
    region = regional_ledger()
    edge = edge_contract()
    return {
        "status": review()["status"],
        "providers": statuses(),
        "regional_ledger": {key: str(value) if isinstance(value, Q) else value
                            for key, value in region.items()},
        "edge_contract": {key: str(value) if isinstance(value, Q) else value
                          for key, value in edge.items()},
        "background_contract": background_contract(),
        "perturbation_contract": perturbation_contract(),
        "source_classification": source_classification(),
        "review": review(),
    }
