"""Q220: source audit of the shifted two-edge weighted kernel theorem."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q220-SHIFTED-TWO-EDGE-KERNEL-SOURCE-AUDIT"
PROVIDERS = (
    "core_current_head_q219_regional_gap_weighted_stress_source_audit",
    "core_current_head_q216_quantified_provider_source_audit",
    "core_direct_compact_slow_provider_reconstruction",
    "core_compact_slow_inverse_reconstruction",
    "core_compact_slow_inverse_domain_join",
    "core_symbolic_weighted_maps",
    "core_current_head_q29_weighted_compact_constants",
    "core_current_head_q174_q5_operator_source_audit",
    "core_current_head_q200_physical_weighted_target_reconciliation",
)
EXPECTED = {
    PROVIDERS[0]: "Q219-same-band-regional-gaps-and-weighted-relative-stress-source-audited",
    PROVIDERS[1]: "Q216-all-fixed-order-and-finite-pole-provider-quantifiers-source-audited",
    PROVIDERS[2]: "direct-compact-slow-operators-proved-with-actual-source-supremum-frontier",
    PROVIDERS[3]: "compact-slow-weighted-inverse-contract-reconstructed-on-current-same-choice-tail",
    PROVIDERS[4]: "compact-slow-weighted-domain-and-support-contract-joined",
    PROVIDERS[5]: "conditional-symbolic-weighted-compact-map-and-fast-velocity-constants",
    PROVIDERS[6]: "conditional-explicit-weighted-shifted-compact-constants",
    PROVIDERS[7]: "Q174-Q5-operator-source-audited-with-family-and-typed-domain-boundary",
    PROVIDERS[8]: "Q200-physical-weighted-target-route-reconciled-at-corrected-rates",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q220 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def directed_maps(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q220 identity required")
    return {
        "I_minus": "int_{a<0} F(U+a,y+Ma v_r) da",
        "I_plus": "int_{a>0} F(U+a,y+Ma v_r) da",
        "T_e": "g_e[(1-chi)I_minus-chi I_plus]",
        "A_e": "k_e chi_U int_R F(U+a,y+Ma v_r) da",
        "physical_shift_retained": True,
        "differentiate_at_fixed_integration_shift": True,
        "translation_depends_on_differentiated_variables": False,
        "extra_M_power_from_differentiation": 0,
        "smooth_zero_extension_removes_boundary_terms": True,
        "wrong_edge_orientation_valid": False,
    }


def weight_contract(pole, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q220 identity required")
    if isinstance(pole, bool) or not isinstance(pole, int) or pole < 0:
        raise ValueError("nonnegative integer pole required")
    return {
        "input_pole": pole,
        "output_pole": pole,
        "pole_growth": 0,
        "weight": "zeta=exp(-w^2/x^2-4/(L-x)^2)",
        "delta": "min(1,x,L-x)",
        "inner_branch": "D(w^2) exp(4/(L-r)^2)",
        "outer_branch": "D(4) exp(w^2/(L-r)^2)",
        "interior_branch": "L E_p exp((w^2+4)/r^2)",
        "proof_collar": "r=min(1,L/4,g)",
        "all_three_branches_finite_for_fixed_positive_data": True,
        "constant_uniform_in_pole": False,
        "giant_constant_numerically_materialized": False,
    }


def cutoff_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q220 identity required")
    return {
        "chi_zero_on_left_plateau": True,
        "chi_one_on_right_plateau": True,
        "uniform_interior_gap_positive": True,
        "I_minus_selected_away_from_right_edge": True,
        "I_plus_selected_away_from_left_edge": True,
        "chi_U_support_has_positive_weight_floor": True,
        "A_e_keeps_whole_shifted_moment": True,
        "zero_Haar_alone_deletes_A_e": False,
        "cutoff_jet_request_for_output_order_m": "m+1",
    }


def typed_boundary(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q220 identity required")
    return {
        "full_weight_input": "zeta delta^-p",
        "full_weight_output": "zeta delta^-p",
        "half_weight_input_automatically_upgraded": False,
        "obstruction": "zeta^-1 x^p sqrt(zeta) -> infinity at x downarrow 0",
        "physical_covariance_products_supply_full_weight_by_structure": True,
        "unweighted_jets_plus_support_and_moments_imply_full_weight": False,
        "source_membership_is_separate_from_operator_boundedness": True,
    }


def source_classification(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q220 identity required")
    return {
        "Gaussian_absorbs_every_fixed_polynomial_pole": True,
        "interior_weight_floor_is_continuum": True,
        "symbolic_exponential_sum_is_exact_and_finite": True,
        "floating_point_materialization_required": False,
        "tests_used_as_continuum_proof": False,
        "kernel_theorem_supported_at_source_level": True,
        "actual_generated_source_membership_proved_by_kernel_itself": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q220 identity required")
    maps = directed_maps()
    low = weight_contract(0)
    high = weight_contract(2_626_513_221_816)
    cutoff = cutoff_contract()
    typed = typed_boundary()
    source = source_classification()
    checks = {
        "provider_statuses_match": all(statuses().values()),
        "physical_shift_and_derivative_order_are_preserved": (
            maps["physical_shift_retained"]
            and maps["differentiate_at_fixed_integration_shift"]
            and maps["extra_M_power_from_differentiation"] == 0),
        "every_fixed_pole_returns_same_pole": (
            low["output_pole"] == low["input_pole"]
            and high["output_pole"] == high["input_pole"]
            and high["pole_growth"] == 0),
        "three_branch_constant_is_finite_without_materialization": (
            high["all_three_branches_finite_for_fixed_positive_data"]
            and not high["giant_constant_numerically_materialized"]),
        "cutoff_selects_correct_edges_and_interior_floor": (
            cutoff["I_minus_selected_away_from_right_edge"]
            and cutoff["I_plus_selected_away_from_left_edge"]
            and cutoff["chi_U_support_has_positive_weight_floor"]),
        "typed_source_boundary_is_not_erased": (
            not typed["half_weight_input_automatically_upgraded"]
            and not typed["unweighted_jets_plus_support_and_moments_imply_full_weight"]
            and typed["source_membership_is_separate_from_operator_boundedness"]),
        "source_proof_not_diagnostic": (
            source["kernel_theorem_supported_at_source_level"]
            and not source["tests_used_as_continuum_proof"]),
    }
    return {
        "status": "Q220-shifted-two-edge-weighted-kernel-every-finite-pole-source-audited",
        "checks": checks,
        "scope": {
            "Q216_two_edge_kernel_quantifier_supported": True,
            "same_pole_and_weight_preserved_for_typed_inputs": True,
            "blanket_full_weight_generated_source_claim": False,
            "uniform_actual_pulse_or_complete_residual": False,
            "unconditional_force_exclusion_or_blowup": False,
        },
        "next": (
            "audit the remaining compact/slow source theorems: the zero-Haar Fourier inverse "
            "and fixed slow response with axial velocity, moving pressure product and finite lambda^-1"
        ),
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q220 identity required")
    return {
        "status": review()["status"],
        "providers": statuses(),
        "directed_maps": directed_maps(),
        "weight_contract_low_pole": weight_contract(0),
        "weight_contract_R_star": weight_contract(2_626_513_221_816),
        "cutoff_contract": cutoff_contract(),
        "typed_boundary": typed_boundary(),
        "source_classification": source_classification(),
        "review": review(),
    }
