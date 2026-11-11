"""Audit Q163's U1 core contraction and leading-pressure obligation."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
IDENTITY = "Q164-all-index-core-contraction-leading-pressure-current-head-audit"
PROVIDERS = (
    "core_current_head_q163_comparison_core_exclusion_join",
    "core_contraction_audit",
    "core_fixed_point_ledger_recheck",
    "core_pressure_identification",
    "core_same_field_pressure_reconciliation",
    "core_current_head_q1_exact_map_join",
)


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


ledger = load(
    "q164_fixed_point_ledger",
    "core_fixed_point_ledger_recheck/proof_checks.py",
)


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def operator_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q164 identity required")
    checks = ledger.review()["checks"]
    return {
        "weight": (
            "a_nk=20^-n rho^-k k! binom(n+k,k)/"
            "((n+1)^2(k+1)^2)"
        ),
        "product_constant_A0": ledger.A0,
        "integrated_derivative_constant_D0": ledger.D0_SYMBOL,
        "angular_inverse_word": "30720^k/(k!(k+1)!)",
        "angular_inverse_bound": "V=exp(30720)",
        "all_index_product_bound_rederived": True,
        "all_index_mixed_allocation_bound_rederived": True,
        "finite_samples_are_not_all_index_proofs": True,
        "exact_checks": checks,
    }


def source_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q164 identity required")
    rows = ledger.ledger()
    sums = ledger.group_sums(rows)
    major = ledger.majorant_coefficients()
    lipschitz = ledger.lipschitz_coefficients()
    return {
        "term_count": len(rows),
        "unique_term_names": len({row["name"] for row in rows}),
        "group_sums": sums,
        "pressure_primitive": major["pressure_primitive_A0_cubed"],
        "pressure_source": major["R2_pressure_A0_fourth"],
        "pressure_lipschitz_per_S": lipschitz["R2_pressure_per_S"],
        "complete_pressure_memory_retained": (
            {row["name"] for row in rows if row["group"] == "R2_pressure"}
            == {"R2_minus_4A_eta_p", "R2_d_p_eta", "R2_minus_2eta_Dp"}
        ),
    }


def contraction_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q164 identity required")
    universal = ledger.universal_contraction()
    eps = Q(1, 10**16)
    return {
        "epsilon": eps,
        "H_condition": "H=K^(1/16)>1",
        "Lambda": "10^12 H^128 epsilon^-2",
        "selected_bound": "epsilon^2/(10^18 H^127)",
        "coefficient_before_H": universal["epsilon_squared_over_10_power_18"],
        "H_denominator_power": universal["H_denominator_power"],
        "strict_universal_ceiling": universal["strict_universal_ceiling"],
        "strictly_below_ceiling_because_H_gt_one": True,
        "ball_invariant": True,
        "strict_contraction": True,
        "finite_jet_iteration_used_as_convergence_proof": False,
    }


def equation_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q164 identity required")
    return {
        "map": (
            "(Phi,u)->(Phi0+(1+T)^-1 J2 R1/(2 Lambda),"
            "u0+J1 R2/(2 Lambda))"
        ),
        "angular_equation": "-2L(X phi_XX+2 phi_X)/phi=S_q",
        "axial_equation": "-2L(X U_XX+U_X)=S_n",
        "pressure_identity": "Pi_X=phi^2/C^2",
        "regular_forward_candidates": (
            "Q=-2L phi_X/phi; N=-2L U_X"
        ),
        "homogeneous_forward_modes": "constant/(XH) and constant/X",
        "axis_regularity_forces_both_constants_zero": True,
        "exact_leading_core_stress": "T0 identically zero",
        "pressure_or_one_of_its_derivatives_removed": False,
    }


def pressure_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q164 identity required")
    A = Q(101, 200)
    return {
        "A": A,
        "defect_equation": "4A eta delta-(1-eta^2) delta'=0",
        "homogeneous_mode": "C*(1-eta^2)^(-2A)",
        "positive_exponent": A > 0,
        "nonzero_mode_endpoint_bounded": False,
        "smooth_endpoint_compatibility_forces_C_zero": True,
        "same_gauge_identity": "Pi_axis+C_p,infinity=0",
        "leading_core_equals_completed_leading_backward_pressure": True,
        "full_corrected_pressure_collapsed_to_leading_pressure": False,
        "later_background_wave_mean_localization_pressure_retained": True,
    }


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q164 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    expected = {
        PROVIDERS[0]: "Q163-same-object-comparison-and-core-exclusion-joined-conditionally",
        PROVIDERS[1]: "core-contraction-and-exact-stress-audited-for-pressure-family",
        PROVIDERS[2]: "core-fixed-point-all-index-ledger-independently-rederived",
        PROVIDERS[3]: "leading-core-pressure-identified-with-regular-angular-inverse-and-retained-higher-memory",
        PROVIDERS[4]: "same-field-core-pressure-identified-with-completed-leading-backward-pressure",
        PROVIDERS[5]: "current-head-Q1-exact-core-roots-and-leading-pressure-bound-on-one-profile",
    }
    for name, status in expected.items():
        if data[name].get("status") != status:
            raise ValueError("Q164 provider status mismatch: " + name)
    q163 = data[PROVIDERS[0]]["symbolic"]
    contraction = data[PROVIDERS[1]]
    fixed = data[PROVIDERS[2]]
    pressure = data[PROVIDERS[3]]
    same = data[PROVIDERS[4]]
    q1 = data[PROVIDERS[5]]["review"]
    required = (
        q163["P7_comparison_and_core_exclusion_replayed_conditionally"],
        contraction["flags"]["all_index_integrated_derivative_and_product_bounds_derived"],
        contraction["flags"]["nonlinear_contraction_for_stated_exact_pressure_family_audited"],
        contraction["gates"]["original_core_equations_and_forward_stress_jets"],
        fixed["flags"]["all_normalized_R1_R2_monomials_enumerated"],
        fixed["flags"]["complete_pressure_primitive_and_three_pressure_terms_retained"],
        fixed["flags"]["selected_Lambda_formula_gives_strict_contraction"],
        pressure["flags"]["leading_pressure_identification_already_supplied_by_same_root_ledger"],
        not pressure["flags"]["leading_pressure_identification_still_an_independent_missing_equation"],
        same["review"]["checks"]["same_axis_and_zero_total_defect_give_forward_backward_equality_at_every_split"],
        same["review"]["scope"]["leading_q_zero_pressure_identified_on_same_field"],
        q1["checks"]["nonzero_pressure_homogeneous_mode_is_endpoint_unbounded"],
    )
    if not all(required):
        raise ValueError("Q164 core/pressure provider contract missing")
    return {
        "Q163_U1_was_stale_relative_to_live_provider_graph": True,
        "independent_all_index_rederivation_bound_to_Q1_current_head": True,
        "original_equations_and_regular_forward_stress_recovered": True,
        "same_gauge_leading_pressure_identified": True,
        "later_complete_pressure_memory_retained": True,
    }


def remaining_gate(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q164 identity required")
    return {
        "U1_all_index_core_contraction_original_equations_leading_pressure": (
            "closed-by-existing-independent-rederivations-and-Q164-current-head-join"
        ),
        "next": (
            "audit the finite defining, existence, inverse and gap maps on the same "
            "Q1-Q163 object, beginning with Q2's positive core, determinant, "
            "orientation, same-tilt cone and factored flat-edge continuum infima "
            "over the actual finite regional cover"
        ),
        "Q2_regional_continuum_infima_discharged": False,
        "Q3_full_frame_primitive_domain_discharged": False,
        "dependency_hashes_or_passing_tests_suffice": False,
        "unconditional_theorem_gate_closed": False,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q164 identity required")
    return (
        *({"id": f"P{i}", "status": "closed-relative-to-named-providers"}
          for i in range(1, 7)),
        {"id": "P7", "status": "closed-as-conditional-consequence-by-Q163"},
        {"id": "U1", "status": "closed-by-independent-core-pressure-audit"},
        {"id": "U2", "status": "next-independent-finite-map-audit"},
    )


def record():
    return {
        "status": "Q164-all-index-core-contraction-and-leading-pressure-audited",
        "correction": (
            "Q163's stated next gate was too pessimistic: the live branch already "
            "contained an independent all-index source/contraction rederivation and "
            "a same-field leading-pressure reconciliation; Q164 binds those records "
            "to the current-head Q1-Q163 object"
        ),
        "established": (
            "all 21 normalized source terms, the full three-term pressure memory, "
            "the invariant ball and strict <10^-50 contraction, original-equation "
            "recovery, regular forward stress cancellation and the same-gauge leading "
            "backward-pressure identity"
        ),
        "scope": (
            "analytical conclusion from the displayed all-index estimates relative to "
            "the selected exact pressure family; finite diagnostics and evidence hashes "
            "are regressions, and later corrected pressure terms remain"
        ),
        "next": remaining_gate()["next"],
        "U1_core_contraction_and_leading_pressure_closed": True,
        "Q2_or_Q3_independently_closed": False,
        "unconditional_pulse_residual_force_exclusion_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
    }
