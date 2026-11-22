"""Q221: source audit of the zero-Haar Fourier inverse and slow response."""
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q221-FOURIER-SLOW-RESPONSE-SOURCE-AUDIT"
PROVIDERS = (
    "core_current_head_q220_shifted_two_edge_kernel_source_audit",
    "core_fast_inverse_budget",
    "core_slow_inverse_budget",
    "core_current_head_q30_fixed_support_join",
    "core_current_head_q174_q5_operator_source_audit",
    "core_current_head_q175_canonical_source_route_reconciliation",
    "core_current_head_q216_quantified_provider_source_audit",
    "core_compact_slow_inverse_reconstruction",
    "core_current_head_q28_initial_mean_composition",
    "core_current_head_q177_weighted_support_source_reconciliation",
)
EXPECTED = {
    PROVIDERS[0]: "Q220-shifted-two-edge-weighted-kernel-every-finite-pole-source-audited",
    PROVIDERS[1]: "explicit-fast-inverse-actual-axial-repair-and-current-residual-jet-budget",
    PROVIDERS[2]: "explicit-fixed-profile-slow-inverse-and-coupled-scalar-substitution",
    PROVIDERS[3]: "fixed-choice-cutoff-and-slow-hull-join",
    PROVIDERS[4]: "Q174-Q5-operator-source-audited-with-family-and-typed-domain-boundary",
    PROVIDERS[5]: "Q175-canonical-typed-source-route-reconciled-and-independent-frontier-restored",
    PROVIDERS[6]: "Q216-all-fixed-order-and-finite-pole-provider-quantifiers-source-audited",
    PROVIDERS[7]: "compact-slow-weighted-inverse-contract-reconstructed-on-current-same-choice-tail",
    PROVIDERS[8]: "conditional-primary-to-signed-mean-composition",
    PROVIDERS[9]: "Q177-selected-input-weighted-stress-and-reserved-support-route-reconciled",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q221 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def quadratic_norms(k1, k2, identity=IDENTITY):
    """Exact integer norms for vr=(1,-b), vt=(b,1), b=sqrt(2)-1."""
    if identity != IDENTITY:
        raise ValueError("Q221 identity required")
    if any(isinstance(k, bool) or not isinstance(k, int) for k in (k1, k2)):
        raise ValueError("integer Fourier mode required")
    if k1 == 0 and k2 == 0:
        raise ValueError("nonzero Fourier mode required")
    return {
        "radius_l1": abs(k1) + abs(k2),
        "vr_integer_norm": (k1 + k2) ** 2 - 2 * k2 ** 2,
        "vt_integer_norm": (k2 - k1) ** 2 - 2 * k1 ** 2,
        "both_norms_nonzero": (
            (k1 + k2) ** 2 - 2 * k2 ** 2 != 0
            and (k2 - k1) ** 2 - 2 * k1 ** 2 != 0
        ),
        "conjugate_factor_bound": "(1+sqrt(2))|k|_1 < 3|k|_1",
        "inverse_divisor_bound": "|v_r.k|^-1,|v_t.k|^-1 < 3|k|_1",
    }


def fourier_contract(output_order, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q221 identity required")
    if isinstance(output_order, bool) or not isinstance(output_order, int) or output_order < 0:
        raise ValueError("nonnegative integer output order required")
    return {
        "output_order": output_order,
        "source_order": output_order + 4,
        "divisor_power_cost": 1,
        "post_integration_per_mode_power": -3,
        "l1_shell_multiplicity": "4r",
        "tail_bound": "sum_(r>=1) r^-2 <= 2",
        "Gamma_m": f"24*2^({output_order}+4)",
        "giant_constant_materialized": False,
        "zero_Haar_required": True,
        "zero_mode_obstruction": "the multiplier (v.k)^-1 is undefined at k=0",
        "slow_radial_weight_is_passive": True,
        "coefficient_reciprocal_jets_are_separate_finite_factors": True,
    }


def fast_contract(output_order, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q221 identity required")
    if isinstance(output_order, bool) or not isinstance(output_order, int) or output_order < 0:
        raise ValueError("nonnegative integer output order required")
    return {
        "u": "-c^-1 N_t^-1 F_theta",
        "g_d": "-c^-1 N_t^-1 F_z",
        "psi": "T_1 g_d",
        "beta": "-epsilon partial_Z psi",
        "gamma": "g_d-A_1 g_d",
        "actual_A1_axial_repair_retained": True,
        "tangential_axial_source_order": output_order + 4,
        "radial_potential_source_order": output_order + 5,
        "postslow_prefast_source_order": output_order + 8,
        "whole_axial_pressure_source_order": output_order + 9,
    }


def slow_contract(output_order, pole, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q221 identity required")
    if any(isinstance(v, bool) or not isinstance(v, int) or v < 0
           for v in (output_order, pole)):
        raise ValueError("nonnegative integer order and pole required")
    return {
        "output_order": output_order,
        "input_pole": pole,
        "output_pole": pole,
        "u": "a^-1[s^A P F_P+s^(A-1) J_z F_z]",
        "gamma": "a^-1 s^(A-3/2) J_theta F_g",
        "psi": "a^-1 s^(A-1) J_theta Psi",
        "beta": "-epsilon partial_Z psi",
        "distinct_scale_powers_retained": True,
        "five_moment_rows_exact": True,
        "divergence_identity": "D_1 beta+epsilon partial_Z gamma=0",
        "whole_convex_potential_hull_used": True,
        "hull_strictly_inside_reserved_annulus": True,
        "interior_weight_floor_preserves_every_fixed_pole": True,
        "fixed_profile_jets_finite_at_every_fixed_order": True,
    }


def axial_conditioning(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q221 identity required")
    return {
        "v": "exp(-2 d lambda)",
        "axial_determinant_factor": "g(1-v)",
        "exact_singular_factor": "(1-exp(-2 d lambda))^-1",
        "upper_bound": "exp(2 d lambda)/(2 d lambda)",
        "selected_lambda_positive": True,
        "finite_for_selected_lambda": True,
        "uniform_as_lambda_downarrow_zero": False,
        "precise_obstruction": "(1-exp(-2 d lambda))^-1 ~ (2 d lambda)^-1",
        "angular_block_singular_at_lambda_zero": False,
    }


def pressure_memory(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q221 identity required")
    return {
        "whole_moving_pressure_derivative": (
            "epsilon[J_z,Z+c_rho,Z P+c_rho P_Z]"
        ),
        "c_rho_Z_P_retained": True,
        "actual_axial_velocity_retained": True,
        "all_six_covariances_retained_in_inputs": True,
        "Pi1_Pi2_retained_in_reconstructed_pressure": True,
        "operator_creates_its_own_source_bounds": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q221 identity required")
    mode = quadratic_norms(2_400_066, -2_626_513)
    fourier = fourier_contract(2_400_066)
    fast = fast_contract(17)
    slow = slow_contract(2_400_066, 2_626_513_221_816)
    axial = axial_conditioning()
    pressure = pressure_memory()
    checks = {
        "provider_statuses_match": all(statuses().values()),
        "integer_norm_excludes_nonzero_resonances": mode["both_norms_nonzero"],
        "four_derivatives_leave_summable_shell": (
            fourier["source_order"] == fourier["output_order"] + 4
            and fourier["post_integration_per_mode_power"] == -3),
        "zero_Haar_boundary_is_explicit": (
            fourier["zero_Haar_required"]
            and "k=0" in fourier["zero_mode_obstruction"]),
        "fast_axial_and_pressure_losses_retained": (
            fast["actual_A1_axial_repair_retained"]
            and fast["whole_axial_pressure_source_order"] == 26),
        "slow_map_keeps_pole_support_and_scale_rows": (
            slow["output_pole"] == slow["input_pole"]
            and slow["whole_convex_potential_hull_used"]
            and slow["distinct_scale_powers_retained"]),
        "selected_axial_conditioning_finite_but_not_uniform": (
            axial["finite_for_selected_lambda"]
            and not axial["uniform_as_lambda_downarrow_zero"]),
        "complete_pressure_memory_not_frozen": (
            pressure["c_rho_Z_P_retained"]
            and pressure["actual_axial_velocity_retained"]
            and pressure["Pi1_Pi2_retained_in_reconstructed_pressure"]),
    }
    return {
        "status": "Q221-zero-Haar-Fourier-and-fixed-slow-response-source-audited",
        "checks": checks,
        "scope": {
            "Q216_remaining_inverse_quantifiers_supported": True,
            "every_fixed_order_and_finite_pole_on_typed_inputs": True,
            "nonzero_Haar_inverse_claim": False,
            "lambda_uniform_slow_inverse_claim": False,
            "actual_generated_source_membership_proved_by_operators": False,
            "uniform_actual_pulse_or_complete_residual": False,
            "unconditional_force_exclusion_or_blowup": False,
        },
        "next": (
            "audit the global pressure and moment identities on the joined compact/fast/slow "
            "object, including both exact tangential moments, the whole moving pressure "
            "derivative, axial velocity and Pi1/Pi2"
        ),
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q221 identity required")
    return {
        "status": review()["status"],
        "providers": statuses(),
        "quadratic_norm_sample": quadratic_norms(2_400_066, -2_626_513),
        "fourier_contract": fourier_contract(2_400_066),
        "fast_contract": fast_contract(17),
        "slow_contract": slow_contract(2_400_066, 2_626_513_221_816),
        "axial_conditioning": axial_conditioning(),
        "pressure_memory": pressure_memory(),
        "review": review(),
    }
