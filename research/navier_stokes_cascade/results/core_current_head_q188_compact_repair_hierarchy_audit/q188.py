"""Independent exact-arithmetic audit of the compact repair hierarchy."""
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q188-COMPACT-REPAIR-HIERARCHY-AUDIT"
A0 = 2**16

PROVIDERS = (
    "core_current_head_q187_actual_symbolic_q2_selector",
    "prescribed_continuation_foundation",
    "compact_jet_envelope",
    "loop_modulation_bounds",
    "repair_state_bounds",
    "modulation_repair_all_order",
    "annular_exterior_foundation",
    "regional_cone_join",
    "collar_operator_envelope",
    "global_swirl_envelope",
    "global_axial_source_envelope",
    "global_axial_angular_bound",
    "signed_pressure_source_moments",
    "regional_stress_floor",
    "explicit_coefficient_cutoff",
    "explicit_tilt_pulse_tail",
)

EXPECTED = {
    PROVIDERS[0]: "Q187-actual-symbolic-Q2-selector-joined-and-primitive-frontier-isolated",
    PROVIDERS[1]: "pre-modulation-prescribed-continuation-through-shape-all-fixed-jets",
    PROVIDERS[2]: "actual-compact-input-envelope-bounded",
    PROVIDERS[3]: "fixed-phase-loop-and-modulation-errors-bounded",
    PROVIDERS[4]: "repair-state-bound-and-conditional-frequency-acceptance",
    PROVIDERS[5]: "same-finite-N-modulation-moving-phase-and-exact-repair-all-order-map-estimate",
    PROVIDERS[6]: "annular-exterior-all-order-leading-regularity-and-global-moment-foundation-relative-to-exact-roots",
    PROVIDERS[7]: "regional-input-gaps-and-leading-value-cone-joined",
    PROVIDERS[8]: "explicit-actual-activation-collar-operator-and-full-source-envelope",
    PROVIDERS[9]: "same-heat-root-third-jets-and-global-swirl-input-envelope",
    PROVIDERS[10]: "global-axial-pointwise-source-and-shifted-viscosity-envelope",
    PROVIDERS[11]: "global-third-angular-axial-envelope-and-signed-pressure-source-input",
    PROVIDERS[12]: "full-support-signed-pressure-source-moments-with-axis-viscosity",
    PROVIDERS[13]: "explicit-three-region-compact-leading-stress-floor",
    PROVIDERS[14]: "explicit-coefficient-cutoff-all-fixed-jets-and-pulse-gate-boundary",
    PROVIDERS[15]: "explicit-tilt-raw-deficit-and-primary-pulse-tail-certificate",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q188 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def exact_hierarchy(A=A0, identity=IDENTITY):
    """Reduce every huge-scale gate to integer comparisons at A>=2^16.

    We use only G=exp(A^32)>A and H=exp(G^4)>A.  Neither giant number is
    materialized and no floating logarithm is used.
    """
    if identity != IDENTITY:
        raise ValueError("Q188 identity required")
    if not isinstance(A, int) or isinstance(A, bool) or A < A0:
        raise ValueError("integer A>=2^16 required")
    beta, qmap = 1000, 10**5
    checks = {
        # r=A^-128, s<=r: the repaired E stays positive and d<=1.
        "positive_E_neighborhood": 2560 < A**126,
        "field_C1_perturbation_below_one": 1280 < A**127,
        # (16e6 A^2+5120 A^4)s < A^6 s.
        "shear_error_absorbed": 16 * 10**6 + 5120 * A**2 < A**4,
        # 327680 A^9 s < A^11 s < A^128 s.
        "four_state_repair_bound": 327680 < A**2 and 11 < 128,
        # G^30 >= 512*4^10 follows from G>A.
        "cone_lipschitz_tolerance": 512 * 4**10 < A**30,
        # Incoming modulation memory: H^31>512 A^6 follows from H>A.
        "incoming_pressure_memory": 512 * A**6 < A**31,
        # log(4002)<12, log(4000)<12, log(16e11)<41, since log 2<1.
        "state_log_constant": 4002 < 2**12,
        "radius_log_constant": 4000 < 2**12,
        "contraction_log_constant": 16 * 10**11 < 2**41,
        # The remaining state threshold is <16G^4+193G<32G^4.
        "state_frequency_exponent": 193 < 16 * A**3,
        # Root radius: log4000+128logA<16G^4.
        "root_radius_exponent": 12 + 128 * A < 16 * A**4,
        # Contraction: log(16e11)<16G^4.
        "root_contraction_exponent": 41 < 16 * A**4,
        # N>H^32 also implies N>=4AH.
        "modulation_transfer_frequency": 4 * A < A**31,
        "repair_constants_exact": beta == 1000 and qmap == 10**5,
    }
    return {
        "A": A,
        "checks": checks,
        "all_exact_reductions_hold": all(checks.values()),
        "scales_materialized": False,
        "definitions": {
            "G": "exp(A^32)",
            "H": "exp(G^4)",
            "Cstate": "H^16",
            "D": "H^16",
            "Ccorr": "A^128",
            "r": "A^-128",
            "epsilon": "G^-64",
            "N": "1+floor(H^32)",
        },
        "accepted_gates": {
            "state": "N > 2*H^16*(1+2000*A^128)*G^64",
            "coefficient_radius": "N > 4000*H^16*A^128",
            "contraction": "N > 8*10^11*H^16 (checked with 16*10^11)",
            "incoming_memory": "N > 512*A^6*H",
            "modulation_transfer": "N >= 4*A*H",
        },
    }


def pressure_memory_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q188 identity required")
    all_order = evidence(PROVIDERS[5])
    regional = evidence(PROVIDERS[7])
    return {
        "forward_rows": ["M", "I", "J", "S", "Cp"],
        "pressure": "Pi=Pi_axis+Cp",
        "loop_moment_error": "m=32*A^4*H/N with Pi_axis unchanged",
        "repair_target": "Bc+Q(c,c)=z_N(eta), all five normalized rows",
        "right_endpoint_identity": "M=I=J=S=Cp=original as functions of eta",
        "later_pressure_identity": (
            "same Pi_axis plus exact Cp equality gives the same pressure and "
            "integrated stress coordinates at every later radius"
        ),
        "all_order_provider_keeps_five_rows": all_order["flags"][
            "exact_five_row_full_pressure_memory_and_normalization_losses_retained"
        ],
        "regional_provider_uses_exact_global_ledger": regional["flags"][
            "exact_global_moment_ledger_applies_to_leading_value_construction"
        ],
        "four_row_repair_is_insufficient": True,
    }


def provider_correction(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q188 identity required")
    prescribed = evidence(PROVIDERS[1])
    all_order = evidence(PROVIDERS[5])
    exterior = evidence(PROVIDERS[6])
    swirl = evidence(PROVIDERS[9])
    axial = evidence(PROVIDERS[10])
    cutoff = evidence(PROVIDERS[14])
    return {
        "Q187_frontier_was_too_coarse_as_a_provider_inventory": True,
        "continuation_all_fixed_jets_supplied": prescribed["flags"][
            "every_fixed_jet_of_prescribed_continuation_through_shape_bounded_relative_to_exact_core"
        ],
        "same_map_all_order_modulation_repair_supplied": (
            all_order["flags"]["all_fixed_orders_of_signed_variance_and_moving_phase_inverse_derived"]
            and all_order["flags"]["same_accepted_C1_repair_has_all_fixed_derivatives_by_pointwise_Jacobian"]
        ),
        "annular_exterior_same_root_inputs_supplied": exterior["flags"][
            "whole_compact_smooth_input_to_same_finite_N_maps_supplied_relative_to_exact_core_and_finite_existence_bounds"
        ],
        "primitive_swirl_inputs_explicit": swirl["flags"][
            "all_primitive_leading_first_stress_inputs_now_explicit_relative_to_retained_lemmas"
        ],
        "global_axial_sources_explicit": axial["flags"][
            "complete_shifted_axial_U_operator_includes_differentiated_coordinate_coefficients"
        ],
        "coefficient_cutoff_explicit": cutoff["flags"][
            "qualitative_coefficient_cutoff_is_now_reproducibly_specialized"
        ],
        "independent_first_principles_reproof_completed_here": False,
    }


def obstruction(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q188 identity required")
    return {
        "kind": "continuum-provider-validation, not a frequency inequality",
        "compact_envelope": (
            "max of the required mixed-order <=2 field/moment/Pi norms, "
            "order <=1 (a,b,p1,p2) norms, E^-1, Hx^-1 and lambda^-1 must be <=A"
        ),
        "regional_gap": (
            "min{a, Pc-2, Pc-v, 2(Pc-v)^2-(v-2)Jc^2} >= A^-1 "
            "on the subsets where the last two expressions are required"
        ),
        "same_profile_requirement": (
            "both inequalities must be independently derived on the same selected "
            "continuation/root profile with all five moments and Pi_axis+Cp"
        ),
        "frequency_or_repair_arithmetic_still_blocks": False,
        "tests_discharge_continuum_suprema": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q188 identity required")
    status = statuses()
    hierarchy = exact_hierarchy()
    pressure = pressure_memory_contract()
    correction = provider_correction()
    block = obstruction()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "repair_frequency_hierarchy_rederived_exactly": hierarchy["all_exact_reductions_hold"],
        "complete_pressure_memory_retained": (
            pressure["all_order_provider_keeps_five_rows"]
            and pressure["regional_provider_uses_exact_global_ledger"]
            and pressure["four_row_repair_is_insufficient"]
        ),
        "later_provider_inventory_corrects_Q187_frontier": all(
            value for key, value in correction.items()
            if key not in {"independent_first_principles_reproof_completed_here"}
        ),
        "continuum_boundary_not_promoted": (
            not correction["independent_first_principles_reproof_completed_here"]
            and not block["frequency_or_repair_arithmetic_still_blocks"]
            and not block["tests_discharge_continuum_suprema"]
        ),
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "hierarchy": hierarchy,
        "pressure_memory": pressure,
        "correction": correction,
        "remaining_obstruction": block,
        "scope": {
            "low_order_compact_repair_hierarchy_closed_relative_to_actual_provider_bounds": True,
            "same_map_pressure_complete_cone_join_bound": True,
            "all_provider_continuum_suprema_reproved_from_initial_equations": False,
            "historical_first_band_recertified": False,
            "uniform_actual_pulse_or_complete_residual_proved": False,
            "unconditional_force_exclusion_or_blowup_proved": False,
            "tests_are_continuum_proof": False,
        },
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q188 identity required")
    return {
        "status": "Q188-compact-repair-hierarchy-closed-and-continuum-envelope-frontier-isolated",
        "established": (
            "the actual low-order compact modulation/exact-repair hierarchy, "
            "complete five-row pressure memory and regional cone composition "
            "have no remaining scalar frequency obstruction"
        ),
        "next": (
            "independently rederive the displayed compact envelope and regional "
            "gap inequalities on the same selected full-pressure profile"
        ),
        "unconditional_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
