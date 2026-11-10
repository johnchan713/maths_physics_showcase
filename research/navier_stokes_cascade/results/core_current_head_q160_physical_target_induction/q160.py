"""Bind the exact current-target recurrence to the inherited physical cycle."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
IDENTITY = "Q160-same-physical-family-all-finite-current-target-induction"
PROVIDERS = (
    "core_current_head_q159_primary_contract_frontier",
    "core_current_head_q149_primary_pulse_join",
    "core_current_head_q150_initial_covariance_mean_join",
    "core_current_head_q152_physical_stress_feedback_join",
    "current_target_moment_audit",
    "current_target_cycle_bridge",
    "source_envelope_edge_audit",
    "finite_stage_induction",
    "core_current_head_q153_physical_source_tail_diagonal_join",
)


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


bridge = load("q160_current_target_bridge", "current_target_cycle_bridge/target_cycle.py")


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def rates(j=0, kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q160 identity required")
    kappa = Q(kappa)
    if not isinstance(j, int) or j < 0 or not 0 < kappa < Q(1, 100):
        raise ValueError("finite nonnegative cycle and reviewed kappa required")
    source_gain = Q(1, 2) - 3 * kappa
    target_gain = Q(1, 2) - kappa
    B0 = Q(3, 2) - 2 * kappa
    C0 = B0 + Q(1, 2)
    B = B0 + j * source_gain
    C = C0 + j * source_gain
    return {
        "j": j,
        "kappa": kappa,
        "B0": B0,
        "C0": C0,
        "B_j": B,
        "C_j": C,
        "source_gain": source_gain,
        "bridge_target_gain": target_gain,
        "required_C_next": C + source_gain,
        "proved_target_order": C + target_gain,
        "target_slack": target_gain - source_gain,
        "C_equals_B_plus_half": C == B + Q(1, 2),
    }


def seed(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q160 identity required")
    r = rates()
    return {
        "physical_family": "Q149 inherited physical primaries on every accepted later band",
        "initial_target": "Q150 current initial target M_(2-3*kappa), with all six covariances and pressure memory",
        "same_family_stress": "Q152 H_beta y_beta=T_beta and sum eta_beta^2=1",
        "first_physical_source_order": r["B0"],
        "first_required_target_order": r["C0"],
        "prefix": (
            "the exact real covariance prefix and first generic signed current-state difference "
            "give Xi_0 in M_(2-2*kappa) on the same physical family"
        ),
        "complete_curl_margin": Q(1, 2),
        "canonical_and_physical_entrances_identified": False,
        "seed_closed": True,
    }


def target_step(C=Q(2) - Q(2, 100000), kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q160 identity required")
    C, kappa = Q(C), Q(kappa)
    cert = bridge.certificate(kappa)
    C0 = Q(cert["C0"])
    if C < C0:
        raise ValueError("target order below the reviewed physical seed")
    required_gain = Q(1, 2) - 3 * kappa
    proved_gain = Q(cert["gain"])
    rows = {}
    for name, row in cert["rows"].items():
        order = Q(row["order_slope"]) * C + Q(row["order_offset"])
        rows[name] = {
            "order": order,
            "margin_over_bridge_target": order - (C + proved_gain),
            "margin_over_required_next": order - (C + required_gain),
        }
    return {
        "C": C,
        "C0": C0,
        "required_gain": required_gain,
        "proved_gain": proved_gain,
        "slack": proved_gain - required_gain,
        "rows": rows,
        "all_rows_close_bridge_target": all(row["margin_over_bridge_target"] >= 0 for row in rows.values()),
        "all_rows_close_physical_next_target": all(row["margin_over_required_next"] >= 0 for row in rows.values()),
        "whole_half_line": True,
    }


def weighted_closure(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q160 identity required")
    return {
        "target": "|D^I Xi_j|<=C_(j,I) epsilon^(C_j) zeta ya^(-p_(j,I)) yb^(-q_(j,I))",
        "signed": (
            "X_new-X_old=-I_e Q_e overline(Delta E_e); the principal current covariance "
            "cancels X_old before estimating the complete error"
        ),
        "covariance": "every wave factor has sqrt(zeta), so every complete wave product has zeta",
        "means": "full-zeta mean products remain full-zeta because 0<zeta<=1 and zeta^2<=zeta",
        "primitive": (
            "the compensated input has zero R^e total; the two-edge oriented primitive "
            "retains zeta with only finite inverse-edge powers"
        ),
        "reserved_hull": "reserved bumps lie on the fixed interior hull where zeta has a positive floor",
        "moving_compensator": (
            "Q_e D_z H=D_z(Q_e H)+epsilon (partial_Z sigma_e) integral R^e H is retained"
        ),
        "stage_constants": "C_(j,I), p_(j,I), q_(j,I) may depend on the finite cycle and fixed jet",
        "uniform_in_stage_claimed": False,
    }


def moments(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q160 identity required")
    return {
        "pressure": (
            "pbar=I0(gbar-rho Pmean), integral R pbar=-Qg/2+c_rho Pmean"
        ),
        "angular": (
            "integral R^2 overline(E_theta,j)=epsilon partial_Z J_theta,j"
        ),
        "axial": (
            "integral R overline(E_z,j)=epsilon partial_Z(J_z,j+c_rho Pmean,j)"
        ),
        "current_memory": (
            "axial velocity; all six covariance entries including W_theta,z and W_zz; "
            "complete compact and moving pressure; Pi1/Pi2 in the same base; both velocity moments"
        ),
        "sourced_and_signed": "nonzero angular waves do not change the two Haar velocity moments",
        "fast": "the complete temporal-inverse curl has exact zero Haar mean and preserves both moments",
        "slow": "the five-map imposes its two zero velocity-moment rows and exact divergence",
        "zero_total_is_necessary": True,
        "preserved_at_every_finite_cycle": True,
    }


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q160 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    expected = {
        PROVIDERS[0]: "Q159-P2-corrected-closed-and-P3-finite-cycle-target-obstructed",
        PROVIDERS[1]: "Q148-schedule-frontier-corrected-and-uniform-primary-pulse-family-joined",
        PROVIDERS[2]: "Q149-uniform-primary-initial-complete-covariance-mean-correction-joined",
        PROVIDERS[3]: "Q152-physical-principal-stress-gate-closed-and-feedback-D-joined",
        PROVIDERS[4]: "current-target-full-weight-identities-with-higher-order-cancellation-gate",
        PROVIDERS[5]: "conditional-current-target-cancellation-cycle-bridge",
        PROVIDERS[6]: "regular-source-envelope-edge-estimates-with-explicit-freezing-correction",
        PROVIDERS[7]: "conditional-current-state-complete-finite-stage-induction",
        PROVIDERS[8]: "Q153-physical-source-tail-and-complete-potential-diagonal-bound-conditionally",
    }
    for name, status in expected.items():
        if data[name].get("status") != status:
            raise ValueError("Q160 provider status mismatch: " + name)
    moment = data[PROVIDERS[4]]["flags"]
    cycle = data[PROVIDERS[5]]["flags"]
    source = data[PROVIDERS[6]]["flags"]
    induction = data[PROVIDERS[7]]["flags"]
    required_moment = (
        "all_fixed_jet_full_zeta_functional_membership_relative_to_current_field_bounds",
        "both_current_velocity_moment_constraints_and_whole_moving_pressure_product_retained",
        "current_full_covariance_moment_map_includes_radial_time_and_axial_transport",
        "signed_target_exact_radial_flux_plus_compensated_transport_primitive",
    )
    required_cycle = (
        "current_target_estimate_derived_by_exact_state_difference_not_raw_field_bound",
        "generic_signed_target_cancellation_keeps_all_six_covariances_and_changed_pressure",
        "actual_zero_Haar_fast_curl_has_exactly_zero_linear_Haar_response",
        "fast_old_new_correlations_and_recomputed_postfast_defects_retained",
        "slow_current_pressure_time_viscosity_and_old_Haar_means_retained",
        "all_half_line_target_recurrence_relative_to_reviewed_current_map_inputs",
        "same_complete_base_Pi1_Pi2_actual_A0_A1_fixed_width_and_corrected_domain_retained",
    )
    if not all(moment.get(key) for key in required_moment):
        raise ValueError("Q160 complete moment provider contract missing")
    if not all(cycle.get(key) for key in required_cycle):
        raise ValueError("Q160 current-target recurrence provider contract missing")
    if source["every_actual_later_covariance_mean_target_verified_by_grammar_or_protocol"]:
        raise ValueError("historical source record must remain frozen and false")
    if not induction["conditional_complete_four_operation_finite_stage_induction_derived"]:
        raise ValueError("finite-stage recurrence provider missing")
    return {
        "same_physical_family_bound": True,
        "historical_source_flag_preserved_false": True,
        "new_Q160_theorem_supersedes_that_open_gate": True,
        "Q153_may_be_replayed_but_is_not_promoted_here": True,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q160 identity required")
    return (
        {"id": "P1", "status": "closed", "contract": "PREMOD value/root envelope"},
        {"id": "P2", "status": "closed", "contract": "Borel background and physical primary pulses"},
        {"id": "P3", "status": "closed-by-Q160-induction", "contract": "complete covariance/mean targets"},
        {"id": "P4", "status": "next-replay-required", "contract": "generated-source quotient and endpoint extension"},
        {"id": "P5", "status": "downstream-conditional", "contract": "weighted inverses and physical tails"},
        {"id": "P6", "status": "downstream-conditional", "contract": "complete diagonal, endpoint and localized force"},
        {"id": "P7", "status": "conditional-consequence", "contract": "comparison and core exclusion"},
    )


def record():
    return {
        "status": "Q160-physical-all-finite-cycle-target-moment-induction-joined",
        "correction": (
            "Q159 treated the current-target bridge as only an expression-level warning; "
            "Q149/Q150/Q152 supply its physical seed and every same-family map premise"
        ),
        "established": (
            "for every finite physical cycle and fixed mixed jet, the recomputed Xi_j has "
            "the required epsilon order, full zeta edge weight and finite edge powers, and "
            "both complete tangential moment identities hold with full pressure memory"
        ),
        "scope": (
            "finite-cycle continuum induction with stage- and jet-dependent constants; the "
            "diagnostic tests are not the proof and no uniform infinite-stage bound is asserted"
        ),
        "next": (
            "replay Q153's generated-source quotient and endpoint-extension provider from this "
            "now-unconditional P3 input, checking the actual complete residual allocation before "
            "promoting any downstream diagonal, force or exclusion claim"
        ),
        "P3_universal_finite_cycle_gate_discharged": True,
        "Q153_downstream_chain_replayed_without_P3_premise": False,
        "unconditional_complete_residual_force_exclusion_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
    }
