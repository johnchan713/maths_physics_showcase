"""Reconcile the physical target and generated-source cycle by simultaneous induction."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
IDENTITY = "Q180-PHYSICAL-TARGET-SOURCE-SIMULTANEOUS-INDUCTION"
PROVIDERS = (
    "core_current_head_q179_covariance_mean_cycle_frontier",
    "core_current_head_q160_physical_target_induction",
    "core_current_head_q161_source_tail_simultaneous_induction",
    "core_current_head_q174_q5_operator_source_audit",
    "current_target_cycle_bridge",
    "finite_stage_induction",
    "source_envelope_edge_audit",
    "regular_source_foundation",
    "sourced_inverse_foundation",
    "compact_mean_foundation",
)
EXPECTED = {
    PROVIDERS[0]: "Q179-initial-and-first-feedback-covariance-means-closed-P3-universal-target-obstructed",
    PROVIDERS[1]: "Q160-physical-all-finite-cycle-target-moment-induction-joined",
    PROVIDERS[2]: "Q161-physical-source-inverse-tail-simultaneous-induction-joined",
    PROVIDERS[3]: "Q174-Q5-operator-source-audited-with-family-and-typed-domain-boundary",
    PROVIDERS[4]: "conditional-current-target-cancellation-cycle-bridge",
    PROVIDERS[5]: "conditional-current-state-complete-finite-stage-induction",
    PROVIDERS[6]: "regular-source-envelope-edge-estimates-with-explicit-freezing-correction",
    PROVIDERS[7]: "regular-source-foundation-sharper-linear-gain-and-actual-initial-supported-input-relative-to-core-and-pulse-lemmas",
    PROVIDERS[8]: "explicit-same-domain-sourced-inverse-and-high-harmonic-parameter-jets-relative-to-reviewed-frame-inputs",
    PROVIDERS[9]: "compact-mean-inverse-foundation-sharp-axial-conditioning-relative-to-core-and-primary-inputs",
}


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


q160 = load("q180_q160", "core_current_head_q160_physical_target_induction/q160.py")
q161 = load("q180_q161", "core_current_head_q161_source_tail_simultaneous_induction/q161.py")


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q180 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def rates(j=0, kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q180 identity required")
    if not isinstance(j, int) or isinstance(j, bool) or j < 0:
        raise ValueError("finite nonnegative cycle required")
    kappa = Q(kappa)
    if not 0 < kappa < Q(1, 100):
        raise ValueError("reviewed kappa range required")
    delta = Q(1, 2) - 3 * kappa
    B0 = Q(3, 2) - 2 * kappa
    B = B0 + j * delta
    C = B + Q(1, 2)
    return {
        "j": j,
        "kappa": kappa,
        "delta": delta,
        "B0": B0,
        "B_j": B,
        "C_j": C,
        "B_next": B + delta,
        "C_next": C + delta,
        "target_bridge_order": C + Q(1, 2) - kappa,
        "target_slack": 2 * kappa,
    }


def coupled_step(j=0, kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q180 identity required")
    r = rates(j, kappa)
    source = q161.coupled_step(j, kappa)
    target = q160.target_step(r["C_j"], kappa)
    source_margins = {
        name: row["margin_over_B_next"] for name, row in source["source_rows"].items()
    }
    target_margins = {
        name: row["margin_over_required_next"] for name, row in target["rows"].items()
    }
    return {
        "rates": r,
        "source_margins": source_margins,
        "target_margins": target_margins,
        "source_limiter": min(source_margins, key=source_margins.get),
        "source_minimum_margin": min(source_margins.values()),
        "target_limiter": min(target_margins, key=target_margins.get),
        "target_minimum_margin": min(target_margins.values()),
        "target_slack": r["target_slack"],
        "closes_source_next": source["all_regular_rows_reach_B_next"],
        "closes_target_next": target["all_rows_close_physical_next_target"],
        "closes_coupled_next": (
            source["all_regular_rows_reach_B_next"]
            and target["all_rows_close_physical_next_target"]
            and r["target_bridge_order"] >= r["C_next"]
        ),
    }


def moment_and_edge_contract(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q180 identity required")
    return {
        "target": "|D^I Xi_j| <= C_(j,I) epsilon^(C_j) zeta ya^(-p_(j,I)) yb^(-q_(j,I))",
        "angular": "int R^2 overline(E_theta,j)=epsilon partial_Z J_theta,j",
        "axial": "int R overline(E_z,j)=epsilon partial_Z(J_z,j+c_rho P_mean,j)",
        "memory": (
            "all six covariance entries, axial velocity, Pi1/Pi2, compact pressure, "
            "the whole moving c_rho P_mean product, zero-Haar data, both velocity moments, "
            "reserved interior hull and unchanged C^-20 width"
        ),
        "wave_weight": "sqrt(zeta) times sqrt(zeta) gives zeta",
        "mean_weight": "zeta times zeta is bounded by zeta because 0<zeta<=1",
        "primitive_gate": "exact zero weighted total is required at both edges",
        "compactness_alone_suffices": False,
        "stage_uniform_constant_claimed": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q180 identity required")
    status = statuses()
    q179 = evidence(PROVIDERS[0])["symbolic"]
    q160_record = evidence(PROVIDERS[1])["symbolic"]
    q161_record = evidence(PROVIDERS[2])["symbolic"]
    q174_record = evidence(PROVIDERS[3])["review"]
    source = evidence(PROVIDERS[6])["flags"]
    regular = evidence(PROVIDERS[7])["flags"]
    inverse = evidence(PROVIDERS[8])["flags"]
    compact = evidence(PROVIDERS[9])["flags"]
    step = coupled_step()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "Q179_supplies_physical_seed_but_leaves_universal_gate_open": (
            not q179["unconditional_complete_covariance_cycle_residual_force_or_blowup_certified"]
        ),
        "Q160_has_exact_target_recurrence": q160_record["P3_universal_finite_cycle_gate_discharged"],
        "Q161_corrects_sequential_order_to_simultaneous_induction": "mutually dependent" in q161_record["correction"],
        "physical_family_not_promoted_to_canonical_family": (
            q174_record["checks"]["physical_join_preserves_nonidentification"]
            and not q174_record["scope"]["physical_source_theorem_promoted_to_canonical_source_theorem"]
        ),
        "source_grammar_retains_full_path_and_current_targets": source["generated_regular_source_grammar_has_full_path_and_entrance_inputs_relative_to_current_target_estimates"],
        "endpoint_and_complete_operator_rows_retained": (
            regular["full_source_endpoint_extension_and_closed_label_supports_preserved"]
            and regular["full_viscosity_pressure_graph_coefficients_and_cylindrical_connections_retained"]
        ),
        "same_domain_sourced_inverse_retained": (
            inverse["explicit_finite_sourced_ODE_existence_uniqueness_from_reviewed_frame"]
            and inverse["same_accepted_domain_for_all_nonzero_integer_modes"]
        ),
        "compact_mean_inverses_retain_all_fixed_orders": compact["all_fixed_order_compact_and_temporal_mean_inverse_estimates_derived"],
        "coupled_step_closes_with_exact_positive_target_slack": step["closes_coupled_next"] and step["target_slack"] == Q(1, 50000),
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "coupled_step": step,
        "moment_and_edge_contract": moment_and_edge_contract(),
        "correction": {
            "Q160_standalone_sequential_P3_then_P4_reading_is_invalid": True,
            "Q160_target_recurrence_is_used_inside_Q161_simultaneous_step": True,
            "Q179_obstruction_is_closed_only_relative_to_named_continuum_maps": True,
            "historical_frozen_evidence_modified": False,
        },
        "established": {
            "physical_P3_target_and_both_moments_for_each_finite_cycle": True,
            "physical_P4_generated_source_quotient_and_endpoint_extension": True,
            "physical_P5_fixed_stage_tail_for_each_fixed_jet_and_power": True,
            "constants_uniform_in_cycle_harmonic_jet_or_power": False,
        },
        "scope": {
            "same_inherited_physical_family": True,
            "canonical_family_identified": False,
            "complete_diagonal_endpoint_force_replayed_here": False,
            "unconditional_complete_residual_force_exclusion_or_blowup": False,
            "tests_are_continuum_proof": False,
        },
        "next": (
            "replay the complete-potential diagonal, punctured endpoint and localized-force "
            "ledger on this coupled physical state, retaining every cutoff curl, pressure "
            "gradient, viscosity term, nonlinear cross and square"
        ),
    }


def record():
    return {
        "status": "Q180-physical-target-source-tail-simultaneous-induction-reconciled",
        "established": (
            "relative to the named current-state source and inverse estimates, one simultaneous "
            "finite-cycle induction proves the physical zeta-weighted Xi_j estimate, both exact "
            "moments, generated-source quotient/endpoint extension and every fixed-stage tail"
        ),
        "correction": (
            "Q160 cannot be read as a complete sequential P3 proof followed by P4; the next "
            "target and next generated source are mutually dependent and close only in the "
            "Q161 simultaneous induction"
        ),
        "limiting_inequalities": (
            "source gain delta=1/2-3*kappa is saturated by the sourced and signed full linear "
            "rows; target gain 1/2-kappa exceeds the required delta by exactly 2*kappa"
        ),
        "historical_frozen_records_modified": False,
        "unconditional_complete_residual_force_exclusion_or_blowup_certified": False,
    }
