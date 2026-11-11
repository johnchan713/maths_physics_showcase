"""Reconcile Q172--Q182 and expose the first undischarged physical estimate."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import json
import math
import sys

HERE = Path(__file__).resolve().parent
IDENTITY = "Q183-CURRENT-HEAD-PROVIDER-DISCHARGE"
PROVIDERS = (
    "core_current_head_q26_source_domain_join",
    "core_current_head_q161_source_tail_simultaneous_induction",
    "core_current_head_q172_q2_q3_source_audit",
    "core_current_head_q173_q4_source_audit",
    "core_current_head_q175_canonical_source_route_reconciliation",
    "core_current_head_q176_modulation_repair_dependency_audit",
    "core_current_head_q177_weighted_support_source_reconciliation",
    "core_current_head_q178_physical_primary_source_audit",
    "core_current_head_q180_target_source_simultaneous_induction",
    "core_current_head_q181_diagonal_endpoint_force_replay",
    "core_current_head_q182_comparison_core_exclusion_replay",
)
EXPECTED = {
    PROVIDERS[0]: "fixed-q-local-source-domain-closure",
    PROVIDERS[1]: "Q161-physical-source-inverse-tail-simultaneous-induction-joined",
    PROVIDERS[2]: "Q172-Q2-source-audited-and-Q3-uniform-frame-corrected",
    PROVIDERS[3]: "Q173-Q4-fixed-profile-source-audited-with-generated-source-boundary",
    PROVIDERS[4]: "Q175-canonical-typed-source-route-reconciled-and-independent-frontier-restored",
    PROVIDERS[5]: "Q176-compact-modulation-exact-repair-audited-with-refined-dependency-order",
    PROVIDERS[6]: "Q177-selected-input-weighted-stress-and-reserved-support-route-reconciled",
    PROVIDERS[7]: "Q178-physical-primary-pulse-source-boundary-and-covariance-obstruction-audited",
    PROVIDERS[8]: "Q180-physical-target-source-tail-simultaneous-induction-reconciled",
    PROVIDERS[9]: "Q181-physical-diagonal-endpoint-localized-force-replayed-conditionally",
    PROVIDERS[10]: "Q182-physical-pressure-comparison-protected-core-replayed-conditionally",
}


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


q153 = load(
    "q183_q153",
    "core_current_head_q153_physical_source_tail_diagonal_join/q153.py",
)


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q183 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def bell_majorant(raw_source_jets, log_gaussian_jets, identity=IDENTITY):
    """Exact Q153 Bell/Leibniz majorant, without manufacturing its inputs."""
    if identity != IDENTITY:
        raise ValueError("Q183 identity required")
    if raw_source_jets is None or log_gaussian_jets is None:
        return {
            "input_suprema_supplied": False,
            "inverse_gaussian": None,
            "quotient": None,
            "conditional_only": True,
        }
    if not raw_source_jets or len(raw_source_jets) != len(log_gaussian_jets):
        raise ValueError("equal nonempty jet arrays required")
    source = [Q(value) for value in raw_source_jets]
    log_jets = [Q(value) for value in log_gaussian_jets]
    if any(value < 0 for value in source + log_jets):
        raise ValueError("nonnegative majorants required")
    order = len(source) - 1
    inverse = [Q(1)] + [Q(0)] * order
    for n in range(1, order + 1):
        inverse[n] = sum(
            Q(math.comb(n - 1, i - 1)) * log_jets[i] * inverse[n - i]
            for i in range(1, n + 1)
        )
    quotient = [
        sum(
            Q(math.comb(n, i)) * inverse[i] * source[n - i]
            for i in range(n + 1)
        )
        for n in range(order + 1)
    ]
    inherited = q153.bell_quotient(source, log_jets)
    if inherited != {"inverse_gaussian": inverse, "quotient": quotient}:
        raise ValueError("Bell recurrence disagrees with frozen Q153")
    return {
        "input_suprema_supplied": True,
        "inverse_gaussian": inverse,
        "quotient": quotient,
        "conditional_only": False,
    }


def physical_generated_source_gate(
    stage=0,
    harmonic=1,
    jet_order=2,
    raw_source_jets=None,
    log_gaussian_jets=None,
    endpoint_flatness_proved=False,
    identity=IDENTITY,
):
    if identity != IDENTITY:
        raise ValueError("Q183 identity required")
    if any(isinstance(value, bool) or not isinstance(value, int) for value in (stage, harmonic, jet_order)):
        raise ValueError("integer indices required")
    if stage < 0 or harmonic == 0 or jet_order < 0:
        raise ValueError("finite stage, nonzero harmonic and nonnegative jet required")
    bell = bell_majorant(raw_source_jets, log_gaussian_jets)
    arrays_cover_jet = (
        bell["input_suprema_supplied"]
        and len(raw_source_jets) == jet_order + 1
    )
    return {
        "stage": stage,
        "harmonic": harmonic,
        "jet_order": jet_order,
        "required_raw_inputs": (
            "finite same-domain weighted suprema A_i for every raw physical generated-source "
            "jet after extracting P, 0<=i<=jet_order"
        ),
        "required_gaussian_inputs": (
            "finite same-domain weighted suprema K_i for D^i log(P), 1<=i<=jet_order"
        ),
        "required_endpoint_input": (
            "all mixed endpoint jets of the actual cut source vanish after pressure, viscosity, "
            "cutoff and nonlinear interaction terms are assembled"
        ),
        "target_inequality": (
            "sup_(ell>=ell_*) epsilon_ell^(-B_j) sup_(Omega_(j,ell)) "
            "zeta^(-1/2) y_a^r y_b^s |D^I(f_(j,m)/P)| < infinity"
        ),
        "bell": bell,
        "arrays_cover_requested_jet": arrays_cover_jet,
        "endpoint_flatness_proved": bool(endpoint_flatness_proved),
        "gate_discharged": arrays_cover_jet and bool(endpoint_flatness_proved),
    }


def frontier_reconciliation(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q183 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    q173 = data[PROVIDERS[3]]["review"]
    q175 = data[PROVIDERS[4]]["review"]
    q177 = data[PROVIDERS[6]]["review"]
    q178 = data[PROVIDERS[7]]["review"]
    q180 = data[PROVIDERS[8]]["review"]
    return {
        "Q171_broad_Q2_Q3_frontier_is_stale": True,
        "Q172_Q177_supply_current_conditional_routes": True,
        "Q26_closes_current_local_later_tail_domain": True,
        "historical_finite_head_determinant_margin_is_not_first_current_obstruction": True,
        "Q173_locates_generated_source_boundary": not q173["first_open_gate"]["fixture"]["gate_closed"],
        "Q175_route_is_not_independent_validation": not q175["scope"]["provider_theorems_independently_reproved_here"],
        "Q177_selected_input_route_is_conditional": not q177["scope"]["provider_theorems_replaced_by_tests_or_numerical_diagnostics"],
        "Q178_primary_bounds_do_not_prove_generated_sources": not q178["scope"]["all_generated_nonlinear_source_inputs_verified"],
        "Q180_reuses_named_source_estimates": q180["correction"]["Q179_obstruction_is_closed_only_relative_to_named_continuum_maps"],
        "Q181_Q182_are_downstream_conditional_replays": True,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q183 identity required")
    status = statuses()
    frontier = frontier_reconciliation()
    missing = physical_generated_source_gate()
    control = physical_generated_source_gate(
        raw_source_jets=(1, 2, 3),
        log_gaussian_jets=(0, 3, 5),
        endpoint_flatness_proved=True,
    )
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "all_frontier_reconciliations_hold": all(frontier.values()),
        "bell_recurrence_does_not_supply_missing_suprema": (
            missing["bell"]["conditional_only"] and not missing["gate_discharged"]
        ),
        "finite_inputs_and_endpoint_flatness_would_close_fixed_jet_gate": control["gate_discharged"],
        "downstream_pressure_comparison_remains_conditional": not evidence(PROVIDERS[10])["symbolic"][
            "unconditional_admissible_force_breakdown_or_blowup_certified"
        ],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "frontier_reconciliation": frontier,
        "blocking_gate": missing,
        "finite_control": control,
        "correction": {
            "Q161_Bell_arithmetic_is_valid": True,
            "Q161_synthetic_source_quotient_fixture_is_not_a_physical_supremum_proof": True,
            "Q161_provider_flags_route_but_do_not_independently_derive_Ai_Ki": True,
            "historical_frozen_evidence_modified": False,
        },
        "scope": {
            "current_head_provider_frontier_reconciled": True,
            "actual_physical_generated_source_supremum_proved": False,
            "endpoint_extension_for_complete_actual_source_proved_here": False,
            "complete_residual_improvement_proved": False,
            "unconditional_pulse_force_exclusion_or_blowup": False,
            "tests_are_continuum_proof": False,
        },
        "next": (
            "derive the fixed-cycle physical A_i and K_i suprema and actual endpoint flatness "
            "directly from the explicit current field and complete source operator, retaining "
            "pressure memory, viscosity, cutoffs and every interaction; then propagate the result "
            "through the finite-cycle simultaneous induction"
        ),
    }


def record():
    return {
        "status": "Q183-current-head-provider-frontier-reconciled-to-generated-source-supremum",
        "established": (
            "Q172--Q182 supersede Q171's broad frontier and conditionally route every later "
            "modulation, repair, covariance, induction, force and pressure step; the exact Bell "
            "recurrence converts supplied raw-source and log-Gaussian jet suprema into quotient "
            "bounds at each fixed physical jet"
        ),
        "correction": (
            "the Q161 executable uses synthetic finite A_i and K_i arrays and its provider flags "
            "only route the corresponding continuum estimates; neither independently derives "
            "those actual physical suprema or complete endpoint flatness"
        ),
        "obstruction": (
            "prove the all-fixed-jet weighted generated-source inequality and smooth endpoint "
            "extension from the explicit physical field/source formulas, with full pressure, "
            "viscosity, cutoff and interaction memory"
        ),
        "historical_frozen_records_modified": False,
        "unconditional_actual_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
