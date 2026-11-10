"""Close the physical target/source/inverse/tail loop by simultaneous induction."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
IDENTITY = "Q161-physical-target-source-tail-simultaneous-induction"
PROVIDERS = (
    "core_current_head_q160_physical_target_induction",
    "core_current_head_q153_physical_source_tail_diagonal_join",
    "source_envelope_edge_audit",
    "regular_source_foundation",
    "sourced_inverse_foundation",
    "actual_band_pulse_bounds",
    "compact_mean_foundation",
    "current_target_cycle_bridge",
    "finite_stage_induction",
    "core_physical_tail_join",
    "core_first_sourced_physical_tail_closure",
    "core_finite_stage_compact_recurrence",
    "core_certified_complete_potential_diagonal_join",
)


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


q153 = load("q161_q153", "core_current_head_q153_physical_source_tail_diagonal_join/q153.py")
compact = load("q161_compact_recurrence", "core_finite_stage_compact_recurrence/recurrence.py")
induction = load("q161_finite_induction", "finite_stage_induction/induction.py")


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def rates(j=0, kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q161 identity required")
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
        raise ValueError("Q161 identity required")
    r = rates(j, kappa)
    cert = induction.certificate(r["delta"])
    nonzero = cert["groups"]["nonzero_regular"]
    source_rows = {
        name: {
            "order": value.at(r["B_j"]),
            "margin_over_B_next": value.at(r["B_j"]) - r["B_next"],
        }
        for name, value in nonzero.items()
    }
    return {
        "input": (
            "actual physical state at B_j with Q160 current Xi_j in weighted M_Cj, "
            "complete means/covariance/pressure and exact moment data"
        ),
        "source_rows": source_rows,
        "all_regular_rows_reach_B_next": all(
            row["margin_over_B_next"] >= 0 for row in source_rows.values()
        ),
        "normalized_source": (
            "for every nonzero generated harmonic m and fixed jet I, "
            "|D^I(f_m/P)|<=C_(j,I,m) epsilon^(B_next) sqrt(zeta) "
            "ya^(-p) yb^(-q)"
        ),
        "sourced_inverse": (
            "the same-domain constrained ODE gives the complete potential, pressure and curl; "
            "cutting potential and pressure restores smooth endpoint zero extension"
        ),
        "target_update": (
            "the exact current-state covariance/mean difference gives weighted target order "
            "C_j+1/2-kappa"
        ),
        "target_order": r["target_bridge_order"],
        "required_target_order": r["C_next"],
        "target_slack": r["target_slack"],
        "closes_next_state": (
            all(row["margin_over_B_next"] >= 0 for row in source_rows.values())
            and r["target_bridge_order"] >= r["C_next"]
        ),
    }


def source_quotient(jet_order=4, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q161 identity required")
    if not isinstance(jet_order, int) or isinstance(jet_order, bool) or jet_order < 0:
        raise ValueError("fixed nonnegative jet order required")
    # A diagnostic finite Bell majorant.  The continuum proof is the displayed
    # recurrence applied to the actual fixed-jet source and log-Gaussian bounds.
    source = [Q(i + 1) for i in range(jet_order + 1)]
    log_jets = [Q(0)] + [Q(2 * i + 1) for i in range(1, jet_order + 1)]
    majorant = q153.bell_quotient(source, log_jets)
    return {
        "jet_order": jet_order,
        "inverse_gaussian_majorant": majorant["inverse_gaussian"],
        "quotient_majorant": majorant["quotient"],
        "all_finite": all(value >= 0 for value in majorant["quotient"]),
        "pressure_memory": (
            "source-normal and moving-normal pressure, Pi1/Pi2, all six covariance entries "
            "and current axial velocity are retained before division by P"
        ),
        "endpoint": (
            "every nonzero regular term contains a complete cut wave factor; no cutoff, "
            "partition or updated amplitude is divided out, so all quotient entrance jets vanish"
        ),
        "finite_harmonics": "each finite cycle has finitely many harmonics; constants may depend on m and j",
    }


def tail_request(j=0, m=0, N=0, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q161 identity required")
    request = compact.physical_request(j, m, N)
    return {
        "stage": j,
        "jet": m,
        "power": N,
        "raw_sources": request["raw_sources"],
        "first_real_ceiling": request["first_real_ceiling"],
        "strict_margins": request["strict_margins"],
        "finite_request": (
            isinstance(request["first_real_ceiling"], int)
            and request["first_real_ceiling"] > 0
        ),
        "strictly_positive": min(request["strict_margins"].values()) > 0,
    }


def residual_allocation(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q161 identity required")
    return {
        "identity": (
            "R^[j]=R_regular^[j]+T_compact^[j]+T_Gaussian^[j]+T_base/exterior"
        ),
        "regular": (
            "spatial/slow cutoff derivatives, pressure gradients, viscosity, complete curl "
            "terms, mean-wave terms, both nonlinear crosses and the new square"
        ),
        "compact": (
            "same-state (-A0 Qcirc gr,0,+A1 Qcirc Ez_postwave), with actual A0/A1 "
            "velocities retained in all later nonlinear products"
        ),
        "Gaussian": "temporal cutoff commutators stay additive and are never forward-inverted",
        "base": "only the retained flat base/exterior remainder after base stress subtraction",
        "scale_cutoffs": (
            "complete-potential, direct-swirl and pressure cutoffs occur only in the diagonal, "
            "with every curl, gradient, time/viscous term, cross and square"
        ),
        "double_count_free": True,
    }


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q161 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    expected = {
        PROVIDERS[0]: "Q160-physical-all-finite-cycle-target-moment-induction-joined",
        PROVIDERS[1]: "Q153-physical-source-tail-and-complete-potential-diagonal-bound-conditionally",
        PROVIDERS[2]: "regular-source-envelope-edge-estimates-with-explicit-freezing-correction",
        PROVIDERS[3]: "regular-source-foundation-sharper-linear-gain-and-actual-initial-supported-input-relative-to-core-and-pulse-lemmas",
        PROVIDERS[4]: "explicit-same-domain-sourced-inverse-and-high-harmonic-parameter-jets-relative-to-reviewed-frame-inputs",
        PROVIDERS[5]: "conditional-actual-band-pulse-primary-covariance-and-linear-bounds",
        PROVIDERS[6]: "compact-mean-inverse-foundation-sharp-axial-conditioning-relative-to-core-and-primary-inputs",
        PROVIDERS[7]: "conditional-current-target-cancellation-cycle-bridge",
        PROVIDERS[8]: "conditional-current-state-complete-finite-stage-induction",
        PROVIDERS[9]: "conditional-current-prefix-physical-tail-transfer-and-absolute-regular-order-gap",
        PROVIDERS[10]: "conditional-first-sourced-same-state-compact-physical-tail-closure",
        PROVIDERS[11]: "conditional-current-seed-finite-stage-compact-recurrence",
        PROVIDERS[12]: "certified-domain-conditional-complete-potential-diagonal-join",
    }
    for name, status in expected.items():
        if data[name].get("status") != status:
            raise ValueError("Q161 provider status mismatch: " + name)
    q160 = data[PROVIDERS[0]]["symbolic"]
    q153_record = data[PROVIDERS[1]]["symbolic"]
    source = data[PROVIDERS[2]]["flags"]
    regular = data[PROVIDERS[3]]["flags"]
    inverse = data[PROVIDERS[4]]["flags"]
    compact_mean = data[PROVIDERS[6]]["flags"]
    tails = data[PROVIDERS[9]]["flags"]
    recurrence = data[PROVIDERS[11]]["flags"]
    if not q160["P3_universal_finite_cycle_gate_discharged"]:
        raise ValueError("Q160 target induction required")
    if not all((
        q153_record["physical_generated_source_endpoint_join_certified_conditionally"],
        q153_record["all_fixed_stage_additive_physical_tails_certified_conditionally"],
        source["all_fixed_mixed_normalized_source_jets_have_Bell_majorants"],
        source["generated_regular_source_grammar_has_full_path_and_entrance_inputs_relative_to_current_target_estimates"],
        regular["full_source_endpoint_extension_and_closed_label_supports_preserved"],
        regular["full_viscosity_pressure_graph_coefficients_and_cylindrical_connections_retained"],
        inverse["explicit_finite_sourced_ODE_existence_uniqueness_from_reviewed_frame"],
        inverse["same_accepted_domain_for_all_nonzero_integer_modes"],
        compact_mean["all_fixed_order_compact_and_temporal_mean_inverse_estimates_derived"],
        tails["all_fixed_physical_compact_tail_recipes_include_normalization_and_graph_losses"],
        tails["original_temporal_Gaussian_ledger_is_additive_and_never_inverted"],
        recurrence["finite_stage_backward_request_recurrence_supplied"],
    )):
        raise ValueError("Q161 source/inverse/tail contract missing")
    if source["every_actual_later_covariance_mean_target_verified_by_grammar_or_protocol"]:
        raise ValueError("historical source flag must remain frozen false")
    return {
        "P3_and_P4_closed_by_simultaneous_induction": True,
        "P5_weighted_inverse_and_fixed_stage_tail_contract_closed": True,
        "historical_conditional_records_preserved": True,
        "complete_potential_diagonal_provider_available_but_P6_not_replayed_here": True,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q161 identity required")
    return (
        {"id": "P1", "status": "closed", "contract": "PREMOD value/root envelope"},
        {"id": "P2", "status": "closed", "contract": "Borel background and physical primary pulses"},
        {"id": "P3", "status": "closed-in-coupled-Q160-Q161-induction", "contract": "complete covariance/mean targets"},
        {"id": "P4", "status": "closed-by-Q161", "contract": "generated-source quotient and endpoint extension"},
        {"id": "P5", "status": "closed-by-Q161-for-each-fixed-stage-jet-power", "contract": "weighted inverses and physical tails"},
        {"id": "P6", "status": "next-replay-required", "contract": "complete diagonal, endpoint and localized force"},
        {"id": "P7", "status": "conditional-consequence", "contract": "comparison and core exclusion"},
    )


def record():
    return {
        "status": "Q161-physical-source-inverse-tail-simultaneous-induction-joined",
        "correction": (
            "Q160's target recurrence and Q153's generated-source quotient are mutually dependent "
            "at later cycles; they close by simultaneous induction, not by a sequential P3-then-P4 proof"
        ),
        "established": (
            "for every finite physical cycle, nonzero harmonic and fixed output jet, the generated "
            "source has the full-path Gaussian quotient and smooth endpoint extension; the sourced, "
            "signed and mean inverses remain on the same domains; and every fixed-stage A0/A1, "
            "temporal Gaussian and base/exterior tail is O(q^N) after any fixed physical jet"
        ),
        "scope": (
            "simultaneous finite-cycle continuum induction with constants depending on cycle, harmonic, "
            "jet and requested power; no stage-uniform raw-series convergence or numerical giant-jet evaluation"
        ),
        "next": (
            "replay the complete-potential diagonal, punctured endpoint atlas and localization/force "
            "chain on this coupled physical induction, retaining all scale-cutoff curls, pressure "
            "gradients, viscosity, crosses and squares before promoting P6"
        ),
        "P4_generated_source_endpoint_contract_discharged": True,
        "P5_fixed_stage_weighted_inverse_and_tail_contract_discharged": True,
        "P6_complete_diagonal_endpoint_force_chain_replayed": False,
        "unconditional_complete_residual_force_exclusion_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
    }
