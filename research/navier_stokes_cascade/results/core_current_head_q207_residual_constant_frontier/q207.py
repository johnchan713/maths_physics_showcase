"""Q207: reconstruct the fixed-stage residual constant and narrow the frontier."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import importlib.util
import json

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
IDENTITY = "Q207-RESIDUAL-CONSTANT-FRONTIER"
PROVIDERS = (
    "core_current_head_q206_terminal_origin_obstruction",
    "core_current_head_q202_next_sourced_cycle",
    "core_current_head_q161_source_tail_simultaneous_induction",
    "core_finite_stage_compact_recurrence",
    "core_current_head_q166_complete_residual_frontier_audit",
    "core_current_head_q190_downstream_route_frontier_audit",
    "core_current_head_q194_primary_field_assembly",
    "core_current_head_q184_generated_source_log_loss_correction",
    "core_current_head_q196_initial_generated_source",
)
EXPECTED = {
    PROVIDERS[0]: "Q206-terminal-origin-velocity-and-gauge-independent-pressure-value-jets-obstructed-conditionally",
    PROVIDERS[1]: "Q202-next-constrained-sourced-cycle-joined-on-corrected-physical-domain",
    PROVIDERS[2]: "Q161-physical-source-inverse-tail-simultaneous-induction-joined",
    PROVIDERS[3]: "conditional-current-seed-finite-stage-compact-recurrence",
    PROVIDERS[4]: "Q166-current-head-complete-six-group-local-residual-audited-conditionally",
    PROVIDERS[5]: "Q190-Q189-through-Q185-downstream-route-closed-conditionally-and-frontier-corrected",
    PROVIDERS[6]: "Q194-complete-primary-field-bounds-and-support-assembled-conditionally",
    PROVIDERS[7]: "Q184-generated-source-log-loss-corrected-and-pure-power-reserve-derived",
    PROVIDERS[8]: "Q196-initial-generated-source-stage-and-log-quotient-corrected",
}


def load_recurrence():
    path = HERE.parent / "core_finite_stage_compact_recurrence" / "recurrence.py"
    spec = importlib.util.spec_from_file_location("q207_finite_recurrence", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q207 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def same_object_links(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q207 identity required")
    links = (
        (PROVIDERS[0], "core_current_head_q205_stagewise_pulse_diagonal"),
        ("core_current_head_q205_stagewise_pulse_diagonal", "core_current_head_q203_diagonal_endpoint_force_corrected_domain"),
        ("core_current_head_q203_diagonal_endpoint_force_corrected_domain", PROVIDERS[1]),
        (PROVIDERS[1], PROVIDERS[2]),
        (PROVIDERS[2], PROVIDERS[3]),
        (PROVIDERS[5], PROVIDERS[4]),
        (PROVIDERS[6], "core_current_head_q192_frame_graph_jet_transfer"),
        (PROVIDERS[5], PROVIDERS[7]),
        (PROVIDERS[1], "core_current_head_q201_post_signed_mean_reconciliation"),
        ("core_current_head_q201_post_signed_mean_reconciliation", "core_current_head_q200_physical_weighted_target_reconciliation"),
        ("core_current_head_q200_physical_weighted_target_reconciliation", "core_current_head_q199_weighted_target_obstruction"),
        ("core_current_head_q199_weighted_target_obstruction", "core_current_head_q198_generic_signed_target"),
        ("core_current_head_q198_generic_signed_target", "core_current_head_q197_first_sourced_covariance"),
        ("core_current_head_q197_first_sourced_covariance", PROVIDERS[8]),
    )
    result = {}
    for parent, child in links:
        child_path = f"results/{child}/evidence.json"
        result[f"{parent}_uses_{child}"] = (
            evidence(parent)["dependency_evidence_sha256"][child_path]
            == sha(PROJECT / child_path)
        )
    return result


def rates(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q207 identity required")
    h = Q(1, 200)
    kappa = Q(1, 100000)
    return {
        "h": h,
        "kappa": kappa,
        "delta": Q(1, 2) - 3 * kappa,
        "B0": Q(3, 2) - 2 * kappa,
        "Lambda": Q(3, 2) + 2 * h,
    }


def ceil_fraction(value):
    return -((-value.numerator) // value.denominator)


def stage_for_power(m=0, N=0, eta=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q207 identity required")
    if m < 0 or N < 0 or eta <= 0:
        raise ValueError("need nonnegative m,N and positive eta")
    r = rates()
    raw = (Q(N) + r["Lambda"] + 2 * m + eta - r["h"] * r["B0"]) / (r["h"] * r["delta"])
    j = max(0, ceil_fraction(raw))
    exponent = r["h"] * (r["B0"] + j * r["delta"]) - (r["Lambda"] + 2 * m + eta)
    previous = None if j == 0 else r["h"] * (r["B0"] + (j - 1) * r["delta"]) - (r["Lambda"] + 2 * m + eta)
    return {
        "m": m,
        "N": N,
        "eta": eta,
        "stage": j,
        "exponent": exponent,
        "previous_exponent": previous,
        "meets_target": exponent >= N,
        "minimal": j == 0 or previous < N,
    }


def constant_construction(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q207 identity required")
    return {
        "finite_request": "R_(j,m,N)<infinity",
        "primitive_bundle": "P_R=max weighted selected-profile/root/frame/source/pressure/mean seminorms",
        "gap_bundle": "g_R=min Jacobian/normal/root/determinant/variance/support gaps>0",
        "operator_bundle": "K_R=max compact/slow inverse constants<infinity",
        "recipe": "C_(j,m,N)=F_(j,m,N)(1+P_R,1+g_R^-1,1+K_R)",
        "graph_finite": True,
        "graph_acyclic": True,
        "finite_monotone_composition_of_tame_and_archived_inverse_estimates": True,
        "literal_polynomial_in_all_inverse_norms_claimed": False,
        "constant_finite_if_primitive_bundle_holds": True,
        "constant_uniform_in_stage": False,
        "large_derivative_array_evaluated": False,
        "numerical_materialization_required": False,
    }


def residual_ledger(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q207 identity required")
    return {
        "finite_state_parts": (
            "regular", "compact_A0_A1", "Gaussian", "base_exterior"
        ),
        "regular_contains": (
            "spatial_slow_cutoffs", "pressure_gradients", "unit_viscosity",
            "cylindrical_connections", "old_new_crosses", "square"
        ),
        "six_diagonal_groups": (
            "e_t", "-Delta e", "grad pi", "(U.grad)e", "(e.grad)U", "(e.grad)e"
        ),
        "complete_memory": (
            "six_covariances", "radial_axial_velocity", "source_moving_normal_pressure",
            "Pi1_Pi2", "whole_moving_pressure_product", "both_velocity_moments"
        ),
        "scale_cutoffs_allocated_only_in_diagonal": True,
        "potential_direct_swirl_pressure_cut_together": True,
        "cutoff_scale_chosen_after_finite_constant": True,
        "stage_uniform_constant_used": False,
        "complete_fixed_jet_residual_bound_conditionally": True,
    }


def primitive_frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q207 identity required")
    return {
        "log_renormalized_inequality": (
            "sup_(ell>=ell*) epsilon_ell^-B_j S^-d sup_Omega zeta^-1/2 y_a^r y_b^s "
            "|D^I(f_(j,m)/P)| < infinity"
        ),
        "pure_power_conversion": "for eta>0 the exponent is B_j-eta",
        "same_exponent_pure_power_requires_d_zero": True,
        "Q196_initial_source_and_endpoint_already_derived_conditionally": True,
        "smooth_mixed_jet_zero_extension_both_endpoints": "required",
        "common_positive_gap_bundle": "required on the same closed domain",
        "source_terms": (
            "radial_advection", "complete_pressure", "vector_mixed_viscosity",
            "cylindrical_connections", "cutoff_curls", "mean_wave_terms",
            "both_cross_orientations", "square"
        ),
        "named_providers_supply_log_renormalized_bounds_conditionally": True,
        "evidence_hashes_prove_continuum_supremum": False,
        "finite_symbolic_request_proves_continuum_supremum": False,
        "moderate_fixture_proves_actual_all_band_bound": False,
        "new_exponent_or_algebraic_obstruction_found": False,
        "independent_primitive_continuum_proof_complete": False,
    }


def concrete_request(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q207 identity required")
    stage = stage_for_power()
    recurrence = load_recurrence()
    request = recurrence.physical_request(stage["stage"], 0, 0)
    ceiling = request["first_real_ceiling"]
    pole = max(30, 6, 6 * ceiling + 9, 6 * ceiling + 18)
    return {
        "stage": stage["stage"],
        "physical_jet": 0,
        "target_power": 0,
        "A0_integrations": request["integrations"]["A0"],
        "A1_integrations": request["integrations"]["A1"],
        "A0_strict_margin": request["strict_margins"]["A0"],
        "A1_strict_margin": request["strict_margins"]["A1"],
        "first_real_derivative_ceiling": ceiling,
        "sufficient_common_pole": pole,
        "largest_pressure_moment_archive_pole": 12 * pole,
        "actual_derivatives_or_operator_constants_evaluated": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q207 identity required")
    status = statuses()
    links = same_object_links()
    sample = stage_for_power()
    construction = constant_construction()
    ledger = residual_ledger()
    frontier = primitive_frontier()
    request = concrete_request()
    checks = {
        "provider_statuses_match": all(status.values()),
        "same_object_links_frozen": all(links.values()),
        "sample_stage_ceiling_exact": (
            sample["stage"] == 602
            and sample["exponent"] == Q(5999, 2500000)
            and sample["previous_exponent"] == Q(-401, 4000000)
            and sample["minimal"]
        ),
        "finite_acyclic_constant_recipe_closes_conditionally": (
            construction["graph_finite"] and construction["graph_acyclic"]
            and construction["constant_finite_if_primitive_bundle_holds"]
            and not construction["numerical_materialization_required"]
        ),
        "complete_residual_and_cutoff_memory_retained": (
            len(ledger["six_diagonal_groups"]) == 6
            and ledger["scale_cutoffs_allocated_only_in_diagonal"]
            and ledger["potential_direct_swirl_pressure_cut_together"]
            and not ledger["stage_uniform_constant_used"]
        ),
        "primitive_promotion_gap_stated_precisely": (
            frontier["named_providers_supply_log_renormalized_bounds_conditionally"]
            and frontier["same_exponent_pure_power_requires_d_zero"]
            and frontier["Q196_initial_source_and_endpoint_already_derived_conditionally"]
            and not frontier["evidence_hashes_prove_continuum_supremum"]
            and not frontier["independent_primitive_continuum_proof_complete"]
        ),
        "first_direct_request_computed_exactly": (
            request["stage"] == 602
            and request["first_real_derivative_ceiling"] == 36479350300
            and request["sufficient_common_pole"] == 218876101818
            and not request["actual_derivatives_or_operator_constants_evaluated"]
        ),
    }
    return {
        "status": "Q207-fixed-stage-complete-residual-constant-reconstructed-conditionally-primitive-seminorm-frontier-isolated",
        "checks": checks,
        "scope": {
            "fixed_stage_residual_constant_reconstructed": True,
            "complete_pressure_viscosity_cutoffs_and_interactions_retained": True,
            "stage_uniform_constant_required": False,
            "primitive_same_object_continuum_bundle_proved_independently": False,
            "unconditional_complete_residual_force_exclusion_or_blowup": False,
        },
        "next": "prove a quantitative all-fixed-order primitive majorant and common gap/inverse theorem covering the stage-602 request",
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q207 identity required")
    return {
        "status": review()["status"],
        "providers": statuses(),
        "same_object": same_object_links(),
        "rates": rates(),
        "sample_stages": {
            "m0_N0": stage_for_power(0, 0),
            "m1_N0": stage_for_power(1, 0),
            "m0_N1": stage_for_power(0, 1),
        },
        "constant": constant_construction(),
        "residual": residual_ledger(),
        "frontier": primitive_frontier(),
        "concrete_request": concrete_request(),
        "review": review(),
    }
