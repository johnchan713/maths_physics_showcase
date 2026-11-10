"""Bind the physical Q152 seed to source, tail, and diagonal contracts."""
from fractions import Fraction as Q
from pathlib import Path
import json
import math

HERE = Path(__file__).resolve().parent
IDENTITY = "Q152-physical-source-tail-complete-potential-diagonal-join"
PROVIDERS = (
    "source_envelope_edge_audit",
    "regular_source_foundation",
    "core_physical_tail_join",
    "core_first_sourced_physical_tail_closure",
    "finite_stage_induction",
    "core_finite_stage_compact_recurrence",
    "core_certified_complete_potential_diagonal_join",
)


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def rates(kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q153 identity required")
    kappa = Q(kappa)
    if not 0 < kappa < Q(1, 6):
        raise ValueError("Inherited kappa range required")
    B = 1 - kappa
    C = B + Q(1, 2)
    D = C - kappa
    delta = Q(1, 2) - 3 * kappa
    return {"B": B, "C": C, "D": D, "delta": delta,
            "physical_seed": D, "next_stage": D + delta}


def bell_quotient(source, log_jets, identity=IDENTITY):
    """Raw one-variable Bell/Leibniz majorant for D^j(f/P)."""
    if identity != IDENTITY or not source or len(source) != len(log_jets):
        raise ValueError("Q153 identity and equal nonempty jet arrays required")
    source = [Q(x) for x in source]
    log_jets = [Q(x) for x in log_jets]
    if any(x < 0 for x in source + log_jets):
        raise ValueError("Nonnegative majorants required")
    R = len(source) - 1
    H = [Q(1)] + [Q(0)] * R
    for j in range(1, R + 1):
        H[j] = sum(Q(math.comb(j - 1, i - 1)) * log_jets[i] * H[j - i]
                   for i in range(1, j + 1))
    G = []
    for j in range(R + 1):
        G.append(sum(Q(math.comb(j, i)) * H[i] * source[j - i]
                     for i in range(j + 1)))
    return {"inverse_gaussian": H, "quotient": G}


def compact_request(m, N, h=Q(1, 200), kappa=Q(1, 100000),
                    identity=IDENTITY):
    """The retained strict physical A0/A1 request at the Q152 seed."""
    if identity != IDENTITY:
        raise ValueError("Q153 identity required")
    if m < 0 or N < 0:
        raise ValueError("Nonnegative physical request required")
    h, kappa = Q(h), Q(kappa)
    if not 0 < h < Q(1, 100) or not 0 < kappa < Q(1, 6):
        raise ValueError("Inherited parameter range required")
    r = rates(kappa)
    Lambda = Q(3, 2) + 2 * h
    target = (Q(N) + Lambda + 2 * m) / h

    def strict(a):
        return max(1, math.floor((target - a) / kappa) + 1)

    p0 = strict(r["B"])
    p1 = strict(r["D"])
    J0 = m + 2 * p0 + 3
    J1 = m + 2 * p1 + 3
    ceiling = max(J1 + 200004, J0 + 200016)
    margin0 = h * (r["B"] + p0 * kappa) - Lambda - 2 * m - N
    margin1 = h * (r["D"] + p1 * kappa) - Lambda - 2 * m - N
    return {"Lambda": Lambda, "target": target, "p0": p0, "p1": p1,
            "J0": J0, "J1": J1, "first_real_ceiling": ceiling,
            "radial_margin": margin0, "axial_margin": margin1}


def cutoff_allocation(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q153 identity required")
    return {
        "per_stage_regular": (
            "spatial/slow cutoff derivatives, pressure gradients, viscosity, "
            "complete curl terms, and all nonlinear interactions"
        ),
        "per_stage_additive": (
            "temporal Gaussian commutators, compact A0/A1 residuals, "
            "and the retained base/exterior flat remainder"
        ),
        "diagonal": (
            "complete-potential, direct-swirl, and pressure scale cutoffs; "
            "their curls, gradients, viscosity, crosses, and square"
        ),
        "double_count_free": True,
    }


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q153 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    required = {
        "source_envelope_edge_audit": (
            "all_fixed_mixed_normalized_source_jets_have_Bell_majorants",
            "generated_regular_source_grammar_has_full_path_and_entrance_inputs_relative_to_current_target_estimates",
            "full_source_normal_moving_normal_and_Pi1_Pi2_pressure_retained",
        ),
        "regular_source_foundation": (
            "full_source_endpoint_extension_and_closed_label_supports_preserved",
            "full_viscosity_pressure_graph_coefficients_and_cylindrical_connections_retained",
            "physically_flat_cutoff_tails_kept_additive_not_forward_inverted",
        ),
        "core_physical_tail_join": (
            "all_fixed_physical_compact_tail_recipes_include_normalization_and_graph_losses",
            "original_temporal_Gaussian_ledger_is_additive_and_never_inverted",
            "same_domain_polynomial_and_band_to_natural_log_Gaussian_constants_derived",
        ),
        "core_first_sourced_physical_tail_closure": (
            "changed_covariance_pressure_and_axial_velocity_retained",
            "first_sourced_compact_graph_acyclic",
            "fixed_m_N_post_sourced_compact_flatness_proved_conditionally",
        ),
        "finite_stage_induction": (
            "every_generated_regular_source_support_Gaussian_envelope_and_endpoint_hypothesis_preserved",
            "actual_flat_velocities_and_exact_excluded_Gaussian_tail_ledger_retained",
            "conditional_complete_four_operation_finite_stage_induction_derived",
        ),
        "core_finite_stage_compact_recurrence": (
            "all_fixed_stage_compact_tails_flat_conditionally",
            "finite_stage_backward_request_recurrence_supplied",
        ),
        "core_certified_complete_potential_diagonal_join": (
            "certified_fixed_stage_domains_joined_to_complete_potential_diagonal",
            "complete_pressure_viscosity_cutoff_crosses_and_square_retained",
            "one_locally_finite_divergence_free_field_constructed_conditionally",
        ),
    }
    for name, keys in required.items():
        flags = data[name].get("flags", {})
        if not all(flags.get(key) for key in keys):
            raise ValueError("Missing provider contract: " + name)
    return {
        "physical_seed": "Q152 same-family D=3/2-2kappa source with pole 9p",
        "source": "all fixed full-path f/P jets and smooth endpoint extension",
        "tails": "all fixed per-stage additive physical tails O(q^N)",
        "stage": "finite acyclic derivative/pole/support ledger for each fixed stage",
        "diagonal": "conditional locally finite complete-potential field with residual flat in every fixed Cartesian jet",
        "independent_continuum_validation": False,
    }


def record():
    return {
        "status": "Q153-physical-source-tail-and-complete-potential-diagonal-bound-conditionally",
        "correction": "Q152's requested source and tail estimates were already proved as family-local conditional theorems; spatial cutoff terms belong to the regular operator and scale cutoffs to the complete-potential diagonal, not to the excluded per-stage tail ledger",
        "established": "Q149-Q152's same physical family satisfies the inherited full-path Gaussian quotient and endpoint source contracts at every fixed requested jet; its temporal Gaussian, compact A0/A1, and base/exterior tails are conditionally flat after every fixed physical derivative, and the corrected D seed binds to the finite-stage recurrence and complete-potential diagonal",
        "scope": "conditional analytical composition of named continuum premises; no giant seminorm is numerically evaluated, no stage-uniform uncut series is claimed, and tests do not prove the continuum estimates",
        "next": "bind this same physical locally finite field, cutoff sequence, support and pressure gauge to the closed-time endpoint atlas and complete one-sided jets before localization, compact forcing, energy, comparison and the core lower bound",
        "physical_generated_source_endpoint_join_certified_conditionally": True,
        "all_fixed_stage_additive_physical_tails_certified_conditionally": True,
        "complete_potential_diagonal_applicable_conditionally": True,
        "closed_time_forcing_energy_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
    }
