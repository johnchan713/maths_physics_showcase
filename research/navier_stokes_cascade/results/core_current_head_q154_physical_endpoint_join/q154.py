"""Bind the Q153 physical diagonal to its closed-time punctured endpoint."""
from fractions import Fraction as Q
from pathlib import Path
import json
import math

HERE = Path(__file__).resolve().parent
IDENTITY = "Q153-physical-diagonal-closed-time-punctured-endpoint-join"
PROVIDERS = (
    "core_current_head_q7_corner_endpoint_join",
    "core_certified_complete_potential_endpoint_join",
    "endpoint_foundation",
    "current_endpoint_jets",
    "physical_pulse_foundation",
    "finite_stage_induction",
)


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def selected_parameters(h=Q(1, 200), kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q154 identity required")
    h, kappa = Q(h), Q(kappa)
    if not 0 < h < Q(1, 100) or not 0 < kappa < Q(1, 6):
        raise ValueError("Inherited parameter range required")
    A = Q(1, 2) + h
    D = Q(1, 2) - h
    return {
        "h": h, "kappa": kappa, "A": A, "D": D,
        "implicit_floor": 1 - 2 * h,
        "stage_seed": Q(3, 2) - 2 * kappa,
        "stage_gain": Q(1, 2) - 3 * kappa,
    }


def band_window(c_exp, qmax_exp, ell0, identity=IDENTITY):
    """Closed dyadic window for c=2^-c_exp, qmax=2^-qmax_exp."""
    if identity != IDENTITY or min(c_exp, qmax_exp, ell0) < 0:
        raise ValueError("Q154 identity and nonnegative exponents required")
    if c_exp < qmax_exp:
        raise ValueError("Need c <= qmax")
    # 1/4 <= q/Q <= 4 gives Q in [q/4,4q].
    lo = max(ell0, qmax_exp - 2)
    hi = c_exp + 2
    return {"lower": lo, "upper": hi, "bands": tuple(range(lo, hi + 1))}


def observation_margin(qstar, ell0, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q154 identity required")
    qstar = Q(qstar)
    if not 0 < qstar < 1 or ell0 < 1:
        raise ValueError("Positive inherited qstar and ell0 required")
    qe = min(qstar / 16, Q(1, 2 ** (ell0 + 2)))
    qmax = 3 * qe / 8
    return {
        "qe": qe,
        "qmax": qmax,
        "sixteen_qmax": 16 * qmax,
        "inside_qstar": 16 * qmax < qstar,
        "first_band_margin": 4 * qmax < Q(1, 2 ** ell0),
    }


def endpoint_requests(spatial, time, identity=IDENTITY):
    if identity != IDENTITY or spatial < 0 or time < 0:
        raise ValueError("Q154 identity and nonnegative orders required")
    m = spatial + time
    return {
        "m": m,
        "potential": m + 4,
        "direct_swirl": m + 3,
        "pressure": m + 2,
        "corner_residual": m + 1,
    }


def rising(x, n):
    x = Q(x)
    out = Q(1)
    for i in range(n):
        out *= x + i
    return out


def exterior_coefficients(j, h=Q(1, 200), identity=IDENTITY):
    if identity != IDENTITY or j < 0:
        raise ValueError("Q154 identity and nonnegative order required")
    h = Q(h)
    if not 0 < h < Q(1, 100):
        raise ValueError("Inherited h range required")
    vals = [rising(h, k) * rising(1 + h, k) for k in range(j + 1)]
    Aj = vals[j]
    Sj = sum(Q(math.comb(j, k)) * vals[k] * vals[j - k]
             for k in range(j + 1))
    return {
        "A_j": Aj,
        "S_j": Sj,
        "swirl_radial_power": 1 + 2 * h + 2 * j,
        "pressure_radial_power": 2 + 4 * h + 2 * j,
        "pressure_denominator": 2 + 4 * h + 2 * j,
    }


def taylor_ratio(j, r=Q(1), h=Q(1, 200), identity=IDENTITY):
    if identity != IDENTITY or j < 0 or Q(r) <= 0:
        raise ValueError("Q154 identity, nonnegative order and positive radius required")
    h, r = Q(h), Q(r)
    return Q(4, 1) * (h + j) * (1 + h + j) / ((j + 1) * r * r)


def provider_contracts(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q154 identity required")
    data = {name: evidence(name) for name in PROVIDERS}
    required = {
        "endpoint_foundation": (
            "same_current_single_realization_is_the_field_under_audit",
            "complete_uniform_next_time_derivative_bounds_joined_from_current_foundations",
            "full_enlarged_closed_band_window_one_quarter_to_four_retained",
            "compatible_potential_direct_swirl_velocity_and_pressure_terminal_jets_supplied_relative_to_retained_lemmas",
            "both_background_pressure_tails_and_fixed_gauge_memory_retained",
        ),
        "current_endpoint_jets": (
            "conditional_complete_potential_direct_swirl_and_pressure_terminal_jets_derived",
            "conditional_uniform_Cauchy_and_finite_Taylor_remainder_bounds_derived",
            "complete_exterior_heat_and_pressure_time_jets_derived",
            "same_frozen_labels_absolute_pullbacks_and_physical_mean_bumps_retained",
        ),
        "physical_pulse_foundation": (
            "one_starting_band_for_all_later_scales_and_every_fixed_physical_derivative_order",
            "exact_closed_supports_and_one_sided_terminal_jets_preserved",
            "uniform_primary_complete_curl_pressure_and_fixed_potential_bounds_derived",
        ),
        "finite_stage_induction": (
            "every_generated_regular_source_support_Gaussian_envelope_and_endpoint_hypothesis_preserved",
            "same_fixed_inverse_domain_and_stage_independent_physical_jet_loss_derived_conditionally",
        ),
    }
    for name, keys in required.items():
        flags = data[name].get("flags", {})
        if not all(flags.get(key) for key in keys):
            raise ValueError("Missing provider contract: " + name)
    q7 = data["core_current_head_q7_corner_endpoint_join"]
    certified = data["core_certified_complete_potential_endpoint_join"]
    if q7["status"] != "current-head-Q7-cartesian-corner-and-punctured-endpoint-joined-conditionally":
        raise ValueError("Q7 endpoint obstruction record required")
    if certified["status"] != "certified-domain-conditional-complete-potential-punctured-endpoint-join":
        raise ValueError("Certified endpoint theorem required")
    return {
        "q7_blocker": "uniform all-fixed-order (A,B,p) bound on each punctured closed-time compact",
        "q7_blocker_supplied_conditionally_by_later_physical_chain": True,
        "same_Q153_field_cutoffs_supports_and_pressure_gauge": True,
        "singular_origin_included": False,
        "future_solution_constructed": False,
    }


def record():
    return {
        "status": "Q154-physical-diagonal-punctured-endpoint-bound-conditionally",
        "correction": "Q7's closed-time inequality was genuinely open at that head, but Q149-Q153 now supply it on the same physical field: one admitted later band tail, one-sided physical pulse jets, finite generated-source/compact graphs and a locally finite complete-potential diagonal",
        "established": "on every punctured terminal compact, every fixed jet of the Q153 complete potential, divergence-safe direct swirl and fixed-gauge pressure has a uniform one-sided bound and compatible terminal limit; full exterior heat and pressure jets retain the complete Leibniz convolution",
        "scope": "conditional analytical join of named continuum providers; the singular origin is excluded, exterior time Taylor series have zero radius when the heat constant is nonzero, and tests are diagnostics rather than continuum proofs",
        "next": "replay the exact axisymmetric space-time localization and zero initialization on this same field, retaining all potential cutoff derivatives, full pressure, unit viscosity, both nonlinear crosses and the cutoff-curl square before smooth force extension",
        "punctured_closed_time_endpoint_certified_conditionally": True,
        "singular_origin_regularized": False,
        "compact_force_energy_comparison_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
    }
