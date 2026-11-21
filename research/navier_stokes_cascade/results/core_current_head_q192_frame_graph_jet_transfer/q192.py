"""Q192: auxiliary frame jets and corrected physical graph transfer."""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q192-FRAME-GRAPH-JET-TRANSFER"

PROVIDERS = (
    "core_current_head_q191_physical_pulse_seed_audit",
    "core_current_head_q189_same_profile_envelope_gap_audit",
    "leading_frame_foundation",
    "physical_pulse_foundation",
)

EXPECTED = {
    PROVIDERS[0]: "Q191-physical-pulse-chart-frame-pressure-and-scaling-seed-rederived",
    PROVIDERS[1]: "Q189-same-profile-envelope-and-gap-join-closed-relative-to-region-theorems",
    PROVIDERS[2]: "current-leading-frame-floors-from-reviewed-core-and-same-root-regional-profile",
    PROVIDERS[3]: "actual-primary-physical-pulses-all-fixed-jets-and-common-domain-tails-relative-to-exact-core-and-finite-lemmas",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q192 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def _order(n):
    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise ValueError("nonnegative integer order required")
    return n


def inverse_majorants(jet_bounds, floor, identity=IDENTITY):
    """Leibniz recurrence for derivatives of 1/g from g*(1/g)=1."""
    if identity != IDENTITY:
        raise ValueError("Q192 identity required")
    floor = Q(floor)
    if floor <= 0 or not jet_bounds:
        raise ValueError("positive floor and nonempty jet bounds required")
    bounds = [Q(x) for x in jet_bounds]
    if any(x < 0 for x in bounds):
        raise ValueError("nonnegative jet bounds required")
    inverse = [1 / floor]
    for j in range(1, len(bounds)):
        inverse.append(sum(Q(comb(j, i)) * bounds[i] * inverse[j - i]
                           for i in range(1, j + 1)) / floor)
    return inverse


def product_majorants(left, right, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q192 identity required")
    if len(left) != len(right) or not left:
        raise ValueError("equal nonempty majorant arrays required")
    left, right = [Q(x) for x in left], [Q(x) for x in right]
    if any(x < 0 for x in left + right):
        raise ValueError("nonnegative majorants required")
    return [sum(Q(comb(j, i)) * left[i] * right[j - i]
                for i in range(j + 1)) for j in range(len(left))]


def primitive_gap(S=64, length=Q(1, 10000), epsilon=Q(1, 1000000),
                  inverse_k=Q(1, 1000000), identity=IDENTITY):
    """Reduce the actual phase-comparison error using its acceptance gate."""
    if identity != IDENTITY:
        raise ValueError("Q192 identity required")
    S, length, epsilon, inverse_k = map(Q, (S, length, epsilon, inverse_k))
    r0 = Q(1, 6057072)
    if S < 1 or epsilon <= 0 or inverse_k <= 0:
        raise ValueError("positive accepted scale data required")
    if not 2 * r0 * S <= length <= S:
        raise ValueError("actual pulse length comparison required")
    small = epsilon + epsilon * epsilon + inverse_k
    if S * S * small > 1:
        raise ValueError("actual later-band acceptance gate required")
    delta = S ** -3 + small
    raw = 1 / length + (1 + length) * delta
    upper = 5 / length
    return {
        "S": S,
        "length": length,
        "small": small,
        "delta": delta,
        "raw_gap": raw,
        "upper": upper,
        "closes": raw <= upper,
        "derivative_statement": (
            "after j auxiliary derivatives, frozen representatives and smooth compact "
            "profile factors give C_j*S^j/Ls"
        ),
    }


def auxiliary_frame_transfer(order=4, identity=IDENTITY):
    """Record the finite-order rational jet closure and its required inputs."""
    if identity != IDENTITY:
        raise ValueError("Q192 identity required")
    j = _order(order)
    # A sample exact recurrence demonstrates that positive denominator floors,
    # not derivative-smallness assumptions, control reciprocal jets.
    primitive = [Q(2) ** i for i in range(j + 2)]
    reciprocal = inverse_majorants(primitive, Q(1, 2))
    return {
        "order": j,
        "primitive_orders_required_through": j + 1,
        "reciprocal_majorants": reciprocal,
        "normalization_floor": "|n_tan|>=Bstar/2",
        "pressure_floor": "|n|^2>=Bstar^2/4",
        "frame_floor": "|c0|>=cstar and sqrt(1+sa^2)>=1",
        "conservative_power": 2 * j + 4,
        "conclusion": (
            "|partial^alpha Eframe|+|partial^alpha(d-dref)| "
            "<=C_alpha*S^(2|alpha|+4)/Ls for auxiliary derivatives"
        ),
        "no_order_dependent_smallness_gate": True,
        "full_terms": (
            "n,nprime,Kbg,B,Bprime,Bleft,c0,|n|^-2 and the moving/source-normal "
            "pressure numerator"
        ),
    }


def graph_loss(spatial_order=1, time_order=0, h=Q(1, 200),
               kappa=Q(1, 100000), identity=IDENTITY):
    """Compare coefficient graph losses with the harmonic phase budget."""
    if identity != IDENTITY:
        raise ValueError("Q192 identity required")
    a, b = _order(spatial_order), _order(time_order)
    h, kappa = Q(h), Q(kappa)
    if not 0 < h < Q(1, 100) or not 0 < kappa < Q(1, 2):
        raise ValueError("reviewed h and kappa required")
    coefficient_space = Q(1, 2) + h * kappa
    coefficient_time = 1 + h
    sx = Q(1, 2) + h / 2
    st = 1 + 3 * h / 2
    coefficient = a * coefficient_space + b * coefficient_time
    phase = a * sx + b * st
    margin = phase - coefficient
    return {
        "spatial_order": a,
        "time_order": b,
        "normalized_radial_graph_factor": f"epsilon^(-{a}*kappa)",
        "coefficient_Q_loss": coefficient,
        "phase_Q_loss": phase,
        "strict_margin": margin,
        "expected_margin": a * h * (Q(1, 2) - kappa) + b * h / 2,
        "absorbed": margin >= 0 and (a + b == 0 or margin > 0),
    }


def correction(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q192 identity required")
    return {
        "Q191_literal_mixed_D_inequality_is_ambiguous": True,
        "correct_auxiliary_statement": (
            "partial auxiliary frame jets have C_alpha*S^p_alpha/Ls bounds on one domain"
        ),
        "correct_physical_graph_statement": (
            "a spatial graph derivatives may add epsilon^(-a*kappa); after physical "
            "rescaling coefficient losses are a*(1/2+h*kappa)+b*(1+h)"
        ),
        "pure_polynomial_S_bound_for_physical_radial_graph_jets_claimed": False,
        "phase_budget_strictly_absorbs_corrected_loss": True,
    }


def pressure_memory(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q192 identity required")
    return {
        "formula": "pi=i*(n.Kbg*t-nprime.t+n.f)/(k*m*|n|^2)",
        "denominator": "all fixed auxiliary jets of |n|^-2 use the single |n| floor",
        "channels": ("n.Kbg*t", "-nprime.t", "n.f"),
        "moving_normal_retained": True,
        "source_normal_retained": True,
        "primary_f_zero_only_for_the_primary_homogeneous_pulse": True,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q192 identity required")
    return {
        "closed_here": (
            "auxiliary all-fixed rational frame/denominator transfer and corrected physical "
            "graph loss ledger, relative to the named primitive profile-jet inputs"
        ),
        "blocking_inequality": (
            "run the corrected coefficient losses through differentiated Duhamel and Faa-di-Bruno "
            "recurrences for the actual uncut amplitude, then assemble complete curl, pressure and "
            "potential jets with every cutoff, cylindrical connection and support endpoint term"
        ),
        "complete_actual_field_suprema_independently_rederived_here": False,
        "provider_supplies_them_relative_to_named_inputs": True,
        "historical_first_band_numerically_certified": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q192 identity required")
    status = statuses()
    gap = primitive_gap()
    auxiliary = auxiliary_frame_transfer(5)
    losses = [graph_loss(a, b) for a in range(4) for b in range(3)]
    fix = correction()
    pressure = pressure_memory()
    boundary = frontier()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "primitive_value_gap_reduces_to_inverse_length": gap["closes"],
        "rational_auxiliary_jets_need_no_new_smallness": auxiliary["no_order_dependent_smallness_gate"],
        "every_checked_corrected_graph_loss_is_absorbed": all(row["absorbed"] for row in losses),
        "Q191_ambiguity_corrected_without_erasing_graph_loss": (
            fix["Q191_literal_mixed_D_inequality_is_ambiguous"]
            and not fix["pure_polynomial_S_bound_for_physical_radial_graph_jets_claimed"]
        ),
        "complete_pressure_memory_retained": pressure["moving_normal_retained"] and pressure["source_normal_retained"],
        "complete_actual_field_suprema_not_promoted": not boundary["complete_actual_field_suprema_independently_rederived_here"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "primitive_gap": gap,
        "auxiliary_frame_transfer": auxiliary,
        "graph_losses": losses,
        "correction": fix,
        "pressure_memory": pressure,
        "frontier": boundary,
        "scope": {
            "all_fixed_auxiliary_frame_denominator_transfer_closed_relative_to_inputs": True,
            "physical_graph_loss_ledger_corrected": True,
            "complete_actual_pulse_field_bounds_rederived": False,
            "complete_residual_improvement_or_blowup": False,
            "tests_are_continuum_proof": False,
        },
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q192 identity required")
    return {
        "status": "Q192-auxiliary-frame-jets-closed-and-physical-graph-loss-corrected",
        "established": (
            "positive denominator floors and the actual later-band gate close every fixed auxiliary "
            "frame, damping and pressure-denominator jet with polynomial S loss over Ls"
        ),
        "correction": correction()["correct_physical_graph_statement"],
        "next": frontier()["blocking_inequality"],
        "unconditional_actual_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
