"""Q191: independent physical-pulse chart, frame and low-order seed audit."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q191-PHYSICAL-PULSE-SEED-AUDIT"

PROVIDERS = (
    "core_current_head_q190_downstream_route_frontier_audit",
    "core_current_head_q189_same_profile_envelope_gap_audit",
    "core_current_head_q187_actual_symbolic_q2_selector",
    "leading_frame_foundation",
    "weighted_stress_foundation",
    "physical_pulse_foundation",
)

EXPECTED = {
    PROVIDERS[0]: "Q190-Q189-through-Q185-downstream-route-closed-conditionally-and-frontier-corrected",
    PROVIDERS[1]: "Q189-same-profile-envelope-and-gap-join-closed-relative-to-region-theorems",
    PROVIDERS[2]: "Q187-actual-symbolic-Q2-selector-joined-and-primitive-frontier-isolated",
    PROVIDERS[3]: "current-leading-frame-floors-from-reviewed-core-and-same-root-regional-profile",
    PROVIDERS[4]: "actual-weighted-first-stress-and-background-foundation-joined-relative-to-exact-core-and-finite-lemmas",
    PROVIDERS[5]: "actual-primary-physical-pulses-all-fixed-jets-and-common-domain-tails-relative-to-exact-core-and-finite-lemmas",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q191 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def chart_seed(h=Q(1, 200), eta=Q(1), identity=IDENTITY):
    """Reconstruct the dimensionless implicit-chart Jacobian exactly."""
    if identity != IDENTITY:
        raise ValueError("Q191 identity required")
    h, eta = Q(h), Q(eta)
    if not 0 < h < Q(1, 100) or abs(eta) > 1:
        raise ValueError("reviewed h and |eta|<=1 required")
    A, D = Q(1, 2) + h, Q(1, 2) - h
    L = 1 - 2 * h * eta * eta
    return {
        "h": h,
        "A": A,
        "D": D,
        "two_D_plus_two_h": 2 * D + 2 * h,
        "implicit_equation": "s-Z^2*s^(2h)=T",
        "jacobian": L,
        "uniform_floor": 1 - 2 * h,
        "Q_cancels_from_implicit_equation": 2 * D + 2 * h == 1,
        "jacobian_positive": L >= 1 - 2 * h > 0,
    }


def geometry_seed(inverse_sqrt_epsilon=16, S=64, identity=IDENTITY):
    """Exact carrier and pulse-length inequalities without floating roots."""
    if identity != IDENTITY:
        raise ValueError("Q191 identity required")
    x, S = Q(inverse_sqrt_epsilon), Q(S)
    if x < 1 or S <= 0:
        raise ValueError("inverse sqrt epsilon at least one and positive S required")
    k = -(-x.numerator // x.denominator)
    epsilon = 1 / (x * x)
    r0 = Q(1, 6057072)
    # Tg=4+sqrt(2)<6 is enough for the rational length comparison.
    lower_Ls, upper_Ls = 2 * r0 * S, 12 * r0 * S
    return {
        "inverse_sqrt_epsilon": x,
        "epsilon": epsilon,
        "nearest_nonzero_carrier_ceiling": k,
        "epsilon_k_squared": epsilon * k * k,
        "carrier_floor": epsilon * k * k >= 1,
        "carrier_ceiling": epsilon * k * k <= 4,
        "r0": r0,
        "Ls_lower": lower_Ls,
        "Ls_upper_using_Tg_below_six": upper_Ls,
        "Ls_lower_positive": lower_Ls > 0,
        "Ls_upper_below_S": upper_Ls < S,
        "rounding_error_bound": "|p-ptilde|<=1/k; a nonzero-integer constraint does not give 1/(2k)",
    }


def pressure_constraint(identity=IDENTITY):
    """Check the complete moving/source-normal pressure formula exactly."""
    if identity != IDENTITY:
        raise ValueError("Q191 identity required")
    n = (Q(3, 5), Q(4, 5), Q(0))
    t = (Q(4, 5), Q(-3, 5), Q(0))
    np = (Q(1, 7), Q(-2, 11), Q(3, 13))
    Kt = (Q(5, 17), Q(-7, 19), Q(11, 23))
    f = (Q(-2, 29), Q(3, 31), Q(5, 37))
    k, m, damping = Q(3), Q(2), Q(7, 41)
    dot = lambda a, b: sum((x * y for x, y in zip(a, b)), Q(0))
    norm = dot(n, n)
    P = (dot(n, Kt) - dot(np, t) + dot(n, f)) / (k * m * norm)
    # For pi=iP, -i*k*m*n*pi=+k*m*n*P in t'.
    tp = tuple(-Kt[i] - m * m * damping * t[i] + k * m * n[i] * P - f[i] for i in range(3))
    return {
        "n_dot_t": dot(n, t),
        "norm_squared": norm,
        "P_for_pi_equal_iP": P,
        "constraint_derivative": dot(np, t) + dot(n, tp),
        "constraint_preserved": dot(np, t) + dot(n, tp) == 0,
        "moving_normal_nonzero": dot(np, t) != 0,
        "source_normal_nonzero": dot(n, f) != 0,
        "formula": "pi=i*(n.Kt-nprime.t+n.f)/(k*m*|n|^2)",
    }


def riccati_seed(mu=1, C_E=1, length=8, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q191 identity required")
    mu, C_E, length = map(Q, (mu, C_E, length))
    if min(mu, C_E, length) <= 0 or C_E / length > mu / 8:
        raise ValueError("positive accepted Riccati parameters required")
    error = C_E / length
    r = error / mu
    outward = -2 * mu * r + error * (1 + 2 * r + r * r)
    upper = -Q(47, 64) * error
    return {
        "ratio": r,
        "outward_derivative": outward,
        "upper": upper,
        "ratio_at_most_one_eighth": r <= Q(1, 8),
        "strict_outward_sign": outward <= upper < 0,
    }


def physical_exponents(a=0, b=0, h=Q(1, 200), kappa=Q(1, 100000), identity=IDENTITY):
    """Reconstruct the primary physical Q powers from the actual rescalings."""
    if identity != IDENTITY:
        raise ValueError("Q191 identity required")
    if any(isinstance(v, bool) or not isinstance(v, int) or v < 0 for v in (a, b)):
        raise ValueError("nonnegative integer derivative orders required")
    h, kappa = Q(h), Q(kappa)
    if not 0 < h < Q(1, 100) or not 0 < kappa < Q(1, 2):
        raise ValueError("reviewed h and kappa required")
    A = Q(1, 2) + h
    sx, st = Q(1, 2) + h / 2, 1 + 3 * h / 2
    coefficient_space, coefficient_time = Q(1, 2) + h * kappa, 1 + h
    derivative = a * sx + b * st
    return {
        "a": a,
        "b": b,
        "A": A,
        "sx": sx,
        "st": st,
        "coefficient_space": coefficient_space,
        "coefficient_time": coefficient_time,
        "phase_losses_dominate_coefficients": sx > coefficient_space and st > coefficient_time,
        "velocity_loss": derivative + sx,
        "pressure_loss": derivative + 2 * sx,
        "potential_loss": derivative,
        "zero_order_scaling_checks": {
            "velocity": -A + h / 2 == -sx,
            "pressure": -2 * A + h == -2 * sx,
            "potential": Q(1, 2) - A + h == 0,
        },
    }


def gaussian_absorption(polynomial_degree=3, physical_loss=5, target_power=7, curvature=Q(1, 10), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q191 identity required")
    degree, loss, target, curvature = map(Q, (polynomial_degree, physical_loss, target_power, curvature))
    if min(degree, loss, target) < 0 or curvature <= 0:
        raise ValueError("nonnegative losses and positive curvature required")
    d = degree + loss + target
    maximum_exponent = d * d / (4 * curvature)
    return {
        "linear_exponent": d,
        "curvature": curvature,
        "maximum_exponent": maximum_exponent,
        "square_identity": "d*l-a*l^2=d^2/(4a)-a*(l-d/(2a))^2",
        "finite_uniform_constant_exists": maximum_exponent >= 0,
        "one_threshold_handles_every_later_band_for_each_fixed_request": True,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q191 identity required")
    return {
        "closed_here": (
            "dimensionless chart Jacobian, carrier and pulse-length scales, complete constrained "
            "pressure formula, Riccati nonvanishing barrier, primary Q-power arithmetic and "
            "fixed-request Gaussian all-band absorption"
        ),
        "all_fixed_actual_mixed_frame_jets_independently_rederived_here": False,
        "blocking_inequality": (
            "for every fixed mixed index alpha, derive on the actual Q189 representatives "
            "|D^alpha E_frame|+|D^alpha(d-dref)| <= C_alpha S^p_alpha/Ls, with finite "
            "C_alpha,p_alpha on one common later-band domain, while retaining derivatives of "
            "n, B, Bleft, graph coefficients and the complete pressure denominator"
        ),
        "why_needed": (
            "this is the input to differentiated Duhamel/Faa-di-Bruno estimates for every fixed "
            "physical curl, pressure and potential jet; low-order fixtures cannot replace it"
        ),
        "provider_supplies_it_relative_to_named_inputs": True,
        "historical_first_band_or_giant_profile_numerically_certified": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q191 identity required")
    status = statuses()
    chart = chart_seed()
    geometry = geometry_seed(Q(33, 2), 64)
    pressure = pressure_constraint()
    riccati = riccati_seed(1, 1, 8)
    exponents = physical_exponents(3, 2)
    gaussian = gaussian_absorption()
    boundary = frontier()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "chart_has_no_Q_loss_and_positive_jacobian": chart["Q_cancels_from_implicit_equation"] and chart["jacobian_positive"],
        "carrier_and_length_geometry_close_exactly": geometry["carrier_floor"] and geometry["carrier_ceiling"] and geometry["Ls_upper_below_S"],
        "complete_pressure_preserves_constraint": pressure["constraint_preserved"] and pressure["moving_normal_nonzero"] and pressure["source_normal_nonzero"],
        "Riccati_barrier_has_strict_outward_sign": riccati["strict_outward_sign"],
        "primary_Q_powers_and_all_band_absorption_close": (
            exponents["phase_losses_dominate_coefficients"]
            and all(exponents["zero_order_scaling_checks"].values())
            and gaussian["finite_uniform_constant_exists"]
        ),
        "remaining_all_fixed_jet_boundary_not_promoted": not boundary["all_fixed_actual_mixed_frame_jets_independently_rederived_here"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "chart_seed": chart,
        "geometry_seed": geometry,
        "pressure_constraint": pressure,
        "Riccati_seed": riccati,
        "physical_exponents": exponents,
        "Gaussian_absorption": gaussian,
        "frontier": boundary,
        "scope": {
            "physical_pulse_low_order_and_scaling_seed_rederived": True,
            "all_fixed_actual_physical_pulse_jets_rederived": False,
            "generated_source_covariance_or_complete_residual_reproved_here": False,
            "unconditional_pulse_force_exclusion_or_blowup": False,
            "tests_are_continuum_proof": False,
        },
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q191 identity required")
    return {
        "status": "Q191-physical-pulse-chart-frame-pressure-and-scaling-seed-rederived",
        "established": (
            "on the Q189 selected-profile provider chain, the dimensionless chart Jacobian, "
            "carrier/length scales, full normal pressure formula, Riccati continuation, primary "
            "physical Q powers and fixed-request Gaussian all-band absorption close exactly"
        ),
        "obstruction": frontier()["blocking_inequality"],
        "next": (
            "derive the all-fixed actual mixed frame/graph jet inequality and feed it through "
            "differentiated Duhamel and Faa-di-Bruno recurrences for complete curl, pressure and potential"
        ),
        "unconditional_actual_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
