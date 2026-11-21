"""Q198: corrected generic complex signed covariance and target audit."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q198-GENERIC-COMPLEX-SIGNED-TARGET"

PROVIDERS = (
    "core_current_head_q197_first_sourced_covariance",
    "core_generic_signed_covariance_join",
    "core_signed_target_edge_budget",
    "core_current_head_q184_generated_source_log_loss_correction",
)

EXPECTED = {
    PROVIDERS[0]: "Q197-first-sourced-complex-covariance-and-current-means-joined",
    PROVIDERS[1]: "conditional-current-generic-signed-covariance-pressure-domain-and-target-join",
    PROVIDERS[2]: "current-signed-target-two-edge-primitive-and-weighted-row-input-budget",
    PROVIDERS[3]: "Q184-generated-source-log-loss-corrected-and-pure-power-reserve-derived",
}


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q198 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def rates(kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q198 identity required")
    kappa = Q(kappa)
    if not 0 < kappa < Q(1, 4):
        raise ValueError("reviewed kappa required")
    B = 1 - kappa
    C = B + Q(1, 2)
    D = C - kappa
    C0 = 2 * B
    return {
        "kappa": kappa,
        "B_current_error_and_signed_principal_log": B,
        "C_target_and_covariance_log": C,
        "D_complete_curl_and_current_source_log": D,
        "C0_covariance_error_and_projected_target_log": C0,
        "target_gain_log": C0 - C,
        "B_pure": B - kappa,
        "C_pure": C - kappa,
        "D_pure": D - kappa,
        "C0_pure": C0 - kappa,
        "target_gain_pure": (C0 - kappa) - (C - kappa),
    }


def covariance_orders(kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q198 identity required")
    row = rates(kappa)
    B = row["B_current_error_and_signed_principal_log"]
    C = row["C_target_and_covariance_log"]
    D = row["D_complete_curl_and_current_source_log"]
    terms = {
        "old_error_signed_principal": 2 * B,
        "complete_wave_signed_curl_remainder": Q(1, 2) + D,
        "complete_signed_square": 2 * B,
    }
    condition = {
        "old_error_margin": B - Q(1, 2),
        "curl_margin": Q(1, 2) - row["kappa"],
        "target_margin": C - 1,
    }
    return {
        "W0_and_U": Q(1, 2),
        "E": B,
        "L": B,
        "R": D,
        "principal_cross": C,
        "error_terms": terms,
        "error_minimum": min(terms.values()),
        "gain_condition_terms": condition,
        "gain_condition_minimum": min(condition.values()),
        "all_errors_tie_at_C0": set(terms.values()) == {row["C0_covariance_error_and_projected_target_log"]},
    }


def _vector(value):
    value = tuple(map(Q, value))
    if len(value) != 3:
        raise ValueError("three components required")
    return value


def _parse(field):
    if any(isinstance(m, bool) or not isinstance(m, int) or m <= 0 for m in field):
        raise ValueError("positive harmonic keys required")
    return {m: (_vector(pair[0]), _vector(pair[1])) for m, pair in field.items()}


def _add(*fields):
    parsed = [_parse(field) for field in fields]
    out = {}
    for m in set().union(*(field.keys() for field in parsed)):
        out[m] = tuple(
            tuple(sum(field.get(m, ((Q(0),) * 3,) * 2)[phase][j] for field in parsed) for j in range(3))
            for phase in range(2)
        )
    return out


def covariance(field_a, field_b=None):
    """Angular zero mode of a tensor b; field_b defaults to field_a."""
    a = _parse(field_a)
    b = a if field_b is None else _parse(field_b)
    matrix = [[Q(0) for _ in range(3)] for _ in range(3)]
    for m in set(a) & set(b):
        ac, asi = a[m]
        bc, bsi = b[m]
        for i in range(3):
            for j in range(3):
                matrix[i][j] += (ac[i] * bc[j] + asi[i] * bsi[j]) / 2
    return tuple(tuple(row) for row in matrix)


def _matrix_add(*matrices):
    return tuple(tuple(sum(matrix[i][j] for matrix in matrices) for j in range(3)) for i in range(3))


def _matrix_sub(a, b):
    return tuple(tuple(a[i][j] - b[i][j] for j in range(3)) for i in range(3))


def generic_covariance_identity(W0, E, L, R, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q198 identity required")
    U = _add(W0, E)
    V = _add(L, R)
    left = _matrix_sub(covariance(_add(U, V)), covariance(U))
    principal = _matrix_add(covariance(W0, L), covariance(L, W0))
    old_error = _matrix_add(covariance(E, L), covariance(L, E))
    curl_error = _matrix_add(covariance(U, R), covariance(R, U))
    square = covariance(V)
    right = _matrix_add(principal, old_error, curl_error, square)
    error = _matrix_add(old_error, curl_error, square)
    return {
        "left": left,
        "right": right,
        "principal": principal,
        "error": error,
        "identity_exact": left == right,
        "all_six_symmetric_entries_present": all(left[i][j] != 0 for i, j in ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))),
        "complex_quadratures_retained": any(pair[1] != (0, 0, 0) for pair in _parse(E).values()),
        "unequal_harmonics_not_crossed": True,
    }


def pressure_projection_orders(kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q198 identity required")
    row = rates(kappa)
    C = row["C_target_and_covariance_log"]
    D = row["D_complete_curl_and_current_source_log"]
    C0 = row["C0_covariance_error_and_projected_target_log"]
    return {
        "complete_covariance_change": C,
        "Haar_pressure_change": C,
        "pointwise_pressure_change": D,
        "scalar_and_radial_Haar_change": C,
        "axial_covariance_pressure_response": C + 1,
        "projected_tangential_target": C0,
        "projection_after_differentiation": True,
        "moving_projection_commutator_retained": True,
        "whole_c_rho_DeltaP_derivative_retained": True,
        "pressure_or_scalar_promoted_to_C0": False,
    }


def moving_derivatives(epsilon, dsigma, moment, sigma, dmoment, dc_rho, delta_p, c_rho, ddelta_p, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q198 identity required")
    values = tuple(map(Q, (epsilon, dsigma, moment, sigma, dmoment, dc_rho, delta_p, c_rho, ddelta_p)))
    epsilon, dsigma, moment, sigma, dmoment, dc_rho, delta_p, c_rho, ddelta_p = values
    return {
        "projected_derivative": epsilon * (dmoment - dsigma * moment - sigma * dmoment),
        "projection_commutator": epsilon * dsigma * moment,
        "whole_pressure_product_derivative": dc_rho * delta_p + c_rho * ddelta_p,
        "moving_terms_nonzero": epsilon * dsigma * moment != 0 and dc_rho * delta_p != 0,
    }


def pressure_domain(n, p, edge_allowance, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q198 identity required")
    if any(isinstance(x, bool) or not isinstance(x, int) or x < 0 for x in (n, edge_allowance)):
        raise ValueError("nonnegative integer output and edge allowance required")
    if isinstance(p, bool) or not isinstance(p, int) or p <= 0:
        raise ValueError("positive integer old pole required")
    q = max(3 * p, 6 * (n + 2) + 9, 6 * (n + 2) + edge_allowance + 6)
    return {
        "common_wave_pole": q,
        "covariance_pressure_source_pole": 2 * q,
        "old_largest_declared_pole": 12 * p,
        "old_domain_suffices": 12 * p >= 2 * q,
        "covariance_input_order": n + 2,
        "target_input_order": n + 3,
    }


def edge_counterexample(N, identity=IDENTITY):
    """Unweighted boundedness does not imply divisibility by zeta."""
    if identity != IDENTITY:
        raise ValueError("Q198 identity required")
    if isinstance(N, bool) or not isinstance(N, int) or N < 2:
        raise ValueError("integer N at least two required")
    x = Q(1, N)
    zeta = x * x * (1 - x) * (1 - x)
    return {
        "unweighted_value": Q(1),
        "zeta": zeta,
        "weighted_quotient": 1 / zeta,
        "grows_with_N": True,
    }


def physical_residual_power(h=Q(1, 200), kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q198 identity required")
    h = Q(h)
    if not 0 < h < Q(1, 100):
        raise ValueError("reviewed h required")
    row = rates(kappa)
    A = Q(1, 2) + h
    return -2 * A - Q(1, 2) + h * row["D_complete_curl_and_current_source_log"]


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q198 identity required")
    return {
        "closed_here": (
            "the corrected Q197 rates close the generic complex covariance algebra and "
            "projected target exponent C0, retaining complete pressure and enlarged finite domains"
        ),
        "blocking_inequality": (
            "for every requested mixed jet I, prove the actual same-family current target bound "
            "|D^I Xi| <= K_I epsilon^C zeta y_a^(-p_I)y_b^(-q_I), together with "
            "the separate original-principal, current-error, root, chart, support and compact "
            "operator bounds through the expanded 2q pressure pole"
        ),
        "actual_weighted_input_inequality_proved_here": False,
        "post_signed_fast_slow_refresh_proved_here": False,
        "complete_residual_or_blowup_proved_here": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q198 identity required")
    status = statuses()
    rate = rates()
    orders = covariance_orders()
    fixture = generic_covariance_identity(
        {1: ((1, 2, 3), (2, 1, 4))},
        {1: ((2, 1, 1), (1, 3, 2)), 2: ((1, 1, 2), (3, 1, 1))},
        {1: ((1, 3, 2), (2, 2, 1))},
        {1: ((1, 1, 3), (1, 2, 2)), 2: ((2, 1, 1), (1, 2, 3))},
    )
    projection = pressure_projection_orders()
    moving = moving_derivatives(2, 3, 5, 7, 11, 13, 17, 19, 23)
    domain = pressure_domain(0, 1, 12)
    edge_small = edge_counterexample(10)
    edge_large = edge_counterexample(100)
    boundary = frontier()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "corrected_Q197_B_not_residual_D_is_old_wave_error": orders["E"] == rate["B_current_error_and_signed_principal_log"],
        "generic_covariance_identity_exact_all_six_complex_entries": fixture["identity_exact"] and fixture["all_six_symmetric_entries_present"] and fixture["complex_quadratures_retained"],
        "all_three_covariance_errors_tie_at_C0": orders["all_errors_tie_at_C0"] and orders["error_minimum"] == rate["C0_covariance_error_and_projected_target_log"],
        "strict_target_gain_is_half_minus_kappa": rate["target_gain_log"] == Q(1, 2) - rate["kappa"] == orders["gain_condition_minimum"],
        "pressure_scalar_and_pointwise_rows_not_promoted": projection["Haar_pressure_change"] == rate["C_target_and_covariance_log"] and projection["pointwise_pressure_change"] == rate["D_complete_curl_and_current_source_log"] and not projection["pressure_or_scalar_promoted_to_C0"],
        "moving_projection_and_pressure_product_terms_retained": moving["moving_terms_nonzero"] and moving["projection_commutator"] != 0 and moving["whole_pressure_product_derivative"] != 0,
        "new_pressure_domain_strictly_exceeds_old_example": domain["common_wave_pole"] == 30 and domain["covariance_pressure_source_pole"] == 60 and not domain["old_domain_suffices"],
        "unweighted_bound_does_not_imply_zeta_divisibility": edge_large["weighted_quotient"] > edge_small["weighted_quotient"] > 1,
        "normalized_current_source_is_not_absolute_decay": physical_residual_power() < 0,
        "actual_weighted_inputs_means_and_complete_residual_not_promoted": not boundary["actual_weighted_input_inequality_proved_here"] and not boundary["post_signed_fast_slow_refresh_proved_here"] and not boundary["complete_residual_or_blowup_proved_here"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "rates": rate,
        "covariance_orders": orders,
        "covariance_fixture": fixture,
        "pressure_projection_orders": projection,
        "moving_derivative_fixture": moving,
        "pressure_domain_fixture": domain,
        "edge_counterexamples": [edge_small, edge_large],
        "physical_residual_power": physical_residual_power(),
        "frontier": boundary,
        "scope": {
            "generic_signed_covariance_and_projected_target_exponents_closed_conditionally": True,
            "complete_pressure_memory_and_finite_domain_expansion_retained": True,
            "actual_signed_pulse_or_post_signed_mean_refresh_certified": False,
            "unconditional_complete_residual_force_exclusion_or_blowup": False,
            "tests_are_continuum_proof": False,
        },
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q198 identity required")
    return {
        "status": "Q198-corrected-generic-complex-signed-target-and-weighted-input-obstruction",
        "established": (
            "at the corrected Q197 rates, the complete generic complex covariance error and "
            "projected current target have order C0=2-2kappa with gain 1/2-kappa over C"
        ),
        "correction": (
            "the C0 cancellation does not promote the Haar pressure, scalar defects or pointwise "
            "pressure, and the next pressure proof requires the enlarged 2q source domain"
        ),
        "next": frontier()["blocking_inequality"],
        "unconditional_actual_signed_pulse_means_residual_force_exclusion_or_blowup_certified": False,
    }
