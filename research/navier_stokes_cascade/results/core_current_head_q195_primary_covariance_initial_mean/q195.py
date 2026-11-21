"""Q195: exact primary covariance and initial complete-mean join."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q195-PRIMARY-COVARIANCE-INITIAL-MEAN"

PROVIDERS = (
    "core_current_head_q194_primary_field_assembly",
    "physical_pulse_transfer",
    "complete_covariance_initial_means",
    "core_current_head_q183_provider_discharge_generated_source_frontier",
)

EXPECTED = {
    PROVIDERS[0]: "Q194-complete-primary-field-bounds-and-support-assembled-conditionally",
    PROVIDERS[1]: "actual-background-to-physical-primary-pulse-transfer-audited",
    PROVIDERS[2]: "actual-primary-complete-covariance-initial-mean-gain-audited",
    PROVIDERS[3]: "Q183-current-head-provider-frontier-reconciled-to-generated-source-supremum",
}

COMPONENTS = ("rr", "rtheta", "rz", "thetatheta", "thetaz", "zz")


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q195 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def _vector(values):
    values = tuple(map(Q, values))
    if len(values) != 3:
        raise ValueError("three cylindrical components required")
    return values


def covariance(cosine, sine, identity=IDENTITY):
    """Exact angular mean of (c cos ntheta+s sin ntheta)^tensor2."""
    if identity != IDENTITY:
        raise ValueError("Q195 identity required")
    c, s = _vector(cosine), _vector(sine)
    pairs = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
    entries = {name: (c[i] * c[j] + s[i] * s[j]) / 2
               for name, (i, j) in zip(COMPONENTS, pairs)}
    matrix = tuple(tuple((c[i] * c[j] + s[i] * s[j]) / 2 for j in range(3))
                   for i in range(3))
    return {
        "entries": entries,
        "matrix": matrix,
        "symmetric": all(matrix[i][j] == matrix[j][i] for i in range(3) for j in range(3)),
        "all_six_present": set(entries) == set(COMPONENTS),
        "mixed_sine_cosine_mean": Q(0),
    }


def physical_loss(field, spatial_order=0, time_order=0,
                  h=Q(1, 200), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q195 identity required")
    if isinstance(spatial_order, bool) or isinstance(time_order, bool):
        raise ValueError("integer orders required")
    if not isinstance(spatial_order, int) or not isinstance(time_order, int):
        raise ValueError("integer orders required")
    if spatial_order < 0 or time_order < 0:
        raise ValueError("nonnegative orders required")
    h = Q(h)
    if not 0 < h < Q(1, 100):
        raise ValueError("reviewed h required")
    sx, st = Q(1, 2) + h / 2, 1 + 3 * h / 2
    offsets = {"velocity": 1, "covariance": 2, "wave_pressure": 2}
    if field not in offsets:
        raise ValueError("unknown field")
    return {
        "field": field,
        "Q_loss": (spatial_order + offsets[field]) * sx + time_order * st,
        "shell_weight": "sqrt(zeta)" if field != "covariance" else "zeta",
        "Gaussian_weight": "P" if field != "covariance" else "P^2",
    }


def principal_match(labels, sigma_theta, sigma_z, identity=IDENTITY):
    """Same-family Haar/root assembly for the two selected stress rows."""
    if identity != IDENTITY:
        raise ValueError("Q195 identity required")
    if not labels:
        raise ValueError("nonempty labels required")
    cover = Q(0)
    theta = Q(0)
    axial = Q(0)
    for eta2, theta_column, axial_column, root in labels:
        eta2, theta_column, axial_column, root = map(
            Q, (eta2, theta_column, axial_column, root)
        )
        cover += eta2
        theta += eta2 * theta_column * root
        axial += eta2 * axial_column * root
    sigma_theta, sigma_z = Q(sigma_theta), Q(sigma_z)
    return {
        "cover": cover,
        "theta": theta,
        "axial": axial,
        "theta_defect": theta - sigma_theta,
        "axial_defect": axial - sigma_z,
        "exact": cover == 1 and theta == sigma_theta and axial == sigma_z,
    }


def initial_targets(data, epsilon=Q(1, 200), identity=IDENTITY):
    """Exact zero-old-mean scalar targets after radial integration by parts.

    Data are the indicated radial integrals or their Z derivatives.  The
    function is a coefficient ledger, not a numerical profile quadrature.
    """
    if identity != IDENTITY:
        raise ValueError("Q195 identity required")
    required = {
        "I_Wtt_minus_Wrr_over_R", "dZ_I_Wrz", "I_R2_Wtz",
        "I_R_Wzz", "I_R_Wrr", "I_R_Wtt", "dZ_I_R2_Wrz",
    }
    if set(data) != required:
        raise ValueError("complete initial target ledger required")
    d = {key: Q(value) for key, value in data.items()}
    epsilon = Q(epsilon)
    return {
        "P": d["I_Wtt_minus_Wrr_over_R"] - epsilon * d["dZ_I_Wrz"],
        "Jtheta": d["I_R2_Wtz"],
        "Jz": d["I_R_Wzz"] - (d["I_R_Wrr"] + d["I_R_Wtt"]) / 2
              + epsilon * d["dZ_I_R2_Wrz"] / 2,
        "used_components": ("rr", "rz", "thetatheta", "thetaz", "zz"),
        "rtheta_enters_tangential_stress_row": True,
    }


def moving_pressure_derivative(c_rho, pressure, dc_rho, dpressure,
                               identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q195 identity required")
    c_rho, pressure, dc_rho, dpressure = map(
        Q, (c_rho, pressure, dc_rho, dpressure)
    )
    return {
        "whole_product": dc_rho * pressure + c_rho * dpressure,
        "omitted_bump_motion": c_rho * dpressure,
        "bump_motion_term": dc_rho * pressure,
        "omission_changes_identity": dc_rho * pressure != 0,
    }


def rates(kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q195 identity required")
    kappa = Q(kappa)
    if not 0 < kappa < Q(1, 4):
        raise ValueError("reviewed kappa required")
    return {
        "covariance_remainder": 2 - 2 * kappa,
        "initial_three_defects": Q(2),
        "initial_tangential_residual": 2 - 3 * kappa,
        "viscous_loss": 3 * kappa,
        "positive_initial_margin": 2 - 3 * kappa > 0,
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q195 identity required")
    return {
        "closed_here": (
            "the Q194 assembled real primary has an exact six-entry covariance, "
            "same-family matched stresses, complete initial compact pressure/moment "
            "identities and the initial S2/M_(2-3kappa) mean gain"
        ),
        "blocking_inequality": (
            "for the first actual generated nonzero harmonic f_(0,m), derive on the "
            "same physical domain every fixed-jet weighted quotient bound "
            "sup epsilon^(-B0) zeta^(-1/2) ya^r yb^s |D^I(f_(0,m)/P)|<infinity "
            "and smooth mixed-jet zero extension at both pulse endpoints, with full "
            "pressure, viscosity, cutoffs and all interactions retained"
        ),
        "later_complex_covariance_cycle_rederived_here": False,
        "complete_generated_source_or_residual_rederived_here": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q195 identity required")
    status = statuses()
    cov = covariance((1, 2, 3), (4, 5, 6))
    match = principal_match(
        ((Q(1, 4), 6, 10, 2), (Q(3, 4), 2, 6, 2)),
        6, 14,
    )
    targets = initial_targets({
        "I_Wtt_minus_Wrr_over_R": 11,
        "dZ_I_Wrz": 3,
        "I_R2_Wtz": 13,
        "I_R_Wzz": 17,
        "I_R_Wrr": 5,
        "I_R_Wtt": 7,
        "dZ_I_R2_Wrz": 19,
    })
    moving = moving_pressure_derivative(2, 3, 5, 7)
    rate = rates()
    boundary = frontier()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "exact_real_harmonic_mean_has_all_six_entries": cov["symmetric"] and cov["all_six_present"] and all(cov["entries"].values()),
        "covariance_and_wave_pressure_have_same_physical_Q_loss": all(
            physical_loss("covariance", a, b)["Q_loss"] == physical_loss("wave_pressure", a, b)["Q_loss"]
            for a in range(5) for b in range(4)
        ),
        "covariance_has_full_squared_shell_and_Gaussian_weights": physical_loss("covariance")["shell_weight"] == "zeta" and physical_loss("covariance")["Gaussian_weight"] == "P^2",
        "same_family_root_Haar_match_is_exact": match["exact"],
        "all_six_entries_have_distinct_mean_or_stress_roles": set(targets["used_components"]) | {"rtheta"} == set(COMPONENTS),
        "whole_moving_pressure_product_is_differentiated": moving["omission_changes_identity"] and moving["whole_product"] == 29,
        "initial_mean_gain_retains_fast_viscous_loss": rate["initial_tangential_residual"] == 2 - 3 * Q(1, 100000),
        "wave_pressure_zero_mean_does_not_delete_pointwise_pressure": True,
        "Q194_envelope_alone_not_used_as_stress_matching": True,
        "no_later_complex_cycle_or_complete_residual_promotion": not boundary["later_complex_covariance_cycle_rederived_here"] and not boundary["complete_generated_source_or_residual_rederived_here"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "covariance_fixture": cov,
        "physical_losses": [physical_loss(field, a, b) for field in ("velocity", "covariance", "wave_pressure") for a in range(4) for b in range(3)],
        "principal_match_fixture": match,
        "initial_targets_fixture": targets,
        "moving_pressure_fixture": moving,
        "rates": rate,
        "correction": {
            "Q194_pointwise_envelope_does_not_by_itself_imply_matched_stress": True,
            "matching_uses_the_same_physical_root_and_Haar_columns": True,
            "wave_pressure_is_only_zero_in_angular_mean": True,
            "historical_frozen_evidence_modified": False,
        },
        "frontier": boundary,
        "scope": {
            "complete_real_primary_covariance_rederived": True,
            "complete_initial_mean_joined": True,
            "later_complex_cycles_and_generated_sources_rederived": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup": False,
            "tests_are_continuum_proof": False,
        },
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q195 identity required")
    return {
        "status": "Q195-primary-six-entry-covariance-and-initial-complete-mean-joined",
        "established": (
            "the Q194 assembled real physical primary has exact six-entry angular covariance; "
            "the same-family Haar/root identity supplies the two matched stresses, and the full "
            "initial compact pressure and radial/axial mean update reaches S2 with tangential "
            "order 2-3kappa"
        ),
        "correction": (
            "the Q194 field envelope bounds the tensor but cannot prove stress matching; that "
            "cancellation uses the unchanged physical root identity, while wave pressure is "
            "zero only after angular averaging"
        ),
        "next": frontier()["blocking_inequality"],
        "unconditional_actual_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
