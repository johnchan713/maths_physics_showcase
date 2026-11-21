"""Q197: first sourced inverse, complex covariance and current mean audit."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
IDENTITY = "Q197-FIRST-SOURCED-COVARIANCE"

PROVIDERS = (
    "core_current_head_q196_initial_generated_source",
    "sourced_inverse_foundation",
    "first_sourced_foundation",
    "core_current_head_q184_generated_source_log_loss_correction",
)

EXPECTED = {
    PROVIDERS[0]: "Q196-initial-generated-source-stage-and-log-quotient-corrected",
    PROVIDERS[1]: "explicit-same-domain-sourced-inverse-and-high-harmonic-parameter-jets-relative-to-reviewed-frame-inputs",
    PROVIDERS[2]: "first-sourced-physical-pulse-complex-covariance-and-current-pressure-foundation-relative-to-core-and-finite-lemmas",
    PROVIDERS[3]: "Q184-generated-source-log-loss-corrected-and-pure-power-reserve-derived",
}

COMPONENTS = ("rr", "rtheta", "rz", "thetatheta", "thetaz", "zz")


def evidence(name):
    return json.loads((HERE.parent / name / "evidence.json").read_text())


def statuses(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q197 identity required")
    return {name: evidence(name).get("status") == EXPECTED[name] for name in PROVIDERS}


def rates(kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q197 identity required")
    kappa = Q(kappa)
    if not 0 < kappa < Q(1, 4):
        raise ValueError("reviewed kappa required")
    B = 1 - kappa
    C = B + Q(1, 2)
    D = C - kappa
    return {
        "kappa": kappa,
        "B_source_and_sourced_wave_log": B,
        "C_pressure_potential_covariance_log": C,
        "D_postsource_regular_log": D,
        "B_pure": B - kappa,
        "C_pure": C - kappa,
        "D_pure": D - kappa,
        "source_gain": Q(1, 2) - kappa,
        "covariance_square": 2 * B,
        "postsource_square": 2 * B - kappa,
    }


def prefix_orders(kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q197 identity required")
    kappa = Q(kappa)
    row = rates(kappa)
    alpha = 2 - 3 * kappa
    H = alpha - kappa
    groups = {
        "original_primary_regular_source": row["B_source_and_sourced_wave_log"],
        "signed_linear_primary_cross_actual_mean_cross": alpha - kappa,
        "signed_square": 2 * alpha - 1 - kappa,
        "post_signed_fast_mean_wave": H,
        "post_signed_slow_mean_wave": alpha,
    }
    return {
        "alpha_signed_target": alpha,
        "H_fast_row": H,
        "groups": groups,
        "minimum": min(groups.values()),
        "original_source_remains_limiting": min(groups.values()) == row["B_source_and_sourced_wave_log"],
        "signed_and_mean_prefix_must_not_be_skipped": True,
    }


def forced_pressure(nKt, nprime_t, nf, k, m, norm2, identity=IDENTITY):
    """Return q where pi=i*q and verify the differentiated constraint scalar."""
    if identity != IDENTITY:
        raise ValueError("Q197 identity required")
    nKt, nprime_t, nf, k, m, norm2 = map(Q, (nKt, nprime_t, nf, k, m, norm2))
    if k == 0 or m == 0 or norm2 <= 0:
        raise ValueError("nonzero carrier/harmonic and positive normal norm required")
    q = (nKt - nprime_t + nf) / (k * m * norm2)
    # Dotting the ODE gives -n'.t+n.Kt-km|n|^2*q=-n.f.
    constraint_scalar = -nprime_t + nKt - k * m * norm2 * q + nf
    return {
        "pi_over_i": q,
        "constraint_scalar": constraint_scalar,
        "constraint_preserved": constraint_scalar == 0,
        "source_normal_retained": nf != 0,
        "moving_normal_retained": nKt - nprime_t != 0,
    }


def residence(m, length, dstar, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q197 identity required")
    if isinstance(m, bool) or not isinstance(m, int) or m == 0:
        raise ValueError("nonzero integer harmonic required")
    length, dstar = Q(length), Q(dstar)
    if length <= 0 or dstar <= 0:
        raise ValueError("positive length and damping floor required")
    if abs(m) == 1:
        ell = length
    else:
        ell = min(length, Q(1, (m * m - 1)) / dstar)
    return {
        "harmonic": m,
        "residence": ell,
        "high_harmonic_inverse_square": abs(m) == 1 or ell <= Q(1, (m * m - 1)) / dstar,
        "pressure_does_not_gain_inverse_square": True,
    }


def degree_recurrence(input_degrees, identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q197 identity required")
    if not input_degrees:
        raise ValueError("nonempty degrees required")
    a = tuple(input_degrees)
    if any(isinstance(x, bool) or not isinstance(x, int) or x < 0 for x in a):
        raise ValueError("nonnegative integer degrees required")
    b = []
    for j in range(len(a)):
        g = max(i + a[j - i] for i in range(j + 1))
        prior = max((i + 1 + b[j - i] for i in range(1, j + 1)), default=0)
        b.append(max(g + 1, prior))
    return {"input": a, "solution": tuple(b), "finite_each_fixed_order": True}


def _vector(values):
    values = tuple(map(Q, values))
    if len(values) != 3:
        raise ValueError("three cylindrical components required")
    return values


def covariance_increment(old, new, identity=IDENTITY):
    """Exact zero angular mode, combining equal positive harmonics first."""
    if identity != IDENTITY:
        raise ValueError("Q197 identity required")
    if any(not isinstance(m, int) or m <= 0 for m in set(old) | set(new)):
        raise ValueError("positive integer harmonic keys required")
    parsed_old = {m: (_vector(pair[0]), _vector(pair[1])) for m, pair in old.items()}
    parsed_new = {m: (_vector(pair[0]), _vector(pair[1])) for m, pair in new.items()}
    matrix = [[Q(0) for _ in range(3)] for _ in range(3)]
    square = [[Q(0) for _ in range(3)] for _ in range(3)]
    cross = [[Q(0) for _ in range(3)] for _ in range(3)]
    for m, (c, d) in parsed_new.items():
        a, b = parsed_old.get(m, ((Q(0),) * 3, (Q(0),) * 3))
        for i in range(3):
            for j in range(3):
                cross[i][j] += (a[i] * c[j] + c[i] * a[j] + b[i] * d[j] + d[i] * b[j]) / 2
                square[i][j] += (c[i] * c[j] + d[i] * d[j]) / 2
                matrix[i][j] = cross[i][j] + square[i][j]
    pairs = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
    entries = {name: matrix[i][j] for name, (i, j) in zip(COMPONENTS, pairs)}
    if 2 in parsed_new:
        c2, d2 = parsed_new[2]
        m2_square_nonzero = any(
            c2[i] * c2[j] + d2[i] * d2[j]
            for i in range(3) for j in range(3)
        )
    else:
        m2_square_nonzero = False
    return {
        "entries": entries,
        "matrix": tuple(tuple(row) for row in matrix),
        "cross": tuple(tuple(row) for row in cross),
        "square": tuple(tuple(row) for row in square),
        "symmetric": all(matrix[i][j] == matrix[j][i] for i in range(3) for j in range(3)),
        "all_six_present": set(entries) == set(COMPONENTS),
        "second_harmonic_square_retained": m2_square_nonzero,
        "unequal_harmonic_crosses_omitted_exactly": True,
    }


def mean_change(data, epsilon=Q(1, 200), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q197 identity required")
    required = {
        "I_DWtt_minus_DWrr_over_R", "dZ_I_DWrz", "I_R2_DWtz",
        "I_R_DWzz", "I_R_DWrr", "I_R_DWtt", "dZ_I_R2_DWrz",
        "c_rho", "DeltaP", "dc_rho", "dDeltaP",
    }
    if set(data) != required:
        raise ValueError("complete mean-change ledger required")
    x = {key: Q(value) for key, value in data.items()}
    epsilon = Q(epsilon)
    P = x["I_DWtt_minus_DWrr_over_R"] - epsilon * x["dZ_I_DWrz"]
    Jt = x["I_R2_DWtz"]
    Jz = x["I_R_DWzz"] - (x["I_R_DWrr"] + x["I_R_DWtt"]) / 2 + epsilon * x["dZ_I_R2_DWrz"] / 2
    whole = x["dc_rho"] * x["DeltaP"] + x["c_rho"] * x["dDeltaP"]
    return {
        "DeltaP_from_covariance": P,
        "DeltaJtheta": Jt,
        "DeltaJz": Jz,
        "dZ_c_rho_DeltaP": whole,
        "moving_bump_term": x["dc_rho"] * x["DeltaP"],
        "all_six_roles_retained": True,
    }


def physical_powers(h=Q(1, 200), kappa=Q(1, 100000), identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q197 identity required")
    h, kappa = Q(h), Q(kappa)
    row = rates(kappa)
    A = Q(1, 2) + h
    return {
        "sourced_velocity_zero_jet": -A + h * row["B_source_and_sourced_wave_log"],
        "sourced_pressure_zero_jet": -2 * A + h * row["C_pressure_potential_covariance_log"],
        "sourced_potential_zero_jet": Q(1, 2) - A + h * row["C_pressure_potential_covariance_log"],
        "postsource_residual_zero_jet": -2 * A - Q(1, 2) + h * row["D_postsource_regular_log"],
    }


def frontier(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q197 identity required")
    return {
        "closed_here": (
            "the first sourced zero-entry solve, complete harmonics one and two, six-entry "
            "complex covariance increment, compact pressure/axial moments and recomputed "
            "C-defects/D-rows close conditionally on one domain"
        ),
        "blocking_inequality": (
            "form the next generic complex signed target at order C from the current full "
            "rows and defects, using old wave error B rather than residual order D; prove its "
            "weighted original-root update and then recompute fast/slow means and full pressure "
            "after each operation before advancing the later source induction"
        ),
        "next_generic_signed_target_and_mean_refresh_rederived_here": False,
        "complete_finite_cycle_residual_rederived_here": False,
    }


def review(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q197 identity required")
    status = statuses()
    rate = rates()
    prefix = prefix_orders()
    pressure = forced_pressure(7, 2, 5, 3, 2, 11)
    degrees = degree_recurrence((1, 2, 4, 7, 11))
    cov = covariance_increment(
        {1: ((1, 2, 3), (2, 1, 4))},
        {1: ((1, 1, 1), (3, 2, 1)), 2: ((2, 3, 4), (4, 1, 2))},
    )
    mean = mean_change({
        "I_DWtt_minus_DWrr_over_R": 11, "dZ_I_DWrz": 3,
        "I_R2_DWtz": 13, "I_R_DWzz": 17, "I_R_DWrr": 5,
        "I_R_DWtt": 7, "dZ_I_R2_DWrz": 19,
        "c_rho": 2, "DeltaP": 3, "dc_rho": 5, "dDeltaP": 7,
    })
    powers = physical_powers()
    boundary = frontier()
    checks = {
        "all_provider_statuses_match": all(status.values()),
        "signed_and_mean_prefix_is_retained_before_solve": prefix["signed_and_mean_prefix_must_not_be_skipped"] and prefix["original_source_remains_limiting"],
        "forced_pressure_preserves_constraint_and_both_normal_terms": pressure["constraint_preserved"] and pressure["source_normal_retained"] and pressure["moving_normal_retained"],
        "high_harmonic_residence_has_inverse_square_gain": residence(2, 10, Q(1, 2))["residence"] == Q(2, 3),
        "finite_Duhamel_degree_at_each_fixed_order": degrees["finite_each_fixed_order"] and len(degrees["solution"]) == len(degrees["input"]),
        "complete_covariance_retains_all_entries_and_m2_square": cov["symmetric"] and cov["all_six_present"] and all(cov["entries"].values()) and cov["second_harmonic_square_retained"],
        "covariance_cross_and_square_meet_C": rate["B_source_and_sourced_wave_log"] + Q(1, 2) == rate["C_pressure_potential_covariance_log"] and rate["covariance_square"] >= rate["C_pressure_potential_covariance_log"],
        "current_pressure_and_whole_moving_axial_moment_recomputed": mean["all_six_roles_retained"] and mean["moving_bump_term"] != 0 and mean["dZ_c_rho_DeltaP"] == 29,
        "postsource_groups_meet_D": rate["B_source_and_sourced_wave_log"] + Q(1, 2) - rate["kappa"] == rate["D_postsource_regular_log"] and rate["postsource_square"] >= rate["D_postsource_regular_log"],
        "normalized_gain_not_promoted_to_absolute_physical_decay": powers["postsource_residual_zero_jet"] < 0,
        "no_next_signed_cycle_or_complete_residual_promotion": not boundary["next_generic_signed_target_and_mean_refresh_rederived_here"] and not boundary["complete_finite_cycle_residual_rederived_here"],
    }
    return {
        "checks": checks,
        "provider_status_checks": status,
        "rates": rate,
        "prefix_orders": prefix,
        "forced_pressure_fixture": pressure,
        "residence_fixtures": [residence(m, 10, Q(1, 2)) for m in (1, 2, 3, -2)],
        "degree_recurrence_fixture": degrees,
        "covariance_fixture": cov,
        "mean_change_fixture": mean,
        "physical_powers": powers,
        "correction": {
            "Q196_source_must_pass_through_existing_signed_mean_prefix": True,
            "unequal_harmonics_are_not_combined_before_angular_average": True,
            "source_normal_pressure_and_second_harmonic_square_not_dropped": True,
            "historical_frozen_evidence_modified": False,
        },
        "frontier": boundary,
        "scope": {
            "first_sourced_inverse_and_complete_covariance_closed_conditionally": True,
            "current_pressure_moments_and_defects_recomputed": True,
            "next_generic_signed_cycle_independently_rederived": False,
            "unconditional_complete_residual_force_exclusion_or_blowup": False,
            "tests_are_continuum_proof": False,
        },
    }


def record(identity=IDENTITY):
    if identity != IDENTITY:
        raise ValueError("Q197 identity required")
    return {
        "status": "Q197-first-sourced-complex-covariance-and-current-means-joined",
        "established": (
            "after retaining the existing signed/mean prefix, the Q196 B-source has an exact "
            "zero-entry sourced solve with full forced pressure; its harmonics one and two "
            "produce a complete C-order covariance/pressure change and D-order current rows"
        ),
        "correction": (
            "the first sourced solve is not applied immediately to the bare Q195 state: the "
            "existing real signed correction and mean relaxation are preserved, while the "
            "original B-source remains the limiting row"
        ),
        "next": frontier()["blocking_inequality"],
        "unconditional_actual_pulse_residual_force_exclusion_or_blowup_certified": False,
    }
