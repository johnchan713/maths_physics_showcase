#!/usr/bin/env python3
"""Exact assembly of the finalized whole-support first U coefficient C_eta^3 bound."""

from fractions import Fraction as F


def rat(value):
    if isinstance(value, (bool, float)):
        raise TypeError("Exact rational inputs required")
    return F(value)


def serial(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: serial(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(item) for item in value]
    return value


def local_value_c3(delta, m_core, m_collar):
    """Tentatively cut F1,U1 value bound; radial cutoff derivatives do not enter."""
    delta, m_core, m_collar = map(rat, (delta, m_core, m_collar))
    if not 0 < delta <= 1 or min(m_core, m_collar) < 0:
        raise ValueError("Positive strip and nonnegative Picard majorants required")
    factor = 1 + 2 / delta + 4 / delta**2 + 8 / delta**3
    return {
        "Ceta3_factor": factor,
        "tentative_value_Ceta3": factor * (m_core + m_collar),
        "cutoff_value_cost": 1,
    }


def reciprocal_cap(h):
    h = rat(h)
    if not 0 < h <= F(1, 100):
        raise ValueError("Retained 0<h<=1/100 required")
    return 1 / (1 - 8 * h)


def signed_source_c3(h, nu, x_velocity, leading_u5):
    """Full-support signed leading-source totals, with the axis viscosity retained."""
    h, nu, x_velocity, leading_u5 = map(
        rat, (h, nu, x_velocity, leading_u5)
    )
    if nu < 0 or x_velocity <= 0 or leading_u5 < 0:
        raise ValueError("Nonnegative viscosity/U bound and positive support required")
    A, D = F(1, 2) + h, F(1, 2) - h
    R = reciprocal_cap(h)
    v4 = R * (4 + 4 * D + 20) * leading_u5
    q03 = (
        2 * nu * v4
        + R
        * x_velocity
        * (8 * D * v4 + (16 + 4 * A) * leading_u5 * v4)
        + x_velocity * v4**2 / 2
    )
    q13 = (
        R
        * x_velocity**2
        / 2
        * ((1 + 8 * D) * v4 + (16 + 4 * D) * leading_u5 * v4)
        + x_velocity**2 * v4**2 / 4
    )
    return {"reciprocal_cap": R, "regular_v_Ceta4": v4, "Q03": q03, "Q13": q13}


def ordinary_discrepancies(x_cut, seed, leading_f3, leading_u3, q03, q13):
    """Five ordinary I3 discrepancies; pressure terms are signed-total bounds."""
    x_cut, seed, leading_f3, leading_u3, q03, q13 = map(
        rat, (x_cut, seed, leading_f3, leading_u3, q03, q13)
    )
    if x_cut <= 0 or min(seed, leading_f3, leading_u3, q03, q13) < 0:
        raise ValueError("Positive cutoff and nonnegative envelopes required")
    return (
        x_cut * seed,
        x_cut**2 * seed,
        2 * leading_f3 * x_cut * seed + q03 / 2,
        x_cut**2 * (leading_u3 + leading_f3) * seed,
        2 * leading_u3 * x_cut * seed
        + leading_f3 * x_cut**2 * seed
        + q13 / 2,
    )


def inverse_800(rho, lam, geometric_k, e_star, discrepancies):
    """The unchanged positive-order linear moment inverse in C_eta^3."""
    rho, lam, geometric_k, e_star = map(rat, (rho, lam, geometric_k, e_star))
    mu = tuple(map(rat, discrepancies))
    if min(rho, lam, geometric_k, e_star) <= 0 or len(mu) != 5 or min(mu) < 0:
        raise ValueError("Positive geometry and five nonnegative discrepancies required")
    m1, m2, m3, m4, m5 = mu
    w3 = 5 / e_star
    rows = (
        m1 / rho,
        (m1 / rho + geometric_k * w3 * m4) / lam,
        geometric_k * m2,
        w3 * m5 / rho,
        w3 * m3,
    )
    return {"W3": w3, "rows": rows, "coefficient_Ceta3": 800 * max(rows)}


def finalized_k13(
    h,
    nu,
    x_velocity,
    leading_u5,
    x_cut,
    local_seed,
    leading_f3,
    leading_u3,
    rho,
    lam,
    geometric_k,
    e_star,
    bump_value=320,
):
    """Assemble Q03,Q13, five targets, inverse-800 and the disjoint U bumps."""
    bump_value = rat(bump_value)
    if bump_value < 0:
        raise ValueError("Nonnegative bump value bound required")
    source = signed_source_c3(h, nu, x_velocity, leading_u5)
    moments = ordinary_discrepancies(
        x_cut,
        local_seed,
        leading_f3,
        leading_u3,
        source["Q03"],
        source["Q13"],
    )
    repair = inverse_800(rho, lam, geometric_k, e_star, moments)
    return {
        "source": source,
        "discrepancies": moments,
        "repair": repair,
        "actual_U_bump_value": bump_value,
        "K13": rat(local_seed) + bump_value * repair["coefficient_Ceta3"],
    }


def selected_formula():
    """Non-materialized formula at the unchanged selected construction parameters."""
    return {
        "delta": "rho_axis/(800 S)",
        "z_core": "2 C_core b_core",
        "z_collar": "2 C_col b_col",
        "Mcore": "(z_core+6 z_core^2/delta) exp(6 z_core^2/delta)",
        "Mcollar": "(z_col+12 z_col^2/delta) exp(12 z_col^2/delta)",
        "J13": "(1+2/delta+4/delta^2+8/delta^3)(Mcore+Mcollar)",
        "leading_inputs": {
            "K05": "H^262",
            "u3": "H^262",
            "f3": "max(Kaxis,2(1+Xa^-1)H^70)",
        },
        "signed_support": "Q03,Q13 integrate over the unchanged full [0,Xv]",
        "repair": "C13=800 max(mu1/rho3,(mu1/rho3+Kg(5/e*)mu4)/lambda,Kg mu2,(5/e*)mu5/rho3,(5/e*)mu3)",
        "final": "K1,3 := J13 + 320 C13",
    }


def review():
    local = local_value_c3(F(1, 4), 7, 11)
    out = finalized_k13(
        F(1, 200), 1, 13, 17, 19, local["tentative_value_Ceta3"], 23, 29,
        31, F(1, 1000), 37, 41,
    )
    no_viscosity = finalized_k13(
        F(1, 200), 0, 13, 17, 19, local["tentative_value_Ceta3"], 23, 29,
        31, F(1, 1000), 37, 41,
    )
    no_pressure = finalized_k13(
        F(1, 200), 1, 13, 0, 19, local["tentative_value_Ceta3"], 23, 29,
        31, F(1, 1000), 37, 41,
    )
    return serial({
        "selected_formula": selected_formula(),
        "exact_fixture": {"local": local, "assembled": out},
        "checks": {
            "local_Ceta3_charges_third_Cauchy_loss": local["Ceta3_factor"] == 585,
            "eta_independent_cutoff_has_unit_value_cost": local["cutoff_value_cost"] == 1,
            "reciprocal_cap_is_25_over_23_at_h_ceiling": reciprocal_cap(F(1, 100)) == F(25, 23),
            "regular_v4_cost_is_below_29": signed_source_c3(F(1, 100), 1, 1, 1)["regular_v_Ceta4"] < 29,
            "axis_viscosity_changes_Q03": out["source"]["Q03"] > no_viscosity["source"]["Q03"],
            "weighted_source_has_no_axis_viscosity_boundary": out["source"]["Q13"] == no_viscosity["source"]["Q13"],
            "pressure_memory_enters_third_and_fifth_targets": (
                out["discrepancies"][2] > no_pressure["discrepancies"][2]
                and out["discrepancies"][4] > no_pressure["discrepancies"][4]
            ),
            "actual_disjoint_U_bump_value_is_320": out["actual_U_bump_value"] == 320,
            "finalized_bound_adds_local_and_repair": out["K13"] == local["tentative_value_Ceta3"] + 320 * out["repair"]["coefficient_Ceta3"],
            "inverse_uses_full_W3": out["repair"]["W3"] == F(5, 41),
        },
    })


if __name__ == "__main__":
    import json
    print(json.dumps(review(), indent=2, sort_keys=True))
