"""Bind finite coefficient supports and reserved pressure to the Q147 profile."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import sys

HERE = Path(__file__).resolve().parent
IDENTITY = "Q147-final-profile-finite-coefficient-support-and-pressure-memory"
GLOBAL = (
    "M_infty=0", "R_infty=0", "J_infty=0", "S_infty=0",
    "Pi_axis+Cp_infty=0",
)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


previous = load("q148_previous", "core_current_head_q147_repair_pressure_join/q147.py")
support = load("q148_support", "core_coefficient_support_recheck/proof_checks.py")


def finalized_profile(closure, identity):
    """Check that the support theorem receives Q147's actual final profile."""
    if identity != IDENTITY:
        raise ValueError("Q148 same-profile identity required")
    if closure.get("stage") != "same-profile-compact-repair-cone-and-global-pressure-join":
        raise ValueError("Q147 compact/global closure required")
    if tuple(closure.get("global_identities", ())) != GLOBAL:
        raise ValueError("All five finalized functional identities required")
    if not closure.get("same_exact_repair_and_inner_outer_value_cone_certified"):
        raise ValueError("The accepted Q147 repair and cone join are required")
    if not closure.get("global_pressure_and_five_moment_identities_certified"):
        raise ValueError("The unchanged-gauge global pressure identity is required")
    if "unchanged gauge" not in closure.get("pressure", ""):
        raise ValueError("Pressure memory may not be reset at an interface")
    return {
        "stage": "Q147-finalized-leading-profile-for-positive-orders",
        "global_identities": GLOBAL,
        "pressure_gauge": closure["pressure"],
        "leading_I_row": "renormalized R_infty=0; positive orders use ordinary compact moments",
        "same_profile": True,
    }


def finite_order(profile, n, norms, k, rho, K, e_star, lam, identity):
    """Return the fixed-order coefficient estimate and exact support ledger."""
    if identity != IDENTITY or profile.get("stage") != "Q147-finalized-leading-profile-for-positive-orders":
        raise ValueError("Finalized Q147 profile required")
    if type(n) is not int or n < 1:
        raise ValueError("Positive finite order required")
    bound = support.coefficient_bound(norms, k, rho, K, e_star, lam)
    return {
        "stage": "fixed-positive-order-linear-five-row-repair",
        "order": n,
        "coefficient_norm_bound": bound,
        "inverse_bound": 800,
        "rows": ("M", "(M-J)/lambda", "I", "-S", "Cp"),
        "support": support.support_edges(n),
        "velocity_on_I4": "leading only, including potential-cutoff transitions",
        "pressure_on_I4": (
            "Pi1=K1(eta)-mu^2/(4X)", "Pi2=K2(eta)", "Pi_n=0 for n>=3"
        ),
        "finite_order_only": True,
    }


def pressure_obstruction(Xa, Xb, mu, chi1, q2h, identity):
    """Sharp two-point lower bound after arbitrary X-independent cancellation."""
    if identity != IDENTITY:
        raise ValueError("Q148 identity required")
    floor = support.pressure_floor(Q(Xa), Q(Xb), Q(mu), Q(chi1), Q(q2h))
    return {
        "necessary_supremum_floor": floor,
        "identity": "P(Xb)-P(Xa)=chi1*q^(2h)*mu^2/4*(Xa^-1-Xb^-1)",
        "independent_of": ("K1", "K2"),
        "sharp_after_arbitrary_X_independent_cancellation": True,
    }


def uniform_schedule_obligation(n, m, Kcut, Mnm, an_to_minus_nh, h, identity):
    """Expose, without claiming it, the retained all-scale Borel inequality."""
    if identity != IDENTITY or type(n) is not int or type(m) is not int or not (1 <= n and 0 <= m <= n):
        raise ValueError("Fixed n>=1 and 0<=m<=n required")
    Kcut, Mnm, an_to_minus_nh, h = map(Q, (Kcut, Mnm, an_to_minus_nh, h))
    if min(Kcut, Mnm, an_to_minus_nh, h) <= 0:
        raise ValueError("Positive schedule inputs required")
    Cnm = max(Q(1), 2 * Kcut * (1 + 2 * n * h) ** m * Mnm)
    # The retained schedule supplies this weight. Accepting it explicitly
    # keeps the diagnostic exact instead of approximating a rational power.
    lhs = Cnm * an_to_minus_nh
    return {
        "C_nm": Cnm,
        "a_n_to_minus_nh": an_to_minus_nh,
        "left_hand_side": lhs,
        "right_hand_side": Q(1, 2 ** n),
        "holds": lhs <= Q(1, 2 ** n),
        "required_uniformly": "for every retained n and every m<=n on the full physical hull",
        "not_proved_by_fixed_order_bound": True,
    }


def record():
    return {
        "status": "Q147-final-profile-finite-coefficient-supports-and-reserved-pressure-obstruction-joined",
        "established": "every fixed positive order has the ordinary five-row linear repair, inverse bound below 800, finite mixed seminorms and exact velocity/pressure/stress support edges",
        "reserved_mean": "positive-order velocity vanishes on I4 but Pi1=K1-mu^2/(4X) and Pi2=K2 retain complete later-source memory",
        "pressure_floor": "sup_I4 |q^(2A)p_B-Pi0| >= |chi(c1q)|q^(2h)mu^2(e^-6-e^-11)/(8rho3)",
        "uniform_blocker": "prove C_nm*a_n^(-nh)<=2^-n for C_nm=max(1,2K_(m+1)(1+2nh)^m M_(n,m)), every retained n and m<=n, with actual full-hull tuple norms",
        "next": "join the already-audited physical primary-pulse transfer only after the actual schedule, inverse, gap and support premises are bound on this same finalized profile",
        "uniform_actual_pulses_complete_covariance_residual_or_blowup_certified": False,
        "historical_frozen_records_modified": False,
    }
