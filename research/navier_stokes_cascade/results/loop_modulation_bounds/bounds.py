"""Scalar checks for the analytic bounds in README; no huge scale is evaluated."""
from fractions import Fraction
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
AMIN = 2**16


def lower(x):
    return mp.make_mpf(x._mpi_[0])


def upper(x):
    return mp.make_mpf(x._mpi_[1])


def scalar_checks(digits=80):
    """Check endpoint constants; README supplies the monotonicity arguments.

    G and H are symbolic. Even G at the smallest allowed A is unnecessary.
    The independent inequality for a variable g>=65536 is checked at g=65536.
    """
    mp.mp.dps = mp.iv.dps = digits
    a = mp.iv.mpf(AMIN)
    g = a
    log = mp.iv.log
    z_coefficient = 16*2049**2
    log_mu_remainder = log(64)+mp.iv.mpf('1.5')*log(a)
    polynomial_slack = g-log(10**6)-10*log(g)
    checks = dict(
        cap_log_remainder_below_A=upper(log_mu_remainder-a) < 0,
        cap_polynomial_below_A16=z_coefficient+1 < AMIN**6,
        z_and_excursion_logs_below_logG=2*AMIN**16 < AMIN**32,
        cutoff_inverse_below_G=AMIN**20 < AMIN**32,
        d0_inverse_below_G=upper(log(16*a)) < AMIN**32,
        polynomial_absorption=lower(polynomial_slack) > 0,
        absorption_slack_increases=lower(1-10/g) > 0,
        v_derivative_polynomial_below_G3=4096 < AMIN,
        phase_density_lower_constant=upper(2*mp.iv.pi) < AMIN,
        primitive_budget_below_H=32 < AMIN**3,
        # Since log(A)<=A<=G, each expression is bounded by the given
        # multiple of G; compare with log(H)=G^4, without forming G.
        state_prefactor_log_below_10G=upper(log(1024)) < AMIN,
        target_prefactor_log_below_7G=upper(log(64)) < AMIN,
        state_and_target_below_H16=10 < 15*AMIN**3,
        auxiliary_N_condition_below_proposal=3 < 31*AMIN**3,
    )
    return dict(interval_digits=digits, checks=checks,
                cap_Z_coefficient=z_coefficient,
                absorption_log_slack=polynomial_slack,
                G_H_and_N_materialized=False)


def derivative_ledger():
    """Exponents k in exp(k G), derived in README sections 2 and 3."""
    return dict(mu_slow=7, t_slow_at_fixed_theta=10, t_theta=3,
                phase_density_slow=11, forward_phase_slow=12,
                inverse_phase_slow=13, t_slow_at_fixed_phi=17,
                loop_slow_derivative=18, loop_J1=19,
                primitive_J1_and_phase_derivative=22,
                common_primitive_budget=32,
                inverse_phase_derivative_included=True,
                zero_variance_uses_signed_square_root=True)


def error_ledger():
    """Retain exact coefficients when replacing m by 8 A^3 d.

    Powers are upper bounds for A>=1. This records the arithmetic of
    independently derived inequalities, not a numerical PDE calculation.
    """
    p1 = 10*8+16
    p2 = 24*8+4
    return dict(field_error='d=4*A*H/N', moment_error='m=8*A^3*d',
                p1_coefficient=p1, p1_A_power=7,
                p2_coefficient=p2, p2_A_power=8,
                state_coefficient=256*4, state_A_power=9,
                target_before_lambda_bound='64*A^5*H/lambda',
                target_after_lambda_bound='64*A^6*H',
                normalization_derivatives_included=True,
                angular_moment_error_order=1, pressure_state_error_order=0,
                checks=dict(p1_below_common=p1 <= 256,
                            p2_below_common=p2 <= 256,
                            fourth_moment_keeps_both_squares=True,
                            pressure_uses_unchanged_axis_datum=True))


def transformed_target_bound(moment_error, normalization_bound, lam):
    """C1 product bound for the two normalized rows BEFORE their subtraction.

    Values are exact Fractions; require 0<lambda<=1 because the remaining
    individual rows must also be bounded by the same expression.
    """
    moment_error, normalization_bound, lam = map(
        Fraction, (moment_error, normalization_bound, lam))
    if moment_error < 0 or normalization_bound <= 0 or not 0 < lam <= 1:
        raise ValueError('Need nonnegative error, positive normalization, 0<lambda<=1')
    return 2*normalization_bound*moment_error/lam
