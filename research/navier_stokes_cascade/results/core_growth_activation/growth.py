"""Core source/growth constants and the frozen activation budget.

Only logarithms of C, R and activation widths are evaluated. Scalar and
atomic-pressure samples are diagnostics, not the giant constructed profile.
The uniform arguments and normalized moment factor counts are in README.md.
"""
from fractions import Fraction as Q
import importlib.util
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('frozen_core_contraction', HERE.parent/'core_contraction_audit/core.py')
core = importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)


def f0(z, derivative=0):
    z = mp.mpf(z)
    if not 0 <= z <= mp.mpf('4.1') or type(derivative) is not int or derivative < 0:
        raise ValueError('Scalar comparison domain or derivative order invalid')
    return (-mp.mpf('.5'))**derivative/mp.factorial(derivative+1) \
        *mp.hyper([], [derivative+2], -z/2)


def comparison(z):
    z = mp.mpf(z)
    f, fp = f0(z), f0(z, 1)
    return dict(f=f, fp=fp, combined=f+z*fp, shear=-2*z*fp/f)


def pressure_axis(eta, P=16, h='.005', j='1e-18', sigma='1e-20', atoms=None):
    """Evaluate an explicit diagnostic measure with the required ideal atom.

    This fixture is not substituted for the defining schedule measure.
    Derivatives and the tiny eta/j scale are retained at arbitrary precision.
    """
    eta, P, h, j, sigma = map(mp.mpf, (eta, P, h, j, sigma))
    atoms = [(mp.mpf(1), mp.mpf('2.5')), (mp.mpf('.3'), mp.mpf('1.2'))] if atoms is None else [tuple(map(mp.mpf, a)) for a in atoms]
    if abs(eta) > 1 or P < 16 or not 0 < h <= mp.mpf('.01') or not 0 < j <= mp.mpf('.05') or sigma <= 0:
        raise ValueError('Invalid axis data')
    if any(not 0 <= theta <= 1 or mass < 0 for theta, mass in atoms) or sum(m for _, m in atoms) >= 5:
        raise ValueError('Invalid positive pressure measure')
    if sum(mass for theta, mass in atoms if theta == 1) < mp.mpf('2.5'):
        raise ValueError('The exit-growth theorem requires the ideal atom')
    A, D, d, L = mp.mpf('.5')+h, mp.mpf('.5')-h, 1-eta*eta, 1-2*h*eta*eta
    U = 4*eta+j
    H = D*eta+d*U
    Hp = D-8*eta*eta+d*4-2*eta*j
    W = 1-4*d-2*D*eta*U
    base = [(theta, mass*(1+eta*eta)**(-2*theta)) for theta, mass in atoms]
    Pi = -P*P*mp.fsum(mass for _, mass in base)
    Pie = P*P*mp.fsum(4*theta*eta/(1+eta*eta)*mass for theta, mass in base)
    Kpressure = mp.fsum((4*d*theta/(1+eta*eta)+4*A)*mass for theta, mass in base)
    velocity = -A*(1-2*eta*U)*U-4*H
    Z = velocity-d*Pie+4*A*eta*Pi
    chi = H*H/(H*H+sigma*sigma)
    chie = 2*H*Hp*sigma*sigma/(H*H+sigma*sigma)**2
    return dict(U=U, H=H, Hp=Hp, W=W, Pi=Pi, Pie=Pie, Kpressure=Kpressure,
                velocity=velocity, Z=Z, chi=chi, chie=chie,
                zeta=-L*H/(H*H+sigma*sigma), L=L)


def algebra_constants():
    return dict(reference_pressure_mass=(5+mp.exp(mp.mpf('.2'))+4*mp.exp(mp.mpf('-.4')))/2,
                small_Z_eta_over_j=(Q(1, 100)+Q(5, 256))/(Q(5, 4)-Q(40, 256)),
                chi_separation=Q(6400, 6401),
                log_radial_comparison=Q(4), log_radial_error=Q(32),
                weighted_eta_comparison=Q(410), weighted_eta_error=Q(2000),
                pure_eta_log_bound=Q(500),
                source_linear_chi=Q(520), source_e_multiplier=Q(2160),
                source_remainder_prefactor=Q(2162, 10**6)+Q(500, 10**12)+Q(9, 10**24),
                p1_over_X_lower=Q(1, 10), p1_upper=Q(10),
                native_exit_v_lower=Q(23, 10), native_collar_v_lower=Q(9, 4),
                activated_v_lower=Q(17, 8), Kact_R_power=80,
                growth_lower_denominator=680)


def budget(scale='Md64', digits=100):
    """Outward checks of the selected source and activation inequalities."""
    if type(digits) is not int or digits < 100:
        raise ValueError('Nested amplitude logarithms require at least 100 digits')
    b = core.core_budget(scale, digits)
    mp.iv.dps = digits
    p, log, exp = mp.iv.mpf, mp.iv.log, mp.iv.exp
    upper = lambda v: mp.make_mpf(v._mpi_[1])
    lower = lambda v: mp.make_mpf(v._mpi_[0])

    def add(*values):
        shift = p(max(upper(v) for v in values))
        return shift+log(sum((exp(v-shift) for v in values), p(0)))

    lp, lk, ll, lm = (b['logs'][k] for k in ('P', 'K', 'Lambda', 'zeta_bound'))
    le, ls, lf = (b['logs'][k] for k in ('real_evaluation_bound', 'ball_norm_bound', 'map_displacement_numerator'))
    lU, lepsilon, lsigma = le+ls, le+lf-ll, log(p('1e-20'))
    source_error = add(log(p(2160))+lepsilon, log(p(2))+lU-ll,
                       log(p(500))+lU-ll-lsigma, log(p(9))+2*lU-ll-2*lsigma)
    lJ = add(p(0), ll, lk)
    lB = log(p('1e20'))+64*lJ
    B = exp(lB)
    logC = 10000*B-100*log(p('1e-16'))
    logR = B/100
    logKact = 2*logC+80*logR
    logt = -20*logC
    logs_c = -log(p(10))-log(logC)/2
    delta_log = log(p('1e-18'))+2*lp-log(p(100))
    axial_ratio = 2*logC+2*delta_log-log(p(80))-ll-4*exp(ll+lm)
    logs = dict(source_remainder=source_error,
                linear_chi_absorption=log(p(26000))-ll-log(p('.98')),
                normalized_shear_error_over_p001=log(p(64000))+lepsilon,
                axial_error_over_half_delta=log(p(4))+lepsilon-delta_log,
                axial_ratio_lower_log=axial_ratio,
                activation_K_over_C3=logKact-3*logC,
                width_times_R3=logt+3*logR,
                native_growth_change=logt+2*logC+8*logR,
                activation_algebra_error=logKact+logt,
                flat_activation_times_V=log(p(2))+p(16)/9-97*logC,
                squared_cone_loss_upper=log(p(32)/81)-31*logC+6*logs_c,
                inner_s=logs_c, inner_width=logt+logs_c,
                lambda_squared_lower_log=-4*exp(ll+lm)-log(p(680))-ll-2*logC,
                basic_envelope_log_gap=logR-20*exp(8*lJ))
    checks = dict(source_remainder_small=upper(source_error) < lower(log(p('.1'))),
                  linear_chi_absorption=upper(logs['linear_chi_absorption']) < 0,
                  endpoint_angular_error_small=upper(logs['normalized_shear_error_over_p001']) < 0,
                  axial_source_separation_survives_error=upper(logs['axial_error_over_half_delta']) < 0,
                  axial_exit_ratio_above_2p3=lower(axial_ratio) > upper(log(p('2.3'))),
                  activation_K_below_C3=upper(logs['activation_K_over_C3']) < 0,
                  small_logarithmic_field_change=upper(logs['width_times_R3']) < lower(log(p('.1'))),
                  native_growth_change_below_p05=upper(logs['native_growth_change']) < lower(log(p('.05'))),
                  activation_algebra_error_below_p01=upper(logs['activation_algebra_error']) < lower(log(p('.01'))),
                  flat_activation_growth_loss_below_one_eighth=upper(logs['flat_activation_times_V']) < lower(log(p(1)/8)),
                  squared_cone_loss_below_half=upper(logs['squared_cone_loss_upper']) < lower(log(p('.5'))),
                  collar_inside_native_core=upper(logt) < lower(log(log(p(41)/40))),
                  inner_collar_inside_flat_quarter=upper(logs_c) < lower(log(p('.25'))),
                  basic_derivative_envelope_fits_R=lower(logs['basic_envelope_log_gap']) > 0,
                  only_logarithms_of_C_R_and_width_evaluated=True)
    return dict(scale=scale, digits=digits, logs=logs, checks=checks)
