"""Actual unmodulated C3 transfer, conditional on frozen core/reference bounds.

Large parameters are logarithms. This does not bound the post-modulation C3
moment: that needs a further derivative of the loop's input coordinates.
"""
from fractions import Fraction as F
import importlib.util
from pathlib import Path
import sys
import mpmath as mp

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent


def module_at(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


core = module_at('first_correction_core_bounds', HERE.parent/'axis_core_attachment'/'bounds.py')
p, lower, upper, logsum = core.point, core.lower, core.upper, core.logsum


def angular_transfer(digits=80):
    old = core.certificate('Md64', digits)
    b = old['log_bounds']
    log, exp = mp.iv.log, mp.iv.exp
    B, lc, eps = exp(b['B']), b['C'], p('1e-16')
    evaluation = sum((F(4, 3)**(k+1)/((k+1)**2) for k in range(4)), F(0))
    ld3 = log(3)+b['ball_norm_bound']-3*b['rho']
    lstar = b['Lambda']+b['zeta_bound']+logsum(log(2), -b['r0'], log(p(4)/3)-2*b['r0'])
    llogphi = log(100)+3*ld3
    L0 = log(p('27.5'))+b['Lambda']
    continuation = logsum(lstar, llogphi, log(p('55e6'))+b['Lambda']+8*b['K'], log(2*L0+10))
    core_amplitude_log = log(2)/2+exp(lstar)+ld3
    native = ld3-b['Lambda']
    axial_change = log(p('1e8'))+8*b['K']+b['activation_width']
    row_log = log(p('2e6'))+412*B-2*lc
    gb = p('3e-18')
    target = p('.255')*eps+gb+gb**2
    coefficient = 2000*target
    zeta = 24*exp(p('1.2'))/256
    checks = dict(inherited_parameters=all(old['checks'].values()),
                  exact_C3_generating_function=evaluation == F(544, 243) and evaluation < 3,
                  native_C3_offset=upper(native) < lower(log(eps/100)),
                  actual_continuation_C3_offset=upper(axial_change) < lower(log(eps/100)),
                  core_logarithm_and_axis_factor=upper(core_amplitude_log-B) < 0,
                  continuation_logarithm=upper(continuation-b['B']) < 0,
                  shape_log_f_C3=upper(201*B+170-204*B) < 0,
                  field_C3_amplitude=upper(204*B-lc) < 0,
                  actual_incoming_rows=upper(row_log) < lower(log(eps/4)),
                  C3_normalization_quadratic=upper(zeta) < 1,
                  same_root_C3_contraction=upper(8*1000**2*100*target) < 1,
                  joining_coefficient_C3=upper(coefficient) < p('5.72e-14'),
                  raw_third_factor_retained=upper(6*coefficient) < p('3.432e-13'))
    return dict(interval_digits=digits, exact_evaluation_constant=str(evaluation),
                log_D3=ld3, log_axis_logarithm_bound=lstar,
                log_continuation_logarithm_bound=continuation,
                log_actual_incoming_row_bound=row_log,
                terminal_offset_C3_bound=gb, target_C3_bound=target,
                joining_coefficient_C3_bound=coefficient,
                joining_raw_third_bound=6*coefficient, zeta_C3_bound=zeta,
                actual_unmodulated_incoming_five_moments_C3_transfer_supplied=True,
                global_third_order_compact_envelope_supplied=False,
                post_modulation_C3_bound_supplied=False, checks=checks)


def norm_ledger():
    return dict(f_inverse_C3=F(5), f_inverse_square_C3=F(24),
                log_f_C3_Cauchy_bound=F(170),
                energy_row_coefficient=24*(101+16*10+64),
                original_moment_angular_order=3,
                integrated_stress_angular_order=2,
                required_loop_angular_order_for_modulated_C3=3,
                required_original_moment_order_for_that_loop=4,
                checks=dict(energy_row_below_old_generous_cap=24*(101+16*10+64) < 10**6,
                            normalization_cost_increases=24 > 20,
                            loop_consumes_one_more_moment_derivative=4 == 3+1))


def inner_operator_certificate(digits=80):
    """Actual core envelope on 0<=X<=3/Lambda, conditional on frozen bounds.

    Outer eta neighborhood rho/4 gives |Phi|,|u|<=2S at |Y|<=4.
    Half that neighborhood supplies up to three angular and two radial
    derivatives. No modulation, activation cutoff or numerical core is used.
    """
    old = core.certificate('Md64', digits)
    b = old['log_bounds']
    log, exp = mp.iv.log, mp.iv.exp
    S, Lambda, rho = exp(b['ball_norm_bound']), exp(b['Lambda']), exp(b['rho'])
    jet = log(p('1e10'))+log(1+S)+2*log(1+Lambda)+3*log(1+1/rho)
    lc = b['C']
    # Complex eta has |eta|<1.01; use h<.01 and retain the denominator.
    pole_factor = p('.02')*p('1.01')**2
    return dict(interval_digits=digits,
                actual_interval='0<=X<=3/Lambda=(3/4)*Xa',
                outer_angular_neighborhood='dist(eta,[-1,1])<rho/4',
                operator_angular_neighborhood='dist(eta,[-1,1])<rho/8',
                jet_majorant='1e10*(1+S)*(1+Lambda)^2*(1+rho^(-1))^3',
                log_jet_majorant=jet, core_jet_envelope='C^4',
                matrix_source_majorant='1e6*K_core^2',
                actual_operator_envelope='C^16',
                actual_inner_interval_numerically_resolved=False,
                checks=dict(inherited_core_parameters=all(old['checks'].values()),
                            inner_interval_inside_exact_core=F(0) < F(3, 4) < F(1),
                            complex_core_generating_function=F(1, 1)/(1-F(1, 5)-F(1, 4)) < 2,
                            no_L_pole_on_neighborhood=upper(pole_factor) < p('.03'),
                            Cauchy_jet_below_C4=upper(jet-4*lc) < 0,
                            operator_bound_below_C16=upper(log(p('1e6'))+8*lc-16*lc) < 0,
                            inner_radius_less_than_one=upper(log(3)-b['Lambda']) < 0))


def third_implicit_terms():
    """Exact scalar differentiation ledger for Bc+Q(eta)[c,c]=z."""
    return dict(equation='L c3 = z3 - 6 Q[c1,c2] - 6 Q1[c,c2] - 6 Q1[c1,c1] - 6 Q2[c,c1] - Q3[c,c]',
                coefficients=[6, 6, 6, 6, 1],
                L='B+2Q[c, dot]',
                # Factorial norm controls the raw third derivative with 3!.
                raw_derivative_factor=6)
