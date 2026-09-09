"""Finite, outward bounds on both intermediate stages; see README derivation.

Only logarithms of the large parameters are needed. The constant 8100
comes from the positive-Q floor and the exact angular maximum, not fitting.
"""
import importlib.util
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
AXIAL = HERE.parent/'axial_stress_audit'
spec = importlib.util.spec_from_file_location('pinned_axial_bounds', AXIAL/'bounds.py')
axial = importlib.util.module_from_spec(spec)
spec.loader.exec_module(axial)
point, lo, hi = axial.point, axial.lo, axial.hi


def elementary_checks(digits=60):
    """Check every numerical constant in the continuum inequalities."""
    if digits < 50:
        raise ValueError('At least 50 interval digits are needed for the margin')
    mp.iv.dps = digits
    p = point
    e, delta = mp.iv.exp(1), p('2e-28')
    lam64, h64 = mp.iv.exp(-256), mp.iv.exp(-512)
    floor = 45/(128*e)
    return {
        'delta_floor': hi(mp.iv.exp(-64)) < lo(delta),
        'axial_Q_angle_floor': lo((1-217*delta)/2) > p(1)/10,
        'axial_Q_memory_floor': lo((1-217*delta)*floor) > p(1)/10,
        'positive_transport': lo(1-4*mp.iv.exp(-11)) > p(1)/2,
        'angular_source_floor': lo((1-2*h64)/2) > p(1)/3,
        'ramp_h_absorption': hi(20*e*mp.iv.exp(-7*64)) < 1,
        'power_h_absorption': hi(4*lam64) < 1,
        'initial_N_bound': hi(2*h64*(1+p('2.5')/64)+180*delta/64) < 1,
        'memory_floor_rounding': lo(1/(20*e)) > p(1)/60,
        'lambda_floor_rounding': lo(1/(4*e)) > p(1)/20,
        'finite_growth_budget': hi(lam64*(1+240*64)) < 1,
        'growth_budget_decreases': hi(240-p(4)*(1+240*64)) < 0,
        'stress_bound_decreases': hi(2*p(241)/(241*64+3)-2) < 0,
        'tail_power_bound': hi(8/e) < 13,
        'ramp_bound_below_power_bound': 900 < 8100,
    }


def decay_bound(md=64, digits=60):
    """All angles and both entire intermediate intervals, for the new family.

    T=exp(Md)+10; log P*=T+1; lambda=exp(-4T); h=lambda^2.
    The claimed universal upper bound uses T>=64, avoiding rounding 2+2lambda
    to 2. Its strict excess is represented separately by log(2lambda).
    """
    if not isinstance(md, int) or md < 4:
        raise ValueError('Md must be an integer >=4')
    if digits < 50:
        raise ValueError('At least 50 interval digits are required')
    mp.iv.dps = digits
    p = point
    T = mp.iv.exp(md)+10
    # H(T)=8100*e^3*(241T+3)^2*exp(-2T) decreases for T>=64.
    log_bound = mp.iv.log(8100)+3+2*mp.iv.log(241*T+3)-2*T
    universal = 8100*mp.iv.exp(3)*(241*64+3)**2*mp.iv.exp(-128)
    claimed = p('1e-42')
    if hi(universal) >= lo(claimed):
        raise ArithmeticError('The chosen outward claim is not justified')
    checks = elementary_checks(digits)
    return dict(
        Md=md, interval_digits=digits,
        Td_interval=[lo(T),hi(T)],
        log_lambda_interval=[lo(-4*T),hi(-4*T)],
        log_h_interval=[lo(-8*T),hi(-8*T)],
        log_a_minus_two_on_power_interval=[lo(mp.iv.log(2)-4*T),hi(mp.iv.log(2)-4*T)],
        Tw_interval=[lo(240*T),hi(240*T)],
        log_stress_upper_interval=[lo(log_bound),hi(log_bound)],
        universal_lambda_w_squared_upper=hi(universal),
        claimed_lambda_w_squared_upper='1e-42',
        first_A24_margin_lower=2,
        second_A24_margin_lower=lo(2-2*claimed),
        pole_tail_pressure_error_power=8,
        Q_ramp_lower='eta^2/10 + exp(-T-x)/20',
        Q_power_lower='eta^2/10 + lambda/20 + exp(-T-(1-lambda)*y)/60',
        status='intermediate-ratios-bounded' if all(checks.values()) else 'failed',
        elementary_checks=checks)


def finite_cone_bound(axial_result,digits=60):
    """A sufficient radius floor on the axial and two intermediate stages.

    This does not select the radius needed by the pending annulus or later
    corrections. Increasing XR preserves the inequalities established here.
    """
    mp.iv.dps=digits
    p=point
    XR=p(100)
    epsilon=p('4.34e-26')
    pa=XR*45/128*(1-epsilon)
    pr=XR*mp.iv.exp(1)/20
    pp=XR*mp.iv.exp(2)/60
    C=p(axial_result['maximum_bsw_absolute_upper'])
    D=p(axial_result['bs_squared_upper'])/2
    # At a=2, (vs-2)*(w+bs/2)^2=(bs*w+bs^2/2)^2/2.
    axial_gap=1-C-D/2-(2+D)/pa
    # At bs=0, the remaining test is alpha*w^2 < (1-a/ps1)^2.
    intermediate_gap=(1-p('2.2')/5)**2-p('1e-42')
    return dict(XR_sufficient_lower=100,
        axial_ps1_lower=lo(pa),ramp_ps1_lower=lo(pr),power_ps1_lower=lo(pp),
        axial_transformed_gap_lower=lo(axial_gap),
        intermediate_quadratic_gap_lower=lo(intermediate_gap),
        all_pass=lo(pa)>5 and lo(pr)>5 and lo(pp)>5 and lo(axial_gap)>0 and lo(intermediate_gap)>0,
        endpoint_scope='Relaxed cone when bs=0 on the axial stage or alpha=0 at the ramp start; strict admissible cone elsewhere on these stages.')
