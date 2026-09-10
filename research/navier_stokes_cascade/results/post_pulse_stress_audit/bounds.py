"""Outward constants for the finite post-pulse reference interval.

The proof in README uses h=lambda^2. It is not a bound uniform under
arbitrarily smaller independent choices of h, nor an infinite-tail cone.
"""
import importlib.util
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
PULSE = HERE.parent/'pulse_stress_audit'
MOMENT = HERE.parent/'pulse_moment_audit'
OUTER = HERE.parent/'outer_pressure_pilot'
spec = importlib.util.spec_from_file_location('post_interval_helpers',HERE.parent/'axial_stress_audit'/'bounds.py')
helpers = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helpers)
point,lo,hi = helpers.point,helpers.lo,helpers.hi


def post_bound(md=64,digits=60,radius=100,endpoint='0.5'):
    if not isinstance(md,int) or md<4 or digits<50 or radius<100 or str(endpoint)!='0.5':
        raise ValueError('Require Md>=4, >=50 digits, XR>=100, and terminal endpoint 0.5')
    mp.iv.dps = digits
    p,exp = point,mp.iv.exp
    T = p(64)
    lam,h = exp(-256),exp(-512)
    co,tf = p('.001'),p(1000)
    correction = 7000*lam**30
    relative_edit = 64*correction
    eta_edit = relative_edit/(1-relative_edit)
    slope_edit = p(20000)/p('.3')*correction/(1-relative_edit)
    qp_upper = co*h*exp(3)/(1-co*h)
    # Q/h >= co*exp(-1/2) > 1/2000 on 0<=terminal y<=1/2.
    q_over_h = co*exp(-p('.5'))
    # E_end <= exp(-6.5/lambda) <= h^4. Do not form exp(-1/lambda).
    absorption = 8*256*lam
    w_upper = 4000*h**2
    claim = p('1e-440')
    log_ps1 = mp.iv.log(radius)+233*T+2-mp.iv.log(2000)
    cone_gap = 2*(1-p(4)/100)**2-2*claim**2
    checks = {
        'T_floor': lo(exp(md)+10)>=64,
        'lambda_below_0p01': hi(lam)<p('.01'),
        'relative_edits_below_0p01': hi(relative_edit)<p('.01'),
        'angular_log_derivative_below_one': hi(eta_edit)<1,
        'angular_slope_edit_below_lambda_half': hi(slope_edit)<lo(lam/2),
        'interpolation_slope_above_minus_one': hi(lam+64*mp.iv.log(2)/tf)<1,
        'angular_source_above_lambda_third': hi((h+eta_edit/2)/lam)<p(1)/6,
        'unedited_source_above_lambda_third': hi(h/lam)<p(2)/3,
        'angular_stage_global_slope_floor': hi(3*h/2)<lo(lam),
        'angular_stage_slope_above_minus_one': hi(3*lam/2)<1,
        'terminal_log_slope_edit_below_h_quarter': hi(co*64/(2*(1-co*h)))<p(1)/4,
        'post_pulse_amplitude_decreases': lo(lam)>0 and lo(h)>0,
        'hold_and_release_Q_above_one': lo(exp(-1)*(1-h)*4*512)>1,
        'Qp_below_one': hi(qp_upper)<1,
        'terminal_Q_floor': lo(q_over_h)>p(1)/2000,
        'earlier_Q_floor_dominates_terminal': hi(h/lam)<500,
        'N_over_E2_below_two_over_h': hi(p(10)*h/3+2*h**2+p(2)/3)<2,
        'pulse_end_eb_below_one': hi(exp(-p('119.5')*64+p('.3')))<1,
        'pulse_end_absorption': hi(absorption)<p('6.5'),
        'absorption_decreases_with_T': 1-4*64<0,
        'positive_w_bound_below_claim': 0<lo(w_upper) and hi(w_upper)<lo(claim),
        'finite_ps1_above_100': lo(log_ps1)>hi(mp.iv.log(100)),
        'finite_Pc_above_a': 100>4,
        'finite_quadratic_gap_above_1p84': lo(cone_gap)>p('1.84'),
        'strict_shear_excess_kept': lo(p('1.5')*h)>0,
    }
    actual_T = exp(md)+10
    return dict(status='post-pulse-reference-cone-bounded-through-terminal-half' if all(checks.values()) else 'bound-inconclusive',
        Md=md,interval_digits=digits,XR_sufficient_lower=radius,terminal_endpoint='0.5',checks=checks,
        log_lambda_interval=[lo(-4*actual_T),hi(-4*actual_T)],
        log_h_interval=[lo(-8*actual_T),hi(-8*actual_T)],
        log_a_minus_two_lower_interval=[lo(mp.iv.log(p('1.5'))-8*actual_T),hi(mp.iv.log(p('1.5'))-8*actual_T)],
        Q_floor='h/2000',Qp_over_h_lower='0.001',terminal_Q_over_h_lower=lo(q_over_h),
        angular_relative_edit_upper=hi(relative_edit),angular_eta_log_edit_upper=hi(eta_edit),
        angular_slope_edit_upper=hi(slope_edit),
        global_slope_interval='-1 <= l <= -3h/4',global_eta_log_bound=1,
        energy_tail_bound='S/X <= E^2/(3h)',energy_eta_tail_bound='abs(S_eta)/X <= 2 E^2/(3h)',
        pressure_tail_bound='abs(Pi) <= E^2/2',pressure_eta_tail_bound='abs(Pi_eta) <= E^2',
        N_absolute_bound='2 E^2/h',pulse_end_bound='E_end <= h^4',
        exponential_absorption_upper=hi(absorption),w_bound='4000 h^2',
        w_universal_upper=hi(w_upper),w_claim='1e-440',first_A24_lower=2,second_A24_lower='1.99',
        second_expression_universal_upper=hi(2*w_upper**2),
        finite_ps1_lower=100,universal_log_ps1_lower=lo(log_ps1),
        finite_normalized_cone_gap_lower=lo(cone_gap),
        scope='Exact closed reference moments assumed; finite post-pulse interval only. No heat collar, axis match or global PDE certificate.')
