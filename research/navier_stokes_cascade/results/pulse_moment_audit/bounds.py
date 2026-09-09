"""Outward constants for the pulse moment argument, not its stress cone.

The derivation in README.md treats the entire angular interval. Exponentials
of 1/lambda are eliminated analytically before any large parameter is used.
Positive Riemann enclosures, not quadrature error estimates, bound K_b.
"""
from fractions import Fraction
import importlib.util
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
AXIAL = HERE.parent/'axial_stress_audit'
spec = importlib.util.spec_from_file_location('pulse_interval_helpers', AXIAL/'bounds.py')
helpers = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helpers)
point, lo, hi, step_at = helpers.point, helpers.lo, helpers.hi, helpers.step_at


def shape_enclosure(boxes=1024, digits=60):
    """Enclose integral exp(-2xi) R0(xi)^2 dxi without sampling inference.

    Phi(.02*s)=.02*integral_0^s sigma. Monotone left/right sums enclose
    this primitive. On [.02,10], symmetry makes Phi(xi)=xi-.01 exactly.
    On [10,11], positive polynomial-exponential weights are integrated
    exactly and multiplied by endpoint bounds for the decreasing step.
    """
    if not isinstance(boxes, int) or boxes < 16 or digits < 40:
        raise ValueError('At least 16 boxes and 40 interval digits are required')
    mp.iv.dps = digits
    p = point
    width, offset = p('.02'), p('.01')
    values = [step_at(Fraction(i, boxes)) for i in range(boxes+1)]
    primitive_lo, primitive_hi = p(0), p(0)
    start_lo, start_hi, cut_lo, cut_hi = [p(0) for _ in range(4)]
    F = lambda x: -mp.iv.exp(-2*x)*((x-offset)**2/2+(x-offset)/2+p('.25'))
    for i in range(boxes):
        a, b = p(Fraction(i, boxes)), p(Fraction(i+1, boxes))
        primitive_hi += values[i+1]/boxes
        mass = (mp.iv.exp(-2*width*a)-mp.iv.exp(-2*width*b))/(2*width)
        start_lo += width**3*mass*primitive_lo**2
        start_hi += width**3*mass*primitive_hi**2
        primitive_lo += values[i]/boxes
        # Squared step is decreasing; the remaining positive weight is exact.
        weight = F(10+b)-F(10+a)
        cut_lo += weight*values[boxes-i-1]**2
        cut_hi += weight*values[boxes-i]**2
    central = F(p(10))-F(width)
    lower = lo(central+start_lo+cut_lo)
    upper = hi(central+start_hi+cut_hi)
    return dict(boxes=boxes, interval_digits=digits, K_lower=lower, K_upper=upper,
                prefix_lower=lo(start_lo), prefix_upper=hi(start_hi),
                cutoff_lower=lo(cut_lo), cutoff_upper=hi(cut_hi))


def closure_bound(boxes=1024, digits=60, md=64):
    """Finite moment closure for the selected reference schedule.

    The smallness comparisons use T>=64. The actual Md only appears in
    logarithmic parameter displays. No numerical Md=4 fixture is promoted
    to the all-angle Md=64 construction.
    """
    if not isinstance(md, int) or md < 4 or digits < 40:
        raise ValueError('Expected integer Md>=4 and at least 40 interval digits')
    shape = shape_enclosure(boxes, digits)
    mp.iv.dps = digits
    p, exp = point, mp.iv.exp
    lam = exp(-256)
    h = lam**2
    T0, tf, co = p(64), p(1000), p('.001')
    # The first two derivatives use coarse analytic bounds, including tails.
    step_bound, second_bound = p(64), p(20000)
    flat_first = (2*p(4)**3+16)*exp(-12)
    flat_second = (4*p(4)**6+6*p(4)**4+64*p(4)**3+352)*exp(-12)
    # Integral of the width-.3 bump is .3, by the fundamental theorem.
    bmin, bmax = p('.3')*exp(p('-.075')), p('.3')*exp(p('.075'))
    axial_inverse_times_lambda = 2*bmax*exp(1)/(bmin**2*exp(p('.98')))
    # I and pressure rows have limiting slopes +1 and -1 respectively.
    i_min, i_max = p('.3')*exp(p('-.15')), p('.3')*exp(p('.15'))
    c_min, c_max = p('.6')*exp(p('-.153')), p('.6')*exp(p('.153'))
    determinant = i_min*c_min*(exp(p('1.98'))-exp(-2))
    angular_inverse = max(hi((c_max*exp(-2)+i_max*exp(2))/determinant),
                          hi((c_max+i_max)/determinant))
    quadratic_norm = p('.3')*step_bound*exp(p('.153'))*(1+exp(-2))
    angular_c = 7000*lam**30
    relative_edit = step_bound*angular_c
    angular_slope_edit = second_bound/p('.3')*angular_c/(1-relative_edit)
    contraction = 8*5**2*32*700*lam**30
    # Absorb exp(-.8/lambda)/lambda into lambda^20, preserving a finite bound.
    exponential_absorption = lam*(mp.iv.log(4000)+21*256)
    remainder = 1000*T0*lam
    amplitude_error = 80*lam**41
    eta_error = lam*(1+512*T0*exp(-26))+80*lam**41
    claimed = p('1e-100')
    C = (1-exp(-26))/4
    left, right = p('1.0100502'), p('1.0100504')
    left_residual = left**2*p(shape['K_upper'])-C+remainder
    right_residual = right**2*p(shape['K_lower'])-C-remainder
    derivative = p('1.8')*p(shape['K_lower'])-amplitude_error
    checks = {
        'T_floor': lo(exp(md)+10) >= 64,
        'lambda_below_0p01': hi(lam) < p('.01'),
        'step_first_derivative_bound': hi(flat_first) < 64 and (256/4) <= 64,
        'step_second_derivative_bound': hi(flat_second) < 20000 and (256**2+3072)/4 < 20000,
        'interpolation_slope_bound': hi(step_bound*mp.iv.log(2)/tf) < p('.1'),
        'terminal_slope_bound': hi(co*step_bound/(2*(1-co*h))) < p('.25'),
        'prefix_m_C1_below_one': hi(24*exp(-p('120.5')*64+p('.7')))<1,
        'prefix_U2_below_one': hi(64*exp(-64+p('.4')))<1,
        'prefix_S_eta_below_one': hi(256*exp(-64+p('.4')))<1,
        'prefix_energy_growth_below_two': hi(2*lam*(240*T0+p('.5'))) < lo(mp.iv.log(2)),
        'prefix_S_below_500T': 241*64+p('3.5') < 500*64,
        'M_J_inverse_bound': hi(axial_inverse_times_lambda) < 20,
        'main_pulse_rhs_bound': hi(exp(p('1.5'))*(1+p('1.2')*11/p('.4'))) < 200,
        'axial_bumps_below_lambda20': hi(exponential_absorption) < p('.8'),
        'angular_discrepancy_bound': hi(10/p('.99')*exp(4)) < 700,
        'angular_inverse_bound': angular_inverse < 5,
        'angular_quadratic_bound': hi(quadratic_norm) < 32,
        'angular_contraction': hi(contraction)<1,
        'angular_edit_below_0p01': hi(relative_edit)<p('.01'),
        'angular_slope_stays_negative': hi(angular_slope_edit)<lo(lam/2),
        'post_density_eta_bound': hi(2+2*relative_edit/(1-relative_edit))<4,
        'post_mass_below_256T': 1002+240*64 <= 256*64,
        'finite_positive_tail_wait': lo(exp(-1)*(1-h)*4*mp.iv.log(1/h))>1 and hi(co*h/(1-co*h)*exp(3))<1,
        'axial_bump_energy_bound': p('.3')*step_bound < 20,
        'S_remainder_below_1000Tlambda': hi(lam*(500*T0+128*T0*exp(-26))+40*lam**41)<lo(remainder),
        'S_eta_below_1000Tlambda': hi(eta_error)<lo(remainder),
        'uniform_remainder_below_claim': hi(remainder)<lo(claimed),
        'amplitude_derivative_error_below_claim': hi(amplitude_error)<lo(claimed),
        'Tlambda_decreases': 1-4*64<0,
        'K_inside_paper_bounds': p('.20')<shape['K_lower']<shape['K_upper']<p('.25'),
        'left_bracket_negative': hi(left_residual)<0,
        'right_bracket_positive': lo(right_residual)>0,
        'unique_amplitude_root': lo(derivative)>p('.35'),
    }
    actual_T = exp(md)+10
    return dict(status='reference-moments-closed-pulse-cone-unverified' if all(checks.values()) else 'bound-inconclusive',
        shape=shape, checks=checks, Md=md, Tf=1000, co='0.001',
        Td_interval=[lo(actual_T),hi(actual_T)],
        log_lambda_interval=[lo(-4*actual_T),hi(-4*actual_T)],
        log_h_interval=[lo(-8*actual_T),hi(-8*actual_T)],
        axial_inverse_times_lambda_upper=hi(axial_inverse_times_lambda),
        angular_inverse_upper=angular_inverse,
        angular_contraction_upper=hi(contraction),
        angular_relative_edit_upper=hi(relative_edit),
        axial_c_C1_bound='lambda^20 at fixed Amp; 2 lambda^20 after substituting the root',
        angular_c_C1_bound='7000 lambda^30',
        post_energy_mass_bound='256 T',
        remainder_bound='1000 T lambda',
        remainder_universal_upper=hi(remainder),
        claimed_remainder_upper='1e-100',
        S_amplitude_derivative_error_upper=hi(amplitude_error),
        amplitude_lower='1.0100502', amplitude_upper='1.0100504',
        left_residual_upper=hi(left_residual),right_residual_lower=lo(right_residual),
        S_amplitude_derivative_lower=lo(derivative),
        amplitude_eta_derivative_bound='3000 T lambda',
        omitted_axial_bump_amplitude_error_bound='120 lambda^41',
        scope='Exact reference moment closure by an analytic reduction and outward constants; no pulse cone or global PDE certificate.')
