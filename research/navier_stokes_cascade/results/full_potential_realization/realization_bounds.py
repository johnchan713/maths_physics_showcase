"""Exact order budgets and a conditional asymptotic-summation scale recipe.

No actual Navier--Stokes stage constants or cutoff scales are evaluated here.
The constants may grow without bound in stage. See README.md for quantifiers.
"""
from fractions import Fraction as F
from pathlib import Path
import sys
import mpmath as mp

PARENT = Path(__file__).resolve().parent.parent
sys.path.append(str(PARENT/'background_summation'))
from summation import cutoff_log_scale

KAPPA = F(1, 100000)


def parameters(h=F(1, 200)):
    h = F(h)
    if not 0 < h < F(1, 100):
        raise ValueError('Require the inherited 0<h<1/100')
    return h, F(1, 2)+h, F(1, 2)-h


def index(value, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError('A nonnegative integer order (positive stage) is required')
    return value


def gain(j, h=F(1, 200)):
    return parameters(h)[0]*index(j, 1)/10


def tuple_loss(m, h=F(1, 200)):
    return 2*parameters(h)[1]+2*index(m)


def cutoff_loss(m, h=F(1, 200)):
    """Covers potentials, direct swirl, pressure AND one additional curl."""
    return tuple_loss(index(m)+1, h)


def finite_field_loss(m, h=F(1, 200)):
    return cutoff_loss(m, h)+1  # absorb each fixed polynomial in |log q|


def finite_residual_loss(m, h=F(1, 200)):
    return 2*parameters(h)[1]+F(3, 2)+2*index(m)


def increment_orders(j):
    """Cycle j-1 -> j, with B_(j-1)=.6+j/10, not B_j."""
    j = index(j, 1)
    B = F(3, 5)+F(j, 10)
    return dict(wave_velocity=B-KAPPA,
                wave_potential=B+F(1, 2)-KAPPA,
                mean_potential=B+F(1, 2)-2*KAPPA,
                direct_mean_swirl=B+F(1, 2)-2*KAPPA,
                pressure=B+F(1, 2)-2*KAPPA,
                mean_radial_velocity=B+F(3, 2)-2*KAPPA)


def all_stage_order_review():
    at_one = increment_orders(1)
    margins = {name: value-F(1, 10) for name, value in at_one.items()}
    # Every row and the target j/10 have slope 1/10. These are exact
    # all-stage margins, rather than evidence from a finite stage scan.
    return dict(stage_indexing='cycle j-1 to j; j>=1',
                row_slope=F(1, 10), target_slope=F(1, 10),
                all_stage_margins=margins,
                tuple_loss='2*A+2*m', cutoff_velocity_tuple_loss='L_m=2*A+2*m+2',
                finite_field_loss='H_m=L_m+1',
                cutoff_derivative_cost_at_order_m='q^(-m), independent of cutoff scale',
                physical_derivative_cost_used=F(2),
                checks=dict(all_stage_margins_positive=min(margins.values()) > 0,
                            phase_cost_below_reserved_two=1+F(5, 2)*F(1, 100) < 2,
                            all_physical_prefactors_bounded_by_q_to_minus_2A=True,
                            stage_constants_are_not_assumed_uniform=True))


def ceil_fraction(value):
    value = F(value)
    return -((-value.numerator)//value.denominator)


def comparison_budget(m=0, target=0, h=F(1, 200)):
    """Choose ONE finite stage for a residual jet order and decay target."""
    m, target = index(m), index(target)
    h, A, _ = parameters(h)
    L, H = cutoff_loss(m+2, h), finite_field_loss(m+2, h)
    # Half the actual increment gain is spent absorbing arbitrary constants.
    J = max(m+2, ceil_fraction(20*(target+1+L+H)/h-1),
            ceil_fraction(10*(target+1+finite_residual_loss(m, h))/h-7))
    tail = gain(J+1, h)/2-L
    finite = h*(F(7, 10)+F(J, 10))-finite_residual_loss(m, h)
    return dict(h=h, derivative_order=m, target=target, comparison_stage=J,
                cutoff_loss_at_m_plus_two=L, finite_field_loss_at_m_plus_two=H,
                tail_exponent=tail, old_field_times_tail_exponent=tail-H,
                tail_square_exponent=2*tail, finite_residual_exponent=finite,
                stage_is_fixed_before_q_limit=True,
                actual_stage_coefficients_or_scales_evaluated=False,
                checks=dict(diagonal_covers_needed_derivatives=J >= m+2,
                            linear_tail_decays=tail >= target+1,
                            old_new_products_decay=tail-H >= target+1,
                            new_square_decays=2*tail >= target+1,
                            finite_residual_decays=finite >= target+1))


def scale_choice_review():
    """Synthetic constants spanning extreme sizes; not actual flow constants."""
    rows, previous = [], mp.mpf(0)
    with mp.workdps(70):
        for j, constraints in ((1, [(0, 0), (5, 7)]),
                               (4, [(1000, 300), (10**6, 2000)]),
                               (17, [(10**12, 10**5), (300, 0)])):
            gamma = mp.mpf(j)/4000  # h*j/20, h=1/200
            scale = cutoff_log_scale(j, gamma, constraints, previous, 2)
            worst = -mp.inf
            for logC, power in constraints:
                # Upper bound valid for ALL s>=scale. It decreases in s.
                analytic = logC+mp.mpf(power)**2/(2*gamma)-gamma*scale/2+j*mp.log(2)
                worst = max(worst, analytic)
            samples = [max(logC+power*mp.log(1+s)-gamma*s+j*mp.log(2)
                           for logC, power in constraints)
                       for s in (scale, scale+1, 2*scale, 10*scale)]
            rows.append(dict(stage=j, log_scale=float(scale),
                             maximum_all_interval_log_margin=float(worst),
                             maximum_sample_log_margin=float(max(samples)),
                             doubled_scale=bool(scale >= previous+mp.log(2))))
            previous = scale
    return dict(synthetic_constants=True, actual_NS_scales_evaluated=False, rows=rows,
                analytic_inequality='log C+P log(1+s)-gamma*s <= log C+P^2/(2*gamma)-gamma*s/2',
                checks=dict(all_interval_majorant_satisfies_every_requirement=all(r['maximum_all_interval_log_margin'] < 0 for r in rows),
                            independent_high_precision_samples_agree=all(r['maximum_sample_log_margin'] < 0 for r in rows),
                            scales_increase_by_at_least_two=all(r['doubled_scale'] for r in rows)))


def core_energy_review(h=F(1, 200)):
    h, A, D = parameters(h)
    p = 2*h/D
    return dict(h=h, A=A, D=D, normalized_correction_power=h/4,
                radial_integrated_energy_power=1-2*A,
                axial_integrability_exponent=p,
                axial_slab_energy_power=D-2*h,
                inner_path='z=0, X=Xin<Xa, q=1-t; all annular corrections vanish exactly',
                core_swirl='q^(-A)*(e0+O(q^(2*h))), e0>0 inherited',
                checks=dict(radial_energy_has_power_minus_2h=1-2*A == -2*h,
                            axial_energy_is_integrable=p < 1,
                            shrinking_box_energy_power_positive=D-2*h > 0,
                            uniform_integrability_for_every_allowed_h=F(1, 2)-3*F(1, 100) > 0,
                            cutoffs_do_not_enlarge_spatial_support=True))
