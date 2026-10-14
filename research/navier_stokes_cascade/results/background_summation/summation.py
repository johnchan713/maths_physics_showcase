"""Cutoffs chosen in logarithmic coordinates; bounds remain conditional.

The scalar Riccati example is a diagnostic of asymptotic summation. It is
not a fluid field, a Navier-Stokes solver, or a blow-up candidate.
"""
from fractions import Fraction
import mpmath as mp


def cutoff_log_scale(index, gain, constraints, previous=0, boundary_log=0):
    """Choose s=log(a) so C*(1+|log q|)^P*q^gain <= 2^-index.

    constraints supplies finite (log(C), P) pairs. The estimate
    log(1+s)<=sqrt(s), followed by Young's inequality, controls ALL q
    below exp(-s), not just the endpoint. Actual coefficient constants
    would have to be supplied to evaluate actual Navier-Stokes scales.
    """
    gain, previous, boundary_log = map(mp.mpf, (gain, previous, boundary_log))
    if type(index) is not int or index < 1 or gain <= 0 or not constraints:
        raise ValueError('Positive index, gain and a nonempty finite constraint list required')
    if not all(mp.isfinite(x) for x in (gain, previous, boundary_log)) or min(previous, boundary_log) < 0:
        raise ValueError('Finite nonnegative previous and boundary logarithms required')
    candidates = [previous+mp.log(2), boundary_log+mp.log(2), mp.mpf(1)]
    for logC, power in constraints:
        logC, power = mp.mpf(logC), mp.mpf(power)
        if not all(mp.isfinite(x) for x in (logC, power)) or min(logC, power) < 0:
            raise ValueError('Pad C to at least one and use a nonnegative logarithmic power')
        candidates.append(2*(logC+index*mp.log(2))/gain+(power/gain)**2)
    # Fixed extra padding keeps extreme diagnostic scales away from rounding
    # at an equality. The continuum proof permits any larger real scale.
    return mp.mpf('1.000001')*(max(candidates)+1)


def loss_budget(h, derivative_order=0, target_power=0):
    """A conservative finite truncation sufficient for the comparison proof."""
    h = Fraction(h)
    if not 0 < h < Fraction(1, 100) or any(type(v) is not int or v < 0 for v in (derivative_order, target_power)):
        raise ValueError('Construction h, derivative order and decay target required')
    m, target, A = derivative_order, target_power, Fraction(1, 2)+h
    K = Fraction(3, 2)+2*h+m
    tail_threshold = (target+4*A+3*m+10)/h-1
    finite_threshold = (target+1+K)/(2*h)-1
    ceil = lambda x: -((-x.numerator)//x.denominator)
    J = max(m+2, ceil(tail_threshold), ceil(finite_threshold))
    return dict(h=h, derivative_order=m, target_power=target, comparison_order=J,
                input_tuple_loss=2*A+m, cutoff_tuple_loss=2*A+m+1,
                differential_polynomial_order=2, differential_polynomial_degree=2,
                stress_coefficient_loss=m+2,
                tail_exponent_at_m_plus_two=h*(J+1)-(2*A+m+3),
                finite_residual_exponent=2*h*(J+1)-K,
                check_finite_margin=2*h*(J+1)-K >= target+1,
                check_tail_margin=h*(J+1) >= target+4*A+3*m+10,
                actual_profile_order_or_scales_evaluated=False)


def scalar_coefficients(count):
    """Formal coefficients of q^2*f'-f+q+f^2=0: 1,2,8,44,..."""
    if type(count) is not int or count < 1:
        raise ValueError('A positive coefficient count is required')
    a = [0, 1]
    for n in range(2, count+1):
        a.append((n-1)*a[n-1]+sum(a[i]*a[n-i] for i in range(1, n)))
    return a


def scalar_residual_coefficients(coefficients, cutoffs, euler_derivatives):
    """Group equal q powers BEFORE evaluating, retaining cutoff derivatives.

    Directly subtracting exponentially tiny series values loses the
    cancellations to rounding; integer polynomial coefficients avoid that.
    euler_derivatives[n] is sigma_n*chi'(sigma_n), not an ordinary derivative.
    """
    N = len(coefficients)-1
    if len(cutoffs) != N+1 or len(euler_derivatives) != N+1:
        raise ValueError('One cutoff and its Euler derivative per coefficient required')
    rows = [mp.mpf(0)]*(2*N+1)
    rows[1] = 1
    for n in range(1, N+1):
        an = coefficients[n]
        rows[n] -= an*cutoffs[n]
        rows[n+1] += an*(n*cutoffs[n]+euler_derivatives[n])
        for j in range(1, N+1):
            rows[n+j] += an*coefficients[j]*cutoffs[n]*cutoffs[j]
    return rows


def proof_ledger():
    return dict(cutoffs='chi(a_n*q), chi=1 on [0,1/2], chi=0 on [1,infinity)',
                summation_applied_to_streamfunction_potentials=True,
                radial_commutator='r*delta_ur=-q^(2nh)*(2*eta/L)*(a_n*q)*chi_prime(a_n*q)*M_n',
                scales_stored_as_logarithms=True,
                scale_condition='a_(n+1)>=2*a_n, with each finite list of seminorms through order n',
                coefficient_bounds_may_grow_arbitrarily_fast=True,
                physical_tail='|U-U^[J]|_m<=2^(-J)*q^(h*(J+1)-(2*A+m+1))',
                residual_comparison_order_chosen_before_neighborhood=True,
                augmented_background_residual_flat_in_all_fixed_physical_derivatives=True,
                bare_background_momentum_residual_flat=False,
                normalized_background_velocity_error='O(q^(2h)) on each fixed enlarged compact annulus',
                stress_error='O(q^(2h)*zeta*delta^(-N_I)); only order one needs an edge inverse power',
                uniform_perturbed_stress_cone_not_inferred=True,
                same_heat_exterior_and_Imean_velocity_preserved=True,
                raw_infinite_power_series_convergence_not_required=True,
                infinite_background_swirl_lower_bound='.132*sqrt(2/Lambda)/C*q^(-A), for sufficiently small q',
                lower_bound_includes_physical_waves=False,
                complete_field_energy_or_smooth_force_not_supplied=True)
