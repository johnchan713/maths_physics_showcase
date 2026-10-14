"""Continuum inequalities used by the conditional actual-band argument."""
from fractions import Fraction
import mpmath as mp


def eventual_band_threshold(h, constant=1, power=4):
    """For every ell>=ell0: C*ell^power*2^(-h*ell/2)<=1.

    log(x)<=x-1 with x=a*ell/(2*power) gives the displayed
    sufficient threshold. This is an all-later-bands argument, not a scan.
    A factor 3 can be included in constant for eps+eps^2+k^-1.
    """
    h, constant, power = map(mp.mpf, (h, constant, power))
    if not all(mp.isfinite(x) for x in (h, constant, power)) or not 0 < h < mp.mpf('.01') or constant < 1 or power <= 0:
        raise ValueError('Construction h, C>=1 and positive power required')
    a = h*mp.log(2)/2
    bracket = mp.log(constant)+power*mp.log(2*power/a)-power
    return int(mp.ceil(max(1, 2*bracket/a)))+1


def phase_majorant(ell, h):
    """S^2(eps+eps^2+k^-1)<=3*ell^4*sqrt(eps)."""
    ell, h = mp.mpf(ell), mp.mpf(h)
    return 3*ell**4*mp.power(2, -h*ell/2)


def barrier_certificate(mu, ce, cd, scale, length_over_scale):
    """Entrywise |E_ab|<=ce/S, |d-dref|<=cd/S; real m=1 pulse.

    At r=+rho the Riccati vector points left, and at -rho it
    points right. log(z_plus/P) is bounded on the entire interval.
    """
    mu, ce, cd, scale, ratio = map(mp.mpf, (mu, ce, cd, scale, length_over_scale))
    if not all(mp.isfinite(x) for x in (mu, ce, cd, scale, ratio)) or min(mu, ce, scale, ratio) <= 0 or cd < 0:
        raise ValueError('Positive spectral gap, entry bound, scale and length ratio required')
    eta = ce/scale
    if eta > mu/8:
        raise ValueError('The invariant-interval smallness condition is unproved')
    rho = eta/mu
    outward_upper = -2*mu*rho+eta+2*eta*rho+eta*rho*rho
    log_amplitude_bound = ratio*(ce+cd+ce*rho)
    return dict(rho=rho, outward_derivative_upper=outward_upper,
                log_growing_amplitude_bound=log_amplitude_bound,
                lower_radial_factor=(1-rho)*mp.exp(-log_amplitude_bound),
                upper_radial_factor=(1+rho)*mp.exp(log_amplitude_bound),
                checks=dict(invariant_interval_strict=outward_upper < 0,
                            positive_radial_lower_bound=rho < 1))


def derivative_exponents(order, coefficient_powers=None, source_powers=None):
    """Integer envelope for differentiated Duhamel at fixed pulse coordinate.

    A derivative of the coefficient of order i costs S^p_i;
    integration costs S once. Transverse/slow derivatives hold v fixed.
    With p_i=i the homogeneous exponent is 2*j.
    """
    if type(order) is not int or order < 0:
        raise ValueError('Nonnegative derivative order required')
    powers = list(range(order+1)) if coefficient_powers is None else list(coefficient_powers)
    if len(powers) != order+1 or any(type(p) is not int or p < 0 for p in powers):
        raise ValueError('One nonnegative integer coefficient power per order required')
    if source_powers is None:
        answer = [0]
        for j in range(1, order+1):
            answer.append(max(powers[i]+1+answer[j-i] for i in range(1, j+1)))
        return answer
    source = list(source_powers)
    if len(source) != order+1 or any(type(p) is not int or p < 0 for p in source):
        raise ValueError('One nonnegative integer source power per order required')
    forcing = [max(powers[i]+source[j-i] for i in range(j+1)) for j in range(order+1)]
    answer = [forcing[0]+1]
    for j in range(1, order+1):
        answer.append(max([forcing[j]+1]+[powers[i]+1+answer[j-i] for i in range(1, j+1)]))
    return answer


def cone_certificate(ac, u, margin, column_error):
    """Neumann-series bound for the two NORMALIZED covariance columns.

    The independent positive scalar weights h_sigma are factored out.
    The proof applies to all errors of norm<=column_error, not just samples.
    """
    ac, u, margin, error = map(mp.mpf, (ac, u, margin, column_error))
    if not all(mp.isfinite(x) for x in (ac, u, margin, error)) or min(ac, u, margin) <= 0 or margin > 1 or error < 0:
        raise ValueError('Positive directions, 0<margin<=1 and nonnegative error required')
    beta = error*(1/ac+1/u)
    if beta > margin/(4+margin):
        raise ValueError('Column error does not certify positive weights')
    return dict(beta=beta, relative_inverse_error=beta/(1-beta),
                minimum_weight_over_target_scale=margin/2-beta/(1-beta),
                minimum_singular_value_over_min_direction=mp.sqrt(2)*(1-error/min(ac, u)),
                checks=dict(positive_weights_survive=margin/2-beta/(1-beta) >= margin/4,
                            invertibility_survives=error < min(ac, u)))


def endpoint_flat_threshold(constant, decay, linear_loss, log_power, target):
    """C*ell^p*exp(B*ell-a*ell^2)<=2^(-target*ell) for all later ell."""
    c, a, b, p, target = map(mp.mpf, (constant, decay, linear_loss, log_power, target))
    if not all(mp.isfinite(x) for x in (c, a, b, p, target)) or c < 1 or a <= 0 or min(b, p, target) < 0:
        raise ValueError('Finite nonnegative losses, C>=1 and positive Gaussian decay required')
    return int(mp.ceil(max(1, 2*(p+b+target*mp.log(2))/a, mp.sqrt(2*mp.log(c)/a))))+1


def exponent_ledger(h=Fraction(1, 200), kappa=Fraction(1, 100000)):
    h, kappa = Fraction(h), Fraction(kappa)
    if not 0 < h < Fraction(1, 100) or not 0 < kappa < Fraction(1, 6):
        raise ValueError('Construction h and small positive radial loss required')
    primary, correction = Fraction(1, 2), 1-kappa
    covariance = 2*correction
    return dict(primary_amplitude=primary, exact_curl_correction=correction,
                real_primary_covariance_remainder=covariance,
                conditional_linear_wave_gain=Fraction(1, 2)-3*kappa,
                primary_linear_remainder=1-3*kappa,
                relative_normalized_stress_gain=covariance-1,
                physical_covariance_remainder_power=-1-2*h*kappa,
                relative_physical_stress_gain=h*(1-2*kappa),
                absolute_covariance_decay=False, full_residual_gain=False)
