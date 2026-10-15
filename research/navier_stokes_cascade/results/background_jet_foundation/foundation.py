"""Exact derivative conversions and finite cutoff budgets; not a field solver."""
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial
import mpmath as mp


def order(value):
    if type(value) is not int or value < 0:
        raise ValueError('A nonnegative integer derivative order is required')
    return value


@lru_cache(None)
def stirling_second_row(m):
    order(m)
    if m == 0:
        return (1,)
    previous = stirling_second_row(m-1)+(0,)
    return tuple((previous[j-1] if j else 0)+j*previous[j] for j in range(m+1))


@lru_cache(None)
def stirling_first_row(m):
    order(m)
    if m == 0:
        return (1,)
    previous = stirling_first_row(m-1)+(0,)
    return tuple((previous[j-1] if j else 0)-(m-1)*previous[j] for j in range(m+1))


def conversion_cost(m):
    order(m)
    return dict(frozen_to_euler=sum(stirling_second_row(m)), euler_to_frozen=factorial(m))


def falling(value, count):
    result = Fraction(1)
    for j in range(order(count)):
        result *= Fraction(value)-j
    return result


def shifted_cutoff_coefficients(lam, m):
    """Coefficients of sigma^j chi^(j) in (lam+sigma*d_sigma)^m chi."""
    lam, m = Fraction(lam), order(m)
    return tuple(sum(Fraction(comb(m, i))*lam**(m-i)*stirling_second_row(i)[j]
                     for i in range(j, m+1)) for j in range(m+1))


def normalized_jet_constant(h, n, m, profile_bound, cutoff_euler_bound):
    """README (10); the cutoff bound must extend through m+1."""
    h = Fraction(h)
    if type(n) is not int or n < 1 or not 0 < h < Fraction(1, 100):
        raise ValueError('Fixed construction h and positive coefficient order required')
    order(m)
    profile_bound, cutoff_euler_bound = map(Fraction, (profile_bound, cutoff_euler_bound))
    if min(profile_bound, cutoff_euler_bound) < 0 or cutoff_euler_bound < 1:
        raise ValueError('Nonnegative profile bound and padded cutoff bound required')
    return max(Fraction(1), 2*cutoff_euler_bound*(1+2*n*h)**m*profile_bound)


def cutoff_log_scale(h, n, log_constants, previous=0, boundary_log=0):
    """Finite list C_(n,0),...,C_(n,n), giving an ALL-smaller-q bound.

    Actual profile constants must be justified before this recipe is used
    for the fluid. Diagnostics supply only manufactured constants.
    """
    h = Fraction(h)
    if type(n) is not int or n < 1 or not 0 < h < Fraction(1, 100):
        raise ValueError('Fixed h in (0,1/100) and positive n required')
    if len(log_constants) != n+1:
        raise ValueError('Supply EVERY normalized derivative bound through order n')
    values = list(map(mp.mpf, log_constants))
    previous, boundary_log = map(mp.mpf, (previous, boundary_log))
    if not all(mp.isfinite(v) and v >= 0 for v in values+[previous, boundary_log]):
        raise ValueError('Finite nonnegative logarithms required')
    gain = mp.mpf(h.numerator)*n/h.denominator
    return 1+max(previous+mp.log(2), boundary_log+mp.log(2),
                 max((v+n*mp.log(2))/gain for v in values))


def axis_average(coefficients):
    """Exact X coefficients of M/X for U=sum u_j X^j, with regular axis value."""
    return tuple(Fraction(c)/(j+1) for j, c in enumerate(coefficients))


def picard_log_term(k, z, strip_loss):
    order(k)
    z, strip_loss = map(mp.mpf, (z, strip_loss))
    if not mp.isfinite(z) or z <= 0 or not mp.isfinite(strip_loss) or not 0 < strip_loss <= 1:
        raise ValueError('Finite positive z and strip loss in (0,1] required')
    p = (k+1)//2
    return (k+1)*mp.log(z)-mp.loggamma(k+2)+(p*mp.log(max(1, p/strip_loss)) if p else 0)


def band_coordinates(R, Z, T, h, seed=1):
    """Diagnostic implicit coordinate solve, keeping h fixed."""
    R, Z, T, h, seed = map(mp.mpf, (R, Z, T, h, seed))
    s = mp.findroot(lambda v: v-Z*Z*v**(2*h)-T, seed)
    return s, R*R/(2*s), Z*s**(-mp.mpf('.5')+h)


def band_coefficients(R, s, eta, h):
    """Rows of (12), in the basis (q*d_q, d_X, d_eta)."""
    R, s, eta, h = map(mp.mpf, (R, s, eta, h))
    D, X, L = mp.mpf('.5')-h, R*R/(2*s), 1-2*h*eta*eta
    return ((mp.mpf(0), R/s, mp.mpf(0)),
            (2*eta*s**(-D)/L, -2*eta*X*s**(-D)/L, (1-eta*eta)*s**(-D)/L),
            (1/(s*L), -X/(s*L), -D*eta/(s*L)))
