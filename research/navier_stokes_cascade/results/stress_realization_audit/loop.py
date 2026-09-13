"""Appendix C loop algebra; numerical fixtures are not a global PDE profile.

The variance uses a positive series, not a difference of nearly equal Bessel
values. Interval evaluations enclose the infinite series with a geometric tail.
"""
from math import comb
import mpmath as mp


def lower(x):
    return mp.make_mpf(x._mpi_[0])


def upper(x):
    return mp.make_mpf(x._mpi_[1])


def point(x):
    return x if hasattr(x, '_mpi_') else mp.iv.mpf(str(x))


def me_minus_one(z):
    """Numerical positive series; preserve the small normalization correction."""
    z = mp.mpf(z)
    term = z*z/4
    total = term
    for k in range(1, 10000):
        term *= z*z/(4*(k+1)**2)
        total += term
        if abs(term) <= mp.eps*max(abs(total), mp.mpf('1e-100000')):
            return total
    raise ArithmeticError('Numerical Me series did not converge')


def positive_series_iv(z, variance_numerator=False):
    """Enclose Me(z), or [Me(2z)-Me(z)^2]/z^2, including z=0.

    The latter coefficient at z^(2k-2) is
    (1 - binom(2k,k)/4^k)/(k!)^2, k>=1. All coefficients are positive.
    The omitted terms are bounded by the factorial series with coefficient 1.
    """
    z = point(z)
    z2 = z*z
    total = point(0)
    term = point(1)
    start = 1 if variance_numerator else 0
    divisor = 1 if variance_numerator else 4
    tolerance = mp.mpf(10)**(-mp.iv.dps+10)
    for k in range(start, 10000):
        coefficient = (point(1)-point(comb(2*k,k))/4**k
                       if variance_numerator else point(1))
        total += coefficient*term
        next_term = term*z2/(divisor*(k+1)**2)
        ratio = z2/(divisor*(k+2)**2)
        if upper(ratio) < 1:
            tail = next_term/(1-ratio)
            if upper(tail) <= tolerance*max(mp.mpf(1), abs(lower(total))):
                return total + mp.iv.mpf([0, upper(tail)])
        term = next_term
    raise ArithmeticError('Interval series did not converge within the safe cap')


def variance_iv(mu, p, d0=1):
    mu, p, d0 = map(point, (mu, p, d0))
    if lower(mu) < 0 or lower(d0) <= 0:
        raise ValueError('Need nonnegative mu and positive d0')
    z = mu*p
    # Squaring an interval across zero requires the interval power operation.
    if lower(z) < 0 < upper(z):
        z = mp.iv.mpf([0, max(abs(lower(z)), abs(upper(z)))])
    return d0*d0*mu*mu*positive_series_iv(z, True)/positive_series_iv(z)**2


def variance(mu, p, d0=1):
    mu, p, d0 = map(mp.mpf, (mu, p, d0))
    if mu < 0 or d0 <= 0:
        raise ValueError('Need nonnegative mu and positive d0')
    z = mu*p
    term = mp.mpf(1)
    total = mp.mpf(0)
    for k in range(1, 10000):
        total += (1-mp.mpf(comb(2*k,k))/4**k)*term
        term *= z*z/(k+1)**2
        if abs(term) <= mp.eps*max(1, abs(total)):
            return d0*d0*mu*mu*total/(1+me_minus_one(z))**2
    raise ArithmeticError('Variance series did not converge')


def shear_ratio(theta, mu, p, ts=0, d0=1):
    theta, mu, p, ts, d0 = map(mp.mpf, (theta, mu, p, ts, d0))
    if p == 0:
        return ts+d0*mu*mp.sin(theta)
    remainder = me_minus_one(mu*p)
    return ts+d0*(mp.expm1(mu*p*mp.sin(theta))-remainder)/(p*(1+remainder))


def root_bracket(target, p, d0=1, cap=64, width='1e-36'):
    """Bracket the unique exact variance root; endpoints checked with intervals.

    Strict monotonicity is proved in README, not inferred from this bisection.
    """
    target, p, d0, cap, width = map(mp.mpf, (target, p, d0, cap, width))
    if target < 0 or d0 <= 0 or cap <= 0 or width <= 0:
        raise ValueError('Invalid root problem')
    if target == 0:
        return (mp.mpf(0), mp.mpf(0))
    if lower(variance_iv(cap, p, d0)) <= target:
        raise ValueError('Cap does not enclose the variance target')
    a, b = mp.mpf(0), cap
    while b-a > width:
        mid = (a+b)/2
        if variance(mid, p, d0) < target:
            a = mid
        else:
            b = mid
    if not (upper(variance_iv(a,p,d0)) < target < lower(variance_iv(b,p,d0))):
        raise ArithmeticError('Endpoint signs unresolved at the active precision')
    return a, b


def cone_gaps(a, b, p1, p2):
    a, b, p1, p2 = map(mp.mpf, (a,b,p1,p2))
    if a <= 0:
        raise ValueError('Cone coordinates require a>0')
    t = -b/a
    v, c, j = a*(1+t*t), p1+p2*t, p2-p1*t
    return [a, v-2, c-v, 2*(c-v)**2-(v-2)*j*j]


def state(theta, a, ts, p1, p2, mu, d0=1):
    a, ts, p1, p2, mu, d0 = map(mp.mpf, (a,ts,p1,p2,mu,d0))
    t = shear_ratio(theta,mu,p2,ts,d0)
    v = a*(1+ts*ts+variance(mu,p2,d0))
    al, bl = v/(1+t*t), -v*t/(1+t*t)
    density = a*(1+t*t)/v  # dphi / d(theta/(2*pi))
    return dict(t=t, v=v, a=al, b=bl, density=density,
                gaps=cone_gaps(al,bl,p1,p2))


def zero_pressure_primitives(theta, a='0.8', v='2.01', E=1):
    """Closed p2=ts=0 test: exact lifted phase and zero-mean antiderivatives.

    This is manufactured shear data, not the joined manuscript profile.
    """
    theta, a, v, E = map(mp.mpf, (theta,a,v,E))
    if not 0 < a < v or E <= 0:
        raise ValueError('Invalid primitive fixture')
    k = mp.sqrt(2*(v/a-1))
    t = k*mp.sin(theta)
    density = a*(1+t*t)/(2*mp.pi*v)
    phi = theta/(2*mp.pi)-a*k*k*mp.sin(2*theta)/(8*mp.pi*v)
    A = -a*a*k*k*mp.sin(2*theta)/(16*mp.pi*v)
    B = E*a*k*mp.cos(theta)/(4*mp.pi)
    return dict(phi=phi, phase_derivative=density, A=A, B=B,
                a=v/(1+t*t), b=-v*t/(1+t*t))
