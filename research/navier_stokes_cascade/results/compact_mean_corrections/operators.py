"""Compact mean inverses and their exact, scoped epsilon bookkeeping.

Numerical functions exercise manufactured sources. The continuum estimates
and their hypotheses are derived in README.md, rather than inferred from
these finite computations. All torus modes use exp(2*pi*i*k.y).
"""
from fractions import Fraction
from functools import lru_cache
import math
import mpmath as mp


def torus_vectors():
    beta = mp.sqrt(2)-1
    return (mp.mpf(1), -beta), (beta, mp.mpf(1))


def divisor(mode, temporal=False):
    if len(mode) != 2 or any(type(k) is not int for k in mode):
        raise ValueError('A two-dimensional integer Fourier mode is required')
    vector = torus_vectors()[int(temporal)]
    return sum(k*v for k, v in zip(mode, vector))


def temporal_inverse(modes, ci):
    """Solve ci*(vt.d_y)u=f with zero mean; reject an uncanceled zero mode."""
    ci = mp.mpf(ci)
    if not mp.isfinite(ci) or ci <= 0:
        raise ValueError('A positive finite fast-time coefficient is required')
    if modes.get((0, 0), 0) != 0:
        raise ValueError('The fast-time inverse cannot cancel the auxiliary mean')
    return {k: mp.mpc(v)/(2j*mp.pi*ci*divisor(k, True))
            for k, v in modes.items() if k != (0, 0)}


def polynomial_value(coefficients, x):
    answer = mp.mpf(0)
    for a in reversed(coefficients):
        answer = answer*x+a
    return answer


def polynomial_derivative(coefficients, order=1):
    answer = list(coefficients)
    for _ in range(order):
        answer = [i*a for i, a in enumerate(answer)][1:]
    return answer


def polynomial_bump(power=8, left=1, right=4):
    """(U-left)^power*(right-U)^power on [left,right], zero outside.

    This diagnostic has finite smoothness (C^(power-1)). The actual
    construction and the all-orders flatness argument require C-infinity.
    """
    a = [mp.mpf(math.comb(power, j))*(-left)**(power-j) for j in range(power+1)]
    b = [mp.mpf(math.comb(power, j))*right**(power-j)*(-1)**j for j in range(power+1)]
    out = [mp.mpf(0)]*(2*power+1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    # Unit integral makes the diagnostic's source scale easy to interpret.
    integral = sum(c*(mp.mpf(right)**(j+1)-mp.mpf(left)**(j+1))/(j+1)
                   for j, c in enumerate(out))
    return [c/integral for c in out]


def polynomial_oscillatory_integral(coefficients, left, right, omega):
    """Exact finite integration-by-parts formula, evaluated at high precision.

    No grid is asked to resolve a frequency proportional to a huge band.
    Cancellation for small omega still requires adequate caller precision.
    """
    left, right, omega = map(mp.mpf, (left, right, omega))
    if omega == 0:
        return sum(c*(right**(j+1)-left**(j+1))/(j+1)
                   for j, c in enumerate(coefficients))
    def endpoint(x):
        values, total = list(coefficients), mp.mpc(0)
        for j in range(len(coefficients)):
            total += (-1)**j*polynomial_value(values, x)/(1j*omega)**(j+1)
            values = polynomial_derivative(values)
        return mp.exp(1j*omega*x)*total
    return endpoint(right)-endpoint(left)


def radial_mode(radius, mode, scale, weight, cutoff, coefficients=None):
    """One mode of T_e for dr=2, shell 1<R<2, with the phase path kept.

    The input is f=R^-e*2R*F(R^2)*exp(2*pi*i*k.y).
    cutoff is chi(R); this returns its coefficient at y=0.
    """
    if weight not in (0, 1, 2):
        raise ValueError('The cylindrical weights are 0, 1 and 2')
    radius, scale = mp.mpf(radius), mp.mpf(scale)
    if radius <= 0 or scale <= 0 or not mp.isfinite(scale):
        raise ValueError('Positive radius and radial evaluation scale required')
    coeff = polynomial_bump() if coefficients is None else coefficients
    u = radius**2
    omega = 2*mp.pi*scale*divisor(mode)
    full = polynomial_oscillatory_integral(coeff, 1, 4, omega)
    upper = min(mp.mpf(4), max(mp.mpf(1), u))
    partial = polynomial_oscillatory_integral(coeff, 1, upper, omega)
    phase = mp.exp(-1j*omega*u)
    return radius**(-weight)*phase*(partial-cutoff(radius)*full)


@lru_cache(maxsize=None)
def bump_moment(power, center='1.12', width='.025', precision=70):
    """A fixed nonnegative C-infinity bump's power moment; a diagnostic only."""
    with mp.workdps(precision):
        p, c, w = map(mp.mpf, (power, center, width))
        def integrand(s):
            return (c+w*s)**p*mp.exp(-1/(1-s*s))*w if abs(s) < 1 else mp.mpf(0)
        return +mp.quad(integrand, [-1, 0, 1])


def moment_blocks(lam, amplitude, spacing='.12'):
    """The three angular rows and two axial rows of the fixed five-map.

    Columns are geometrically scaled bumps eta_j(x)=a_j^-1 eta_0(x/a_j).
    Row order: angular moment, P, Jz; axial flux moment, Jtheta.
    """
    lam, amplitude, spacing = map(mp.mpf, (lam, amplitude, spacing))
    if not all(mp.isfinite(x) for x in (lam, amplitude, spacing)) or lam <= 0 or amplitude == 0 or spacing <= 0:
        raise ValueError('Fixed lambda>0, nonzero patch amplitude and positive spacing required')
    angular_powers = [mp.mpf(2), -2-2*lam, -2*lam]
    axial_powers = [mp.mpf(1), 1-2*lam]
    def build(powers, factors):
        return mp.matrix([[factor*bump_moment(mp.nstr(p, mp.mp.dps))*mp.exp(j*spacing*p)
                           for j in range(len(powers))] for p, factor in zip(powers, factors)])
    return build(angular_powers, [1, 2*amplitude, -amplitude]), build(axial_powers, [1, amplitude])


def solve_moment_targets(pressure, angular_flux, axial_flux, lam='.2', amplitude='1.3'):
    """Linear inverse at the FIXED base, for arbitrary signed targets."""
    angular, axial = moment_blocks(lam, amplitude)
    u = mp.lu_solve(angular, mp.matrix([0, -pressure, -axial_flux]))
    s = mp.lu_solve(axial, mp.matrix([0, -angular_flux]))
    return u, s


def mean_update_budget(alpha=Fraction(9, 10), old=Fraction(9, 10), kappa=Fraction(1, 100000)):
    """Every radial-source term after the slow five-equation mean update.

    These estimates assume the current mean has tangential order old and
    radial order old+1. They do not certify a full wave/mean residual cycle.
    """
    alpha, old, kappa = map(Fraction, (alpha, old, kappa))
    if alpha < old or not 0 < old <= 1 or not 0 < kappa < Fraction(1, 10):
        raise ValueError('alpha>=old in (0,1] and small positive kappa required')
    terms = dict(slow_time=alpha+2,
                 base_radial_flux=alpha+2-kappa,
                 old_radial_flux=alpha+old+2-kappa,
                 new_radial_square=2*alpha+2-kappa,
                 base_axial_flux=alpha+2,
                 old_axial_flux=alpha+old+2,
                 new_axial_product=2*alpha+2,
                 old_swirl_product=alpha+old,
                 new_swirl_square=2*alpha,
                 radial_viscosity=alpha+2-2*kappa,
                 axial_viscosity=alpha+4)
    claimed = alpha+old-2*kappa
    return dict(radial_source_terms=terms, moment_defect_order=claimed,
                moment_defect_gain=claimed-alpha,
                temporal_mean_remainder=alpha+1-2*kappa,
                minimum_margin=min(terms.values())-claimed,
                full_residual_cycle_gain_verified=False)
