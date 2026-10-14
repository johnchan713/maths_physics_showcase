"""Regular scalar loop and raw third-derivative formulas for diagnostics.

Only finite manufactured parameters are evaluated here. Actual profile bounds
are the analytic estimates in higher_bounds.py and README.
"""
from math import comb
import mpmath as mp


def regular_variance_factor(z):
    """R(z)=[I0(2z)-I0(z)^2]/[z^2 I0(z)^2], regular at z=0.

    Evaluate the numerator divided by z^2 using its positive series; never
    subtract nearly equal Bessel values or divide by a small pressure.
    """
    z = mp.mpf(z)
    term, numerator = mp.mpf(1), mp.mpf(0)
    for k in range(1, 10000):
        numerator += (1-mp.mpf(comb(2*k, k))/4**k)*term
        term *= z*z/(k+1)**2
        if abs(term) < mp.eps*max(1, abs(numerator)):
            return numerator/mp.besseli(0, z)**2
    raise ArithmeticError('Diagnostic positive series failed to converge')


def signed_F(mu, pressure, d0):
    return d0*mu*mp.sqrt(regular_variance_factor(mu*pressure))


def root_jet(target, pressure, d0, x0):
    """Solve F(mu,p)=s, then differentiate it without dividing by V_mu."""
    s = [mp.diff(target, x0, k) for k in range(4)]
    p = [mp.diff(pressure, x0, k) for k in range(4)]
    f = lambda mu, q: signed_F(mu, q, d0)
    mu = mp.findroot(lambda value: f(value, p[0])-s[0], s[0]*mp.sqrt(2)/d0)
    derivative = lambda i, j: mp.diff(f, (mu, p[0]), (i, j))
    fm, fp = derivative(1, 0), derivative(0, 1)
    fmm, fmp, fpp = derivative(2, 0), derivative(1, 1), derivative(0, 2)
    m1 = (s[1]-fp*p[1])/fm
    m2 = (s[2]-fp*p[2]-fmm*m1*m1-2*fmp*m1*p[1]-fpp*p[1]**2)/fm
    cubic = (derivative(3, 0)*m1**3+3*derivative(2, 1)*m1*m1*p[1]
             +3*derivative(1, 2)*m1*p[1]**2+derivative(0, 3)*p[1]**3)
    cross = 3*(fmm*m1*m2+fmp*(m2*p[1]+m1*p[2])+fpp*p[1]*p[2])
    m3 = (s[3]-fp*p[3]-cubic-cross)/fm
    return [mu, m1, m2, m3]


def implicit_phase_jet(forward, x0, theta0):
    """Third raw jet of Theta(x,phi0), with phi0 fixed."""
    D = lambda i, j: mp.diff(forward, (x0, theta0), (i, j))
    w = D(0, 1)
    t1 = -D(1, 0)/w
    t2 = -(D(2, 0)+2*D(1, 1)*t1+D(0, 2)*t1*t1)/w
    t3 = -(D(3, 0)+3*D(2, 1)*t1+3*D(1, 2)*t1*t1+D(0, 3)*t1**3
           +3*(D(1, 1)+D(0, 2)*t1)*t2)/w
    return [theta0, t1, t2, t3]


def pulled_back_third(function, x0, phase_jet):
    theta, t1, t2, t3 = phase_jet
    D = lambda i, j: mp.diff(function, (x0, theta), (i, j))
    return (D(3, 0)+3*D(2, 1)*t1+3*D(1, 2)*t1*t1+D(0, 3)*t1**3
            +3*(D(1, 1)+D(0, 2)*t1)*t2+D(0, 1)*t3)


def fixed_frequency_derivative(function, y, eta, N, radial, angular):
    """D_y^r D_eta^s [function(y,eta,N*y)/N], retaining every N power."""
    if type(N) is not int or N < 1 or min(radial, angular) < 0:
        raise ValueError('Positive fixed integer frequency and derivative orders required')
    phase = N*y
    return mp.fsum(comb(radial, fast)*mp.mpf(N)**(fast-1)
                   *mp.diff(function, (y, eta, phase), (radial-fast, angular, fast))
                   for fast in range(radial+1))
