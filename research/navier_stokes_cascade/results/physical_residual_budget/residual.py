"""Similarity-coordinate residual coefficients, independently of physical jets.

The imported arithmetic is reused; the residual assembly uses different
coordinates from the frozen Cartesian checker. All profiles used by this
module's numerical reviews are manufactured, not the actual huge-N profile.
"""
from fractions import Fraction
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
sys.path.insert(0, str(HERE.parent / 'paper_profile_audit'))
from profiles import (Jet, Poly, ManufacturedProfile, compose, coordinates,
                      field, physical_point, require_h)


def parameters(h):
    require_h(h)
    return .5+h, .5-h


def time_operator(k, f, X, eta, h):
    _, D = parameters(h)
    return (-k*f.v + X*f.g[0] + D*eta*f.g[1])/(1-2*h*eta*eta)


def axial_operator(k, f, X, eta, h):
    require_h(h)
    return (2*k*eta*f.v-2*eta*X*f.g[0]+(1-eta*eta)*f.g[1])/(1-2*h*eta*eta)


def axial_second(k, f, X, eta, h, shift_power=True):
    """Z_(k-D) Z_k, including derivatives of the coordinate coefficients."""
    _, D = parameters(h)
    L, d = 1-2*h*eta*eta, 1-eta*eta
    a, b, c = -2*eta*X/L, d/L, 2*k*eta/L
    ax = -2*eta/L
    ae = -2*X*(1+2*h*eta*eta)/L**2
    be = -2*(1-2*h)*eta/L**2
    ce = 2*k*(1+2*h*eta*eta)/L**2
    outer_c = 2*(k-D if shift_power else k)*eta/L
    return (a*a*f.H[0, 0]+2*a*b*f.H[0, 1]+b*b*f.H[1, 1]
            +(a*ax+b*ae+(c+outer_c)*a)*f.g[0]
            +(b*be+(c+outer_c)*b)*f.g[1]
            +(b*ce+outer_c*c)*f.v)


def radial_profile(profile, X, eta, h, lam=0., omit_lambda=False):
    """V=r*u_r/q^lam; V/X is the regular Cartesian coefficient v0."""
    _, D = parameters(h)
    average = profile.average_U(X, eta)
    average_eta = profile.average_U.derivative(1)(X, eta)
    shifted_D = D if omit_lambda else D+lam
    return X*(2*eta*profile.U(X, eta)-2*shifted_D*eta*average
              -(1-eta*eta)*average_eta)/(1-2*h*eta*eta)


def profile_jets(profile, X, eta, h, lam):
    x, e = Jet.variable(X, 0), Jet.variable(eta, 1)
    return dict(F=Jet.cast(profile.F(x, e)), U=Jet.cast(profile.U(x, e)),
                Pi=Jet.cast(profile.Pi(x, e)), V=radial_profile(profile, x, e, h, lam))


def series_coefficients(profiles, X, eta, h, nu=1.):
    """Exact coefficient assembly for any finite q^(2hn) profile sum.

    radial: r*q^(2A)*R_r; theta: q^(A+1)*R_theta/sqrt(2X);
    axial: q^(A+1)*R_z. Arrays retain all quadratic and axial-viscous terms.
    """
    A, D = parameters(h)
    if not profiles or not math.isfinite(nu) or nu <= 0 or X <= 0 or abs(eta) >= 1:
        raise ValueError('Positive viscosity, radius and an interior angle required')
    b, c = -A-.5, -A
    jets = [profile_jets(p, X, eta, h, 2*h*n) for n, p in enumerate(profiles)]
    count = len(jets)

    def products(n):
        return [(i, j) for i in range(count) for j in range(count) if i+j == n]

    omega = []
    for n in range(2*count-1):
        value = 0.
        if n < count:
            V = jets[n]['V']
            value = time_operator(2*h*n, V, X, eta, h)-2*nu*X*V.H[0, 0]
        for i, j in products(n):
            V, U, W = jets[i]['V'], jets[i]['U'], jets[j]['V']
            value += V.v*(W.g[0]-W.v/(2*X))+U.v*axial_operator(2*h*j, W, X, eta, h)
        omega.append(value)

    tangential = [[], []]
    for n in range(max(2*count-1, count+1)):
        theta = axial = 0.
        if n < count:
            F, U, Pi = [jets[n][key] for key in ('F', 'U', 'Pi')]
            theta = time_operator(b+2*h*n, F, X, eta, h)-2*nu*(X*F.H[0, 0]+2*F.g[0])
            axial = (time_operator(c+2*h*n, U, X, eta, h)
                     -2*nu*(X*U.H[0, 0]+U.g[0])
                     +axial_operator(-2*A+2*h*n, Pi, X, eta, h))
        for i, j in products(n):
            V, U, F, W = jets[i]['V'], jets[i]['U'], jets[j]['F'], jets[j]['U']
            theta += V.v*(F.g[0]+F.v/X)+U.v*axial_operator(b+2*h*j, F, X, eta, h)
            axial += V.v*W.g[0]+U.v*axial_operator(c+2*h*j, W, X, eta, h)
        if 0 <= n-1 < count:
            theta -= nu*axial_second(b+2*h*(n-1), jets[n-1]['F'], X, eta, h)
            axial -= nu*axial_second(c+2*h*(n-1), jets[n-1]['U'], X, eta, h)
        tangential[0].append(theta)
        tangential[1].append(axial)

    radial = []
    for n in range(max(2*count, count+2)):
        value = 2*X*jets[n]['Pi'].g[0] if n < count else 0.
        value -= 2*X*sum(jets[i]['F'].v*jets[j]['F'].v for i, j in products(n))
        if 0 <= n-1 < len(omega):
            value += omega[n-1]
        if 0 <= n-2 < count:
            value -= nu*axial_second(2*h*(n-2), jets[n-2]['V'], X, eta, h)
        radial.append(value)
    return dict(radial=radial, theta=tangential[0], axial=tangential[1], omega=omega)


def evaluate_series(coefficients, q, X, h):
    A, _ = parameters(h)
    if not math.isfinite(q) or q <= 0 or X <= 0:
        raise ValueError('Positive finite q and X required')
    eps, R = q**(2*h), math.sqrt(2*X)
    def polynomial(row):
        return sum(v*eps**n for n, v in enumerate(row))
    return np.array([q**(-2*A-.5)*polynomial(coefficients['radial'])/R,
                     R*q**(-A-1)*polynomial(coefficients['theta']),
                     q**(-A-1)*polynomial(coefficients['axial'])])


def physical_series(profiles, point, h, omit_lambda=False):
    """Direct Cartesian ansatz, with physical-coordinate derivative jets."""
    A, _ = parameters(h)
    x, y, _, _, q, X, eta = coordinates(point, h)
    velocity, pressure = [Jet(0.) for _ in range(3)], Jet(0.)
    for n, profile in enumerate(profiles):
        lam = 2*h*n
        average = profile.average_U(X, eta)
        average_eta = profile.average_U.derivative(1)(X, eta)
        Dshift = .5-h+(0 if omit_lambda else lam)
        v0 = (2*eta*profile.U(X, eta)-2*Dshift*eta*average
              -(1-eta*eta)*average_eta)/(1-2*h*eta*eta)
        radial = q**(lam-1)*v0/2
        swirl = q**(-A-.5+lam)*profile.F(X, eta)
        row = [radial*x-swirl*y, radial*y+swirl*x, q**(-A+lam)*profile.U(X, eta)]
        velocity = [a+b for a, b in zip(velocity, row)]
        pressure += q**(-2*A+lam)*profile.Pi(X, eta)
    return velocity, pressure


def momentum(velocity, pressure, nu=1.):
    values = np.array([v.v for v in velocity])
    terms = dict(time=np.array([v.g[3] for v in velocity]),
                 advection=np.array([np.dot(values, v.g[:3]) for v in velocity]),
                 viscosity=-nu*np.array([np.trace(v.H[:3, :3]) for v in velocity]),
                 pressure=pressure.g[:3])
    scales = np.maximum(1., np.max(np.abs(list(terms.values())), axis=0))
    return sum(terms.values()), scales


def cylindrical(vector, theta):
    ct, st = math.cos(theta), math.sin(theta)
    return np.array([ct*vector[0]+st*vector[1], -st*vector[0]+ct*vector[1], vector[2]])


def stress_pair(profile, X, eta, h, nu=1.):
    A, D = parameters(h)
    d, L, R = 1-eta*eta, 1-2*h*eta*eta, (2*X)**.5
    F, U, Pi = profile.F(X, eta), profile.U(X, eta), profile.Pi(X, eta)
    W = 1-2*D*eta*profile.average_U(X, eta)-d*profile.average_U.derivative(1)(X, eta)
    I, J, M, S, H = [getattr(profile, key)(X, eta) for key in ('I', 'J', 'M', 'S', 'H')]
    Ie, Je, Me, Se, Pie = [getattr(profile, key).derivative(1)(X, eta)
                          for key in ('I', 'J', 'M', 'S', 'Pi')]
    Qs = -W+((1-h)*I-D*eta*Ie-d*Je+2*(h-D)*eta*J)/(X*H)
    Ns = -W*U+(D*(M-eta*Me)+4*h*eta*S-d*Se)/X+4*A*eta*Pi-d*Pie
    return (F*X*Qs/L+2*nu*X*profile.F.derivative(0)(X, eta),
            X*Ns/(L*R)+nu*R*profile.U.derivative(0)(X, eta))


def symmetric_cross_stress_divergence(profile, point, h):
    """Divergence of an illustrative symmetric tensor with ONLY r-theta,r-z.

    This is not the manuscript's two-component cylindrical stress operator
    and is not the full positive wave covariance tensor.
    """
    A, _ = parameters(h)
    x, y, _, _, q, X, eta = coordinates(point, h)
    radius = (x*x+y*y)**.5
    er, et, ez = [x/radius, y/radius, Jet(0.)], [-y/radius, x/radius, Jet(0.)], [Jet(0.), Jet(0.), Jet(1.)]
    theta, axial = [q**(-A-.5)*v for v in stress_pair(profile, X, eta, h)]
    tensor = [[theta*(er[i]*et[j]+et[i]*er[j])+axial*(er[i]*ez[j]+ez[i]*er[j])
               for j in range(3)] for i in range(3)]
    return np.array([sum(tensor[i][j].g[j] for j in range(3)) for i in range(3)])


def cutoff_velocity(profile, point, h, scale, preserve_divergence):
    """Diagnostic chi(sigma)=exp(-sigma), applied to velocity or streamfunction.

    The diagnostic cutoff is not compactly supported. Its purpose is to test
    the exact product rule, not to supply the final summation cutoffs.
    """
    velocity, _ = physical_series([profile], point, h)
    x, y, _, _, q, X, eta = coordinates(point, h)
    arg = scale*q
    chi = compose(arg, math.exp(-arg.v), -math.exp(-arg.v), math.exp(-arg.v))
    result = [chi*v for v in velocity]
    if preserve_divergence:
        # -(q_z/r)*scale*chi'(scale*q)*q^D*M in the radial direction.
        # chi'=-chi; M/X=average_U, so the Cartesian extension is regular.
        factor = eta*arg*chi*profile.average_U(X, eta)/((1-2*h*eta*eta)*q)
        result[0] += factor*x
        result[1] += factor*y
    return result


def derivative_loss(h, transverse=0, axial=0, time=0):
    h = Fraction(h)
    if not 0 < h < Fraction(1, 2) or any(type(n) is not int or n < 0 for n in (transverse, axial, time)):
        raise ValueError('Invalid exponent or derivative counts')
    return Fraction(transverse, 2)+(Fraction(1, 2)-h)*axial+time


def minimum_order_for_decay(h, transverse=0, axial=0, time=0):
    """Necessary order for a positive exponent in the stated sufficient budget.

    Assumes every retained normalized coefficient is actually cancelled.
    Does not assert solvability, uniform constants, or a convergent series.
    """
    h = Fraction(h)
    loss = derivative_loss(h, transverse, axial, time)
    threshold = (Fraction(3, 2)+2*h+loss)/(2*h)
    return threshold.numerator//threshold.denominator
