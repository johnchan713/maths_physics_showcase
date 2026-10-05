"""Derived six-variable system for the first inner momentum correction.

Numerical evaluations use manufactured polynomial profiles. The analytic
construction applies to the inherited unmodulated, protected inner profile.
"""
import math
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent/'physical_residual_budget'))
from residual import Jet, axial_operator, axial_second, parameters, profile_jets
sys.path.insert(0, str(Path(__file__).resolve().parent))

DIAGONAL = np.array([0., 0., 2., 0., 3., 1.])


def radial_regular_jet(profile, X, eta, h):
    x, e = Jet.variable(X, 0), Jet.variable(eta, 1)
    _, D = parameters(h)
    return (2*e*profile.U(x, e)-2*D*e*profile.average_U(x, e)
            -(1-e*e)*profile.average_U.derivative(1)(x, e))/(1-2*h*e*e)


def omega_over_X(profile, X, eta, h, nu=1.):
    """Regular at X=0; no subtraction of nearly singular Omega/X values."""
    _, D = parameters(h)
    v = radial_regular_jet(profile, X, eta, h)
    U = profile.U(X, eta)
    d, L = 1-eta*eta, 1-2*h*eta*eta
    return ((v.v+X*v.g[0]+D*eta*v.g[1])/L
            +v.v*(v.v/2+X*v.g[0])
            +U*(-2*eta*(v.v+X*v.g[0])+d*v.g[1])/L
            -2*nu*(2*v.g[0]+X*v.H[0, 0]))


def matrices(profile, xi, eta, h, nu=1.):
    A, D = parameters(h)
    if not math.isfinite(xi) or xi < 0 or abs(eta) > 1 or not math.isfinite(nu) or nu <= 0:
        raise ValueError('Invalid inner system coordinates or viscosity')
    X, lam = xi*xi, 2*h
    b, c, kp = -A-.5, -A, -2*A+lam
    d, L = 1-eta*eta, 1-2*h*eta*eta
    jets = profile_jets(profile, X, eta, h, 0.)
    F, U = jets['F'], jets['U']
    v = radial_regular_jet(profile, X, eta, h).v
    Hc = D*eta+d*U.v
    JF, JU = F.v+X*F.g[0], X*U.g[0]
    ZF, ZU = axial_operator(b, F, X, eta, h), axial_operator(c, U, X, eta, h)
    omega = omega_over_X(profile, X, eta, h, nu)
    CF = axial_second(b, F, X, eta, h)
    CU = axial_second(c, U, X, eta, h)
    a0, a1 = np.zeros((6, 6)), np.zeros((6, 6))
    a0[0, 4] = a0[1, 5] = 1.
    a0[2, 5] = -1.
    a0[3, 0] = 4*xi*F.v
    transport = xi*(.5-eta*U.v)/L+xi*v/2
    a0[4, 0] = 2*((b+lam)*(2*eta*U.v-1)/L+v)/nu
    a0[4, 1] = 2*(2*(A-lam)*eta*JF/L+ZF)/nu
    a0[4, 2] = -4*(D+lam)*eta*JF/(nu*L)
    a0[4, 4] = 2*transport/nu
    a0[5, 0] = -8*eta*X*F.v/(nu*L)
    a0[5, 1] = 2*((c+lam)*(2*eta*U.v-1)/L+ZU+2*(A-lam)*eta*JU/L)/nu
    a0[5, 2] = -4*(D+lam)*eta*JU/(nu*L)
    a0[5, 3] = 4*kp*eta/(nu*L)
    a0[5, 5] = 2*transport/nu
    a1[4, 0] = 2*Hc/(nu*L)
    a1[4, 1] = a1[4, 2] = -2*d*JF/(nu*L)
    a1[5, 1] = 2*(Hc-d*JU)/(nu*L)
    a1[5, 2] = -2*d*JU/(nu*L)
    a1[5, 3] = 2*d/(nu*L)
    source = np.array([0., 0., 0., -xi*omega, -2*CF,
                       2*eta*X*omega/(nu*L)-2*CU])
    return a0, a1, source


def correction_vector(profile, xi, eta):
    """W, its xi derivative and eta derivative for an arbitrary regular test."""
    X = xi*xi
    x, e = Jet.variable(X, 0), Jet.variable(eta, 1)
    F, U, Pi, average = [Jet.cast(f(x, e)) for f in (profile.F, profile.U, profile.Pi, profile.average_U)]
    K = average-U
    w = np.array([F.v, U.v, K.v, Pi.v, 2*xi*F.g[0], 2*xi*U.g[0]])
    wx = np.array([2*xi*F.g[0], 2*xi*U.g[0], 2*xi*K.g[0], 2*xi*Pi.g[0],
                   2*F.g[0]+4*X*F.H[0, 0], 2*U.g[0]+4*X*U.H[0, 0]])
    we = np.array([F.g[1], U.g[1], K.g[1], Pi.g[1],
                   2*xi*F.H[0, 1], 2*xi*U.H[0, 1]])
    return w, wx, we


def axis_slopes(profile, eta, h, nu=1.):
    A, _ = parameters(h)
    jets = profile_jets(profile, 0., eta, h, 0.)
    return dict(F1_X=-axial_second(-A-.5, jets['F'], 0., eta, h)/4,
                U1_X=-axial_second(-A, jets['U'], 0., eta, h)/2,
                Pi1_X=-omega_over_X(profile, 0., eta, h, nu)/2)
