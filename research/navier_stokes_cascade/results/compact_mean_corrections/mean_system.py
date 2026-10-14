"""Exact mean balances. Inputs are Taylor fields, not sampled PDE solutions.

Coordinates are (R,theta,Z,T). In normalized charts Dz=epsilon*d_Z and
tstar=-epsilon*d_T on auxiliary-independent data. Physical fixtures use
Dz=d_z and tstar=d_t. The continuum graph/torus terms are in README.md.
"""
from pathlib import Path
import sys

WAVE_FOLDER = Path(__file__).resolve().parent.parent/'physical_wave_interface'
if str(WAVE_FOLDER) not in sys.path:
    sys.path.append(str(WAVE_FOLDER))
from jets import Jet, arrays, bump, variables
import waves


def mean_correction(coords, epsilon=1):
    """Compact curl plus a zero-angular-moment azimuthal field."""
    r, _, z, time = coords
    radial = bump((r-1.4)/.48)
    psi = .27*radial*(1+.21*z.sin()+.16*time.cos()+.13*z*time)
    h = .19*radial*(1+.17*z.cos()+.23*time.sin())
    return [-epsilon*psi.diff(2), h.diff(0)/r**2, psi.diff(0)+psi/r]


def mean_increment(coords, epsilon=1):
    """A separate compact increment; used to test exact quadratic changes."""
    r, _, z, time = coords
    radial = bump((r-1.38)/.39)
    psi = .13*radial*(1+.31*z.cos()-.19*time.sin()+.11*z*time)
    h = -.11*radial*(1+.16*z.sin()-.14*time.cos())
    return [-epsilon*psi.diff(2), h.diff(0)/r**2, psi.diff(0)+psi/r]


def mean_balances(coords, base, correction, covariance, pressure=None, stress=None,
                  epsilon=1, normalized=False):
    """Return (gr,Etheta,Ez) of the complete conservative mean equations."""
    r = coords[0]
    b, vbase, gbase = base
    beta, swirl, gamma = correction
    w = covariance
    pressure = Jet(0) if pressure is None else pressure
    st, sz = (Jet(0), Jet(0)) if stress is None else stress
    dz = lambda f: epsilon*f.diff(2)
    dt = lambda f: -epsilon*f.diff(3) if normalized else f.diff(3)
    radial = lambda f, e: f.diff(0)+e*f/r
    lap = lambda f: f.diff(0).diff(0)+f.diff(0)/r+epsilon**2*f.diff(2).diff(2)
    gr = (-dt(beta)-radial(2*b*beta+beta*beta+w[0][0], 1)
          -dz(b*gamma+gbase*beta+beta*gamma+w[0][2])
          +(2*vbase*swirl+swirl*swirl+w[1][1])/r
          +epsilon*(lap(beta)-beta/r**2))
    et = (dt(swirl)+radial(b*swirl+beta*vbase+beta*swirl+w[0][1], 2)
          +dz(gbase*swirl+vbase*gamma+gamma*swirl+w[1][2])
          -epsilon*(lap(swirl)-swirl/r**2)-radial(st, 2))
    ez = (dt(gamma)+radial(b*gamma+beta*gbase+beta*gamma+w[0][2], 1)
          +dz(2*gbase*gamma+gamma*gamma+w[2][2]+pressure)
          -epsilon*lap(gamma)-radial(sz, 1))
    return gr, et, ez


def radial_increment_remainder(coords, base, old, increment, epsilon=1, normalized=False):
    """gr(new)-gr(old)-(2V/R)*Delta v, with every other term retained."""
    r = coords[0]
    b, _, g = base
    beta, v, gamma = old
    db, dv, dg = increment
    radial = lambda f: f.diff(0)+f/r
    dz = lambda f: epsilon*f.diff(2)
    dt = lambda f: -epsilon*f.diff(3) if normalized else f.diff(3)
    lap = lambda f: f.diff(0).diff(0)+f.diff(0)/r+epsilon**2*f.diff(2).diff(2)
    return (-dt(db)-radial(2*b*db+2*beta*db+db*db)
            -dz(b*dg+g*db+beta*dg+gamma*db+db*dg)
            +(2*v*dv+dv*dv)/r+epsilon*(lap(db)-db/r**2))


def scaled_data(data, epsilon):
    """Apply the normalized derivative multipliers to value/gradient/Hessian."""
    import numpy as np
    factors = np.array([1., 1., epsilon, -epsilon])
    return data[0], data[1]*factors, data[2]*factors[None, :, None]*factors[None, None, :]
