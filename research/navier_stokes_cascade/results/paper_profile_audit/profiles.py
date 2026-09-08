"""Explicit paper building blocks, not the nonlinear matched blowup construction.

Source: the paper identified by protocol.json. Derivatives below retain all
Cartesian terms, including axial diffusion and cylindrical curvature.
"""
from fractions import Fraction
import math
import warnings

import numpy as np
from scipy.integrate import IntegrationWarning, quad
from scipy.optimize import brentq
from scipy.special import gamma, hyperu, j1, poch


class Jet:
    """Value, gradient and Hessian in Cartesian (x, y, z, t)."""
    def __init__(self, value, gradient=None, hessian=None):
        self.v = float(value)
        self.g = np.zeros(4) if gradient is None else np.asarray(gradient, dtype=float)
        self.H = np.zeros((4, 4)) if hessian is None else np.asarray(hessian, dtype=float)

    @staticmethod
    def variable(value, axis):
        gradient = np.zeros(4)
        gradient[axis] = 1.
        return Jet(value, gradient)

    @staticmethod
    def cast(value):
        return value if isinstance(value, Jet) else Jet(value)

    def __add__(self, other):
        other = Jet.cast(other)
        return Jet(self.v + other.v, self.g + other.g, self.H + other.H)

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.v, -self.g, -self.H)

    def __sub__(self, other):
        return self + (-Jet.cast(other))

    def __rsub__(self, other):
        return Jet.cast(other) - self

    def __mul__(self, other):
        other = Jet.cast(other)
        return Jet(self.v * other.v, self.g * other.v + other.g * self.v,
                   self.H * other.v + other.H * self.v +
                   np.outer(self.g, other.g) + np.outer(other.g, self.g))

    __rmul__ = __mul__

    def __pow__(self, power):
        if power == 0:
            return Jet(1.)
        if power == 1:
            return self
        if not math.isfinite(power) or (self.v <= 0 and (power < 0 or not float(power).is_integer())):
            raise ValueError("Invalid real power in derivative jet")
        return compose(self, self.v**power, power * self.v**(power - 1),
                       power * (power - 1) * self.v**(power - 2))

    def __truediv__(self, other):
        return self * Jet.cast(other)**-1

    def __rtruediv__(self, other):
        return Jet.cast(other) / self


def compose(argument, value, first, second):
    return Jet(value, first * argument.g,
               first * argument.H + second * np.outer(argument.g, argument.g))


class Poly:
    """Small bivariate polynomial in (X, eta), with exact monomial calculus."""
    def __init__(self, coefficients):
        self.c = {tuple(k): float(v) for k, v in coefficients.items() if v != 0}

    @staticmethod
    def cast(value):
        return value if isinstance(value, Poly) else Poly({(0, 0): value})

    def __call__(self, x, eta):
        return sum(value * x**i * eta**j for (i, j), value in self.c.items())

    def __add__(self, other):
        result = self.c.copy()
        for key, value in Poly.cast(other).c.items():
            result[key] = result.get(key, 0.) + value
        return Poly(result)

    __radd__ = __add__

    def __mul__(self, other):
        result = {}
        for (i, j), a in self.c.items():
            for (k, l), b in Poly.cast(other).c.items():
                key = i + k, j + l
                result[key] = result.get(key, 0.) + a * b
        return Poly(result)

    __rmul__ = __mul__

    def derivative(self, axis):
        result = {}
        for key, value in self.c.items():
            if key[axis]:
                reduced = list(key)
                reduced[axis] -= 1
                result[tuple(reduced)] = value * key[axis]
        return Poly(result)

    def integral(self):
        return Poly({(i + 1, j): v / (i + 1) for (i, j), v in self.c.items()})

    def average(self):
        return Poly({(i, j): v / (i + 1) for (i, j), v in self.c.items()})


class ManufacturedProfile:
    """Regular polynomial fixture; it does not assert zero leading stress."""
    def __init__(self, settings):
        def polynomial(name):
            return Poly({(i, j): v for i, j, v in settings[name]})
        self.F, self.U = polynomial("F"), polynomial("U")
        x = Poly({(1, 0): 1.})
        self.Pi = polynomial("axis_pressure") + (self.F * self.F).integral()
        self.H = 2 * x * self.F
        self.M = self.U.integral()
        self.I = self.H.integral()
        self.J = (self.U * self.H).integral()
        self.S = (self.U * self.U + (-1) * x * self.F * self.F).integral()
        self.average_U = self.U.average()


def require_h(h):
    if not math.isfinite(h) or not 0 < h < .5:
        raise ValueError("Coordinate and heat formulas require 0 < h < 1/2")


def q_value(z, time, h):
    require_h(h)
    if not all(math.isfinite(v) for v in (z, time)) or not 0 <= time < 1:
        raise ValueError("Need finite coordinates and 0 <= t < 1")
    tau = 1 - time
    if z == 0:
        return tau
    low = max(tau, abs(z)**(1 / (.5 - h)))
    high = 2 * (tau + low)

    def equation(q):
        return q - z*z*q**(2*h) - tau

    while equation(high) <= 0:
        high *= 2
        if not math.isfinite(high):
            raise ValueError("Similarity root overflow")
    if equation(low) > 0:
        # Only a lower-bracket rounding excursion is allowed.
        if equation(low) > 32 * np.finfo(float).eps * low:
            raise ValueError("Invalid similarity root bracket")
        low *= 1 - 64 * np.finfo(float).eps
    return brentq(equation, low, high, xtol=np.nextafter(0., 1.),
                  rtol=8 * np.finfo(float).eps, maxiter=128)


def coordinates(point, h, differentiated=True):
    if len(point) != 4 or not np.isfinite(point).all():
        raise ValueError("Need a finite Cartesian space-time point")
    q = q_value(point[2], point[3], h)
    D = .5 - h
    if differentiated:
        # Implicit differentiation of G=q-z^2 q^(2h)-(1-t)=0.
        z = point[2]
        L = 1 - 2*h*z*z*q**(2*h - 1)
        gradient = np.array([0., 0., 2*z*q**(2*h)/L, -1/L])
        Gqq = -2*h*(2*h - 1)*z*z*q**(2*h - 2)
        Gqi = np.array([0., 0., -4*h*z*q**(2*h - 1), 0.])
        Gij = np.zeros((4, 4))
        Gij[2, 2] = -2*q**(2*h)
        hessian = -(Gqq*np.outer(gradient, gradient) + np.outer(Gqi, gradient) +
                    np.outer(gradient, Gqi) + Gij) / L
        q = Jet(q, gradient, hessian)
        x, y, z, time = [Jet.variable(v, i) for i, v in enumerate(point)]
    else:
        x, y, z, time = point
    eta = z * q**-D
    X = (x*x + y*y) / (2*q)
    return x, y, z, time, q, X, eta


def physical_point(q, eta, X, theta, h):
    require_h(h)
    if not all(math.isfinite(v) for v in (q, eta, X, theta)) or not (q > 0 and X >= 0 and abs(eta) < 1):
        raise ValueError("Invalid similarity point")
    radius = math.sqrt(2*q*X)
    return np.array([radius*math.cos(theta), radius*math.sin(theta),
                     q**(.5 - h)*eta, 1 - q*(1 - eta*eta)])


def field(point, h, profile, differentiated=True):
    x, y, z, time, q, X, eta = coordinates(point, h, differentiated)
    A, D = .5 + h, .5 - h
    d, L = 1 - eta*eta, 1 - 2*h*eta*eta
    F, U = profile.F(X, eta), profile.U(X, eta)
    average = profile.average_U(X, eta)
    average_eta = profile.average_U.derivative(1)(X, eta)
    v0 = (2*eta*U - 2*D*eta*average - d*average_eta) / L
    radial_factor, angular_factor = v0/(2*q), q**(-A - .5)*F
    velocity = [radial_factor*x - angular_factor*y,
                radial_factor*y + angular_factor*x, q**-A*U]
    pressure = q**(-2*A)*profile.Pi(X, eta)
    return velocity, pressure, (q, X, eta, v0)


def full_residual(velocity, pressure):
    values = np.array([v.v for v in velocity])
    terms = {
        "time": np.array([v.g[3] for v in velocity]),
        "advection": np.array([np.dot(values, v.g[:3]) for v in velocity]),
        "viscosity": -np.array([np.trace(v.H[:3, :3]) for v in velocity]),
        "pressure": pressure.g[:3].copy(),
    }
    residual = sum(terms.values())
    scale = max(1., *(float(np.max(np.abs(v))) for v in terms.values()))
    return residual, scale, terms


def scalar_radial_derivatives(jet, radial):
    return float(np.dot(radial, jet.g[:2])), float(radial @ jet.H[:2, :2] @ radial)


def cylindrical_and_stress(point, h, profile):
    """Independent cylindrical assembly, then integrated-stress identity."""
    x, y, z, time, q, X, eta = coordinates(point, h)
    radius = (x*x + y*y)**.5
    if radius.v == 0:
        raise ValueError("Use the Cartesian regular extension at the axis")
    er = np.array(point[:2]) / radius.v
    et = np.array([-er[1], er[0]])
    A, D = .5 + h, .5 - h
    d, L = 1 - eta*eta, 1 - 2*h*eta*eta
    average = profile.average_U(X, eta)
    W = 1 - 2*D*eta*average - d*profile.average_U.derivative(1)(X, eta)
    U, F, Pi = profile.U(X, eta), profile.F(X, eta), profile.Pi(X, eta)
    v0 = (2*eta*U + W - 1) / L
    ur = radius*v0/(2*q)
    ut, uz = q**-A*(2*X)**.5*F, q**-A*U
    pressure = q**(-2*A)*Pi
    components = [ur, ut, uz]
    derivatives = [scalar_radial_derivatives(v, er) for v in components]
    residual = []
    for j, (v, (vr, vrr)) in enumerate(zip(components, derivatives)):
        value = v.g[3] + ur.v*vr + uz.v*v.g[2] - vrr - vr/radius.v - v.H[2, 2]
        if j == 0:
            value += -ut.v**2/radius.v + ur.v/radius.v**2 + np.dot(er, pressure.g[:2])
        elif j == 1:
            value += ur.v*ut.v/radius.v + ut.v/radius.v**2
        else:
            value += pressure.g[2]
        residual.append(value)
    cartesian = np.array([er[0]*residual[0] + et[0]*residual[1],
                          er[1]*residual[0] + et[1]*residual[1], residual[2]])
    H, M, I, J, S = [p(X, eta) for p in (profile.H, profile.M, profile.I, profile.J, profile.S)]
    Me, Ie, Je, Se, Pie = [p.derivative(1)(X, eta) for p in
                          (profile.M, profile.I, profile.J, profile.S, profile.Pi)]
    Qs = -W + ((1 - h)*I - D*eta*Ie - d*Je + 2*(h - D)*eta*J)/(X*H)
    Ns = -W*U + (D*(M - eta*Me) + 4*h*eta*S - d*Se)/X + 4*A*eta*Pi - d*Pie
    # Algebraically simplified (4.11), with F=E/sqrt(2X).
    Ttheta = q**(-A - .5)*(F*X*Qs/L + 2*X*profile.F.derivative(0)(X, eta))
    Tz = q**(-A - .5)*(X*Ns/(L*(2*X)**.5) + (2*X)**.5*profile.U.derivative(0)(X, eta))
    stress_rhs = -np.array([np.dot(er, Ttheta.g[:2]) + 2*Ttheta.v/radius.v,
                            np.dot(er, Tz.g[:2]) + Tz.v/radius.v])
    axial_diffusion = np.array([ut.H[2, 2], uz.H[2, 2]])
    return {"cartesian": cartesian, "cylindrical": np.array(residual),
            "leading_tangential": np.array(residual[1:]) + axial_diffusion,
            "stress_rhs": stress_rhs, "axial_diffusion": axial_diffusion}


def finite_difference_residual(point, h, profile, relative_step):
    """Fourth-order physical-coordinate stencils; no derivative jets."""
    if not math.isfinite(relative_step) or relative_step <= 0:
        raise ValueError("Finite-difference step must be positive")
    q = q_value(point[2], point[3], h)
    widths = relative_step*np.array([q**.5, q**.5, q**(.5-h), q])

    def values(location):
        u, p, _ = field(location, h, profile, differentiated=False)
        return np.array([*u, p])

    center = values(point)
    first = np.empty((4, 4))
    laplacian = np.zeros(3)
    for axis, step in enumerate(widths):
        samples = []
        for multiple in (-2, -1, 1, 2):
            offset = np.array(point, dtype=float)
            offset[axis] += multiple*step
            samples.append(values(offset))
        a, b, c, d = samples
        first[:, axis] = (a - 8*b + 8*c - d)/(12*step)
        if axis < 3:
            laplacian += (-a[:3] + 16*b[:3] - 30*center[:3] + 16*c[:3] - d[:3])/(12*step*step)
    return first[:3, 3] + first[:3, :3] @ center[:3] - laplacian + first[3, :3]


def heat_derivatives(h, Z, maximum_order=3):
    require_h(h)
    if not math.isfinite(Z) or Z < 0 or not isinstance(maximum_order, int) or not 0 <= maximum_order <= 6:
        raise ValueError("Invalid heat-profile argument or derivative order")
    result, estimates = [], []
    for order in range(maximum_order + 1):
        coefficient = (-1)**order*poch(h, order)/gamma(1+h)
        with warnings.catch_warnings():
            warnings.simplefilter("error", IntegrationWarning)
            value, error = quad(lambda v: math.exp(-v)*v**(h+order)*(1+Z*v)**(-h-order),
                                0., np.inf, epsabs=2e-13, epsrel=2e-13, limit=200)
        result.append(float(coefficient*value))
        estimates.append(float(abs(coefficient)*error))
    return result, estimates


def heat_reference(h, Z, order=0):
    """Independent special-function implementation of the same integral."""
    require_h(h)
    if not math.isfinite(Z) or Z < 0 or not isinstance(order, int) or not 0 <= order <= 6:
        raise ValueError("Invalid heat reference argument")
    coefficient = (-1)**order*poch(h, order)*poch(1+h, order)
    if Z == 0:
        return float(coefficient)
    a = 1 + h + order
    return float(coefficient*Z**-a*hyperu(a, 2., 1/Z))


def heat_exterior(point, h):
    """Unit-amplitude exterior at viscosity one; valid only at r>0, t<=1."""
    require_h(h)
    if len(point) != 4 or not np.isfinite(point).all() or not 0 <= point[3] <= 1 or np.linalg.norm(point[:2]) == 0:
        raise ValueError("Heat exterior is not a regular-axis/global candidate")
    x, y, z, time = [Jet.variable(v, i) for i, v in enumerate(point)]
    s = (x*x + y*y)/2
    radius = (2*s)**.5
    Z = 2*(1-time)/s
    derivatives, _ = heat_derivatives(h, Z.v, 2)
    H = compose(Z, *derivatives)
    K = s**(-.5-h)*H
    velocity = [-K*y/radius, K*x/radius, Jet(0.)]
    # p=-integral_r^infinity K(rho,t)^2/rho drho; FTC gives this gradient.
    pressure = Jet(0., [K.v*K.v*point[0]/radius.v**2,
                       K.v*K.v*point[1]/radius.v**2, 0., 0.])
    residual, scale, terms = full_residual(velocity, pressure)
    er = np.array(point[:2])/radius.v
    Kr, Krr = scalar_radial_derivatives(K, er)
    missing_curvature = K.g[3] - Krr - Kr/radius.v
    return {"residual": residual, "scale": scale, "K": K.v, "Kr": Kr,
            "pressure_gradient": pressure.g[:3],
            "without_swirl_curvature": missing_curvature,
            "expected_curvature_error": -K.v/radius.v**2, "terms": terms}


def heat_pressure(radius, time, h):
    """Independent radial pressure integration; shares H's quadrature evaluator.

    The special-function H reference is not used inside this adaptive integral:
    numerical noise at some arguments swamped the pressure-difference check.
    """
    require_h(h)
    if not (math.isfinite(radius) and radius > 0 and math.isfinite(time) and 0 <= time <= 1):
        raise ValueError("Invalid exterior pressure point")
    s, A = radius*radius/2, .5+h
    Z = 2*(1-time)/s
    with warnings.catch_warnings():
        warnings.simplefilter("error", IntegrationWarning)
        value = quad(lambda y: y**(2*A-1)*heat_derivatives(h, Z*y, 0)[0][0]**2, 0., 1.,
                     epsabs=2e-12, epsrel=2e-12, limit=200)[0]
    return -.5*s**(-2*A)*value


def core_coefficients(terms=24):
    if not isinstance(terms, int) or not 4 <= terms <= 100:
        raise ValueError("Need 4 through 100 scalar-comparison coefficients")
    return [Fraction((-1)**n, 2**n*math.factorial(n)*math.factorial(n+1)) for n in range(terms)]


def core_comparison(z, terms=24):
    if not math.isfinite(z) or not 0 <= z <= 4.1:
        raise ValueError("Only the declared positive comparison interval is audited")
    value = math.fsum(float(c)*z**n for n, c in enumerate(core_coefficients(terms)))
    reference = 1. if z == 0 else 2*j1(math.sqrt(2*z))/math.sqrt(2*z)
    next_term = z**terms/(2**terms*math.factorial(terms)*math.factorial(terms+1))
    return value, float(reference), next_term


def normalized_error(a, b, floor=1.):
    a, b = np.asarray(a), np.asarray(b)
    if not math.isfinite(floor) or floor <= 0 or not a.size or not b.size or a.shape != b.shape:
        raise ValueError("Invalid comparison shape or scale")
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("Nonfinite comparison")
    return float(np.max(np.abs(a-b)) / max(floor, float(np.max(np.abs(a))), float(np.max(np.abs(b)))))
