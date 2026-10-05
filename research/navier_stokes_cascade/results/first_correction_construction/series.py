"""Finite first-correction Taylor construction for manufactured polynomial data.

This is a diagnostic of the exact recurrence, not an evaluation of the actual
large-parameter joined profile. Angular rows are ordinary Taylor coefficients
at eta=0; every radial step consumes angular derivatives explicitly.
"""
from pathlib import Path
import sys
import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent/'nonlinear_axis_pilot'))
from inner import add, cut, derivative, inv, mul, scale
from residual import Poly, ManufacturedProfile
sys.path.insert(0, str(Path(__file__).resolve().parent))


def build(base, h, nu, degree, angular_padding=20):
    if type(degree) is not int or not 1 <= degree <= 12 or angular_padding < 8:
        raise ValueError('Positive radial degree and sufficient angular padding required')
    h, nu = mp.mpf(str(h)), mp.mpf(str(nu))
    if not 0 < h < mp.mpf('.25') or nu <= 0:
        raise ValueError('Invalid diagnostic parameters')
    A, D, lam = mp.mpf('.5')+h, mp.mpf('.5')-h, 2*h
    budget = degree+angular_padding+2
    eta, d, L = [mp.mpf(0), mp.mpf(1)], [mp.mpf(1), 0, -1], [mp.mpf(1), 0, -2*h]
    Li = inv(L, budget)

    def row(poly, n):
        return cut([mp.mpf(str(poly.c.get((n, k), 0))) for k in range(budget+1)], budget)

    def get(rows, n):
        return rows[n] if 0 <= n < len(rows) else []

    def T(k, rows, n, m):
        f = get(rows, n)
        return mul(Li, add(scale(f, n-k, m), scale(mul(eta, derivative(f), m), D, m), m=m), m)

    def Z(k, rows, n, m):
        f = get(rows, n)
        return mul(Li, add(scale(mul(eta, f, m), 2*(k-n), m), mul(d, derivative(f), m), m=m), m)

    def ZZ(k, rows, n, m):
        intermediate = [[] for _ in range(n)]+[Z(k, rows, n, m+1)]
        return Z(k-D, intermediate, n, m)

    def conv(left, right, n, m, weight=lambda j: 1):
        return add(*(scale(mul(get(left, i), get(right, n-i), m), weight(n-i), m)
                     for i in range(n+1)), m=m)

    def mixed(left, right, k, n, m):
        return add(*(mul(get(left, i), Z(k, right, n-i, m), m) for i in range(n+1)), m=m)

    def V(rows, shift):
        values = [[]]
        for n, u in enumerate(rows):
            m = len(u)-2
            values.append(mul(Li, add(scale(mul(eta, u, m), 2-2*(D+shift)/(n+1), m),
                                      scale(mul(d, derivative(u), m), -mp.mpf(1)/(n+1), m), m=m), m))
        return values

    F0, U0 = [[row(poly, n) for n in range(degree+2)] for poly in (base.F, base.U)]
    V0 = V(U0, 0)

    def omega(n, m):
        return add(T(0, V0, n, m), conv(V0, V0, n+1, m, lambda j: j-mp.mpf('.5')),
                   mixed(U0, V0, 0, n, m),
                   scale(get(V0, n+1), -2*nu*n*(n+1), m), m=m)

    F1, U1, Pi1 = [[cut([], budget)] for _ in range(3)]
    retained = []
    for n in range(degree):
        m = budget-n-2
        V1 = V(U1, lam)
        theta = add(T(-A-mp.mpf('.5')+lam, F1, n, m),
                    conv(V0, F1, n+1, m, lambda j: j+1),
                    conv(V1, F0, n+1, m, lambda j: j+1),
                    mixed(U0, F1, -A-mp.mpf('.5')+lam, n, m),
                    mixed(U1, F0, -A-mp.mpf('.5'), n, m),
                    scale(ZZ(-A-mp.mpf('.5'), F0, n, m), -nu, m), m=m)
        axial = add(T(-A+lam, U1, n, m),
                    conv(V0, U1, n+1, m, lambda j: j),
                    conv(V1, U0, n+1, m, lambda j: j),
                    mixed(U0, U1, -A+lam, n, m), mixed(U1, U0, -A, n, m),
                    Z(-2*A+lam, Pi1, n, m), scale(ZZ(-A, U0, n, m), -nu, m), m=m)
        pressure = add(scale(conv(F0, F1, n, m), 2, m), scale(omega(n+1, m), -mp.mpf('.5'), m), m=m)
        F1.append(scale(theta, mp.mpf(1)/(2*nu*(n+1)*(n+2)), m))
        U1.append(scale(axial, mp.mpf(1)/(2*nu*(n+1)**2), m))
        Pi1.append(scale(pressure, mp.mpf(1)/(n+1), m))
        retained.append(m)

    def polynomial(rows):
        return Poly({(n, k): float(value) for n, values in enumerate(rows)
                     for k, value in enumerate(values) if value})

    correction = ManufacturedProfile(dict(F=[], U=[], axis_pressure=[]))
    correction.F, correction.U, correction.Pi = map(polynomial, (F1, U1, Pi1))
    correction.average_U = Poly({(n, k): value/(n+1) for (n, k), value in correction.U.c.items()})
    return correction, dict(manufactured=True, radial_degree=degree, angular_budget=budget,
                            retained_angular_orders=retained,
                            actual_large_parameter_profile_evaluated=False)
