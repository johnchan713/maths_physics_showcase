"""Factorial-normalized, fixed-path stencil for N t' = B t."""
from fractions import Fraction
from math import factorial


def coefficients(normal, motion, matrix, gamma):
    """Positive D and E tables from fixed-v slow bounds already divided by i!."""
    order = len(normal)
    if not order or len(motion) != order or len(matrix) != order or gamma < 0:
        raise ValueError("inconsistent coefficient data")

    def product(left, right):
        out = {}
        for (a,b), u in left.items():
            for (i,j), v in right.items():
                if a+i < order:
                    key = (a+i,b+j)
                    out[key] = out.get(key,0) + u*v
        return out

    n = {(a,j): seq[a] for j,seq in enumerate((normal,motion)) for a in range(order)}
    m = {(a,0): val for a,val in enumerate(motion)}
    K = {(a,0): val for a,val in enumerate(matrix)}
    D = product(n,n)
    # -NK and nn^T K have the same separate norm majorant.
    nk, nm, dd = product(D,K), product(n,m), product(D,D)
    E = {key: 2*nk.get(key,0) + nm.get(key,0) + gamma*dd[key] for key in dd}
    return D,E


def path(d, total, inverse_floor, entrance, den, num):
    """Pure-slow normalized bounds and coefficient majorants; see (2)."""
    if not 0 <= d <= total or inverse_floor <= 0:
        raise ValueError("invalid orders or reciprocal floor")
    out = {(a, 0): entrance(a) for a in range(d + 1)}
    for b in range(1, total + 1):
        for a in range(min(d, total - b) + 1):
            rhs = sum(num(i, j) * out[a-i, b-1-j]
                      for i in range(a+1) for j in range(min(4, b-1)+1))
            lhs = sum((b-j) * den(i, j) * out[a-i, b-j]
                      for i in range(a+1) for j in range(min(2, b-1)+1)
                      if i or j)
            out[a, b] = Fraction(1, b) * inverse_floor * (rhs + lhs)
    return out


def record():
    # Nonconstant N=1+v², B=-(1+v²)². Exact scalar damping fixture at v=0.
    t = path(0, 6, 1, lambda a: 1, lambda i,j: int(j in (0,2)),
             lambda i,j: {0:1, 2:2, 4:1}.get(j,0))
    return {"fixture_raw_bounds": [int(t[0,j]*factorial(j)) for j in range(7)],
            "background_slow_order_offset": 99998,
            "physical_jet_enclosures_verified": False}
