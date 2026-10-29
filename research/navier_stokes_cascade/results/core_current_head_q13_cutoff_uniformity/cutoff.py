"""Exact rational diagnostics for the all-order implicit-coordinate cutoff proof."""
from fractions import Fraction as Q


def add(*ps):
    out = [Q(0)] * max(map(len, ps), default=0)
    for p in ps:
        for i, v in enumerate(p): out[i] += v
    while out and out[-1] == 0: out.pop()
    return tuple(out)


def scale(p, c):
    return add(tuple(c*v for v in p))


def mul(p, q):
    out = [Q(0)] * max(0, len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b in enumerate(q): out[i+j] += a*b
    return add(tuple(out))


def derivative(p):
    return add(tuple(i*p[i] for i in range(1, len(p))))


def coefficients(word, h):
    """Return P_k with derivative = q^beta sum P_k/L^(2n) s^k chi^(k)."""
    if isinstance(h, (float, bool)) or not 0 <= Q(h) <= Q(1, 100):
        raise ValueError('exact 0<=h<=1/100 required')
    if not isinstance(word, str) or any(x not in 'zt' for x in word):
        raise ValueError('word must contain only z and t')
    h = Q(h); D = Q(1, 2)-h
    L = (Q(1), Q(0), -2*h)
    Lp = derivative(L); eta = (Q(0), Q(1))
    beta = Q(0); rows = [(Q(1),)]
    for n, op in enumerate(word):
        out = []
        for k in range(len(rows)+1):
            P = rows[k] if k < len(rows) else ()
            prev = rows[k-1] if k else ()
            deriv = add(mul(derivative(P), L), scale(mul(P, Lp), -2*n))
            if op == 'z':
                out.append(add(scale(mul(mul(eta, P), L), 2*(beta+k)),
                               mul((Q(1), Q(0), Q(-1)), deriv),
                               scale(mul(mul(eta, prev), L), 2)))
            else:
                out.append(add(scale(mul(P, L), -(beta+k)),
                               scale(mul(eta, deriv), D), scale(mul(prev, L), -1)))
        rows = out
        beta -= D if op == 'z' else 1
    return beta, tuple(rows)


def coefficient_bound(word, h):
    _, rows = coefficients(word, h)
    return sum(sum(abs(v) for v in row) for row in rows) / Q(49, 50)**(2*len(word))


def evaluate(poly, eta):
    return sum(v*Q(eta)**i for i, v in enumerate(poly))
