"""Basis-free positive jet budgets for a conditional first forced inverse."""
from functools import lru_cache
from fractions import Fraction
from math import comb
from pathlib import Path
import importlib.util

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("forced_inverse_roots", HERE.parent / "core_primary_root_budget/budgets.py")
roots = importlib.util.module_from_spec(spec)
spec.loader.exec_module(roots)
P, order, exact = roots.P, roots.order, roots.exact


def product(a, b, n):
    """Raw mixed derivative maximum: sum of multibinomial coefficients by size."""
    order(n)
    return sum((comb(n, i) * a(i) * b(n-i) for i in range(n+1)), P())


class ForcedBudget:
    """Inputs bound full-path coefficient jets; they are not actual fields."""
    def __init__(self, normal, moving_normal, matrix, source, profile_rate, tangent, length,
                 inverse_normal_square, harmonic=1, trace=None):
        if type(harmonic) is not int or harmonic == 0:
            raise ValueError("A nonzero integer input harmonic is required")
        self.harmonic = abs(harmonic)
        self.normal, self.moving_normal = normal, moving_normal
        self.matrix, self.source = matrix, source
        self.profile_rate = profile_rate
        self.tangent, self.length = tangent, length
        self.inverse_normal_square = exact(inverse_normal_square)
        if self.inverse_normal_square <= 0:
            raise ValueError("A positive full normal-square reciprocal floor is required")
        for p in (tangent, length):
            if not isinstance(p, P):
                raise TypeError("Positive scale-polynomial propagator and length required")
        if tangent.evaluate(1) < 1:
            raise ValueError("A propagator bound includes the identity at equal times")
        self.trace = trace

    def data(self, name, n):
        order(n)
        if self.trace is not None:
            self.trace.add((name, n))
        value = getattr(self, name)(n)
        if not isinstance(value, P):
            raise TypeError("Positive scale-polynomial jet bounds required")
        return value

    @lru_cache(None)
    def square(self, n):
        return product(lambda j:self.data("normal", j), lambda j:self.data("normal", j), n)

    @lru_cache(None)
    def reciprocal(self, n):
        order(n)
        if not n:
            return P(self.inverse_normal_square)
        return self.inverse_normal_square * sum(
            (comb(n, i)*self.square(i)*self.reciprocal(n-i) for i in range(1, n+1)), P())

    @lru_cache(None)
    def lift(self, n):
        return product(lambda j:self.data("normal", j), self.reciprocal, n)

    @lru_cache(None)
    def pressure_row(self, n):
        """n^T K - n'^T; the moving normal is not dropped."""
        return product(lambda j:self.data("normal", j), lambda j:self.data("matrix", j), n) + self.data("moving_normal", n)

    @lru_cache(None)
    def undamped(self, n):
        return self.data("matrix", n) + product(self.lift, self.pressure_row, n)

    @lru_cache(None)
    def damping(self, n):
        # Fixed epsilon and k per lifted label, epsilon*k^2 <= 4.
        return 4*self.harmonic**2*self.square(n)

    @lru_cache(None)
    def normal_projector(self, n):
        return product(self.lift, lambda j:self.data("normal", j), n)

    @lru_cache(None)
    def generator(self, n):
        # Gaussian-adapted ambient extension: Ahat=A0-lambda*I
        # +(lambda+P'/P)*h*n^T. It agrees with A on tangent vectors.
        return self.undamped(n) + self.damping(n) + product(
            lambda j:self.damping(j)+self.data("profile_rate", j), self.normal_projector, n)

    @lru_cache(None)
    def projection(self, n):
        order(n)
        if not n:
            return P(1)  # orthogonal projection, no component floor
        return self.normal_projector(n)

    @lru_cache(None)
    def projected_source(self, n):
        return product(self.projection, lambda j:self.data("source", j), n)

    @lru_cache(None)
    def ambient(self):
        # Q = A0 h - h', with h=n/|n|^2. It is tangent.
        Q = self.undamped(0)*self.lift(0) + self.lift(1)
        return self.tangent*(P(1) + self.length*Q*self.data("normal", 0)) + self.lift(0)*self.data("normal", 0)

    @lru_cache(None)
    def velocity(self, n):
        order(n)
        parameter = self.ambient()*self.length*(self.projected_source(n) + sum(
            (comb(n, i)*self.generator(i)*self.velocity(n-i) for i in range(1, n+1)), P()))
        if not n:
            return parameter
        pulse = self.projected_source(n-1) + sum(
            (comb(n-1, i)*self.generator(i)*self.velocity(n-1-i) for i in range(n)), P())
        return parameter + pulse

    @lru_cache(None)
    def pressure(self, n):
        row = lambda j:product(self.reciprocal, self.pressure_row, j)
        return Fraction(1, self.harmonic)*(product(row, self.velocity, n) + product(self.lift, lambda j:self.data("source", j), n))

    @lru_cache(None)
    def potential(self, n):
        return Fraction(1, self.harmonic)*product(self.lift, self.velocity, n)

    def cut(self, kind, cutoff, n):
        if kind not in ("velocity", "pressure", "potential"):
            raise ValueError("Known coefficient kind required")
        return product(cutoff, getattr(self, kind), n)


def fixture(trace=None, harmonic=2):
    return ForcedBudget(normal=lambda n:P({0:n+2, 1:1}),
                        moving_normal=lambda n:P({0:n+1, 1:2}),
                        matrix=lambda n:P({0:n+3, 1:1}),
                        source=lambda n:P({0:n+1, 2:1}),
                        profile_rate=lambda n:P({0:n+1, 1:1}),
                        tangent=P({0:2, 1:1}), length=P({1:1}),
                        inverse_normal_square=2, harmonic=harmonic, trace=trace)
