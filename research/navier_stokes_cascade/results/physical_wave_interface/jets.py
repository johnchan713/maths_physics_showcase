"""Third-order Taylor arithmetic in four coordinates, with factorial coefficients.

Only derivatives through order two of a curl are used. Thus order three of
its vector potential suffices; no finite differences enter the residual.
"""
from itertools import product
import math
import numpy as np

DIM, ORDER = 4, 3
INDICES = sorted((a for a in product(range(ORDER+1), repeat=DIM) if sum(a) <= ORDER),
                 key=lambda a: (sum(a), a))
INDEX = {a: i for i, a in enumerate(INDICES)}
PAIRS = [(i, j, INDEX[tuple(x+y for x, y in zip(a, b))])
         for i, a in enumerate(INDICES) for j, b in enumerate(INDICES)
         if sum(a)+sum(b) <= ORDER]
LEFT, RIGHT, TARGET = (np.array(a, dtype=int) for a in zip(*PAIRS))


class Jet:
    __array_priority__ = 1000

    def __init__(self, value=0.0, coefficients=None):
        self.c = np.zeros(len(INDICES)) if coefficients is None else np.array(coefficients, dtype=float)
        if coefficients is None:
            self.c[0] = value

    @property
    def value(self):
        return float(self.c[0])

    @classmethod
    def variable(cls, value, axis):
        out = cls(value)
        a = [0]*DIM
        a[axis] = 1
        out.c[INDEX[tuple(a)]] = 1.0
        return out

    def __add__(self, other):
        other = coerce(other)
        return Jet(coefficients=self.c+other.c)

    __radd__ = __add__

    def __neg__(self):
        return Jet(coefficients=-self.c)

    def __sub__(self, other):
        return self+-coerce(other)

    def __rsub__(self, other):
        return coerce(other)+-self

    def __mul__(self, other):
        other = coerce(other)
        out = np.zeros(len(INDICES))
        np.add.at(out, TARGET, self.c[LEFT]*other.c[RIGHT])
        return Jet(coefficients=out)

    __rmul__ = __mul__

    def compose(self, derivatives):
        """derivatives[j] is the ordinary jth derivative at self.value."""
        delta = self-self.value
        power, out = Jet(1), Jet(derivatives[0])
        for j in range(1, ORDER+1):
            power = power*delta
            out = out+power*(derivatives[j]/math.factorial(j))
        return out

    def __pow__(self, exponent):
        x = self.value
        if x == 0 and exponent in (0, 1, 2, 3):
            out = Jet(1)
            for _ in range(exponent):
                out = out*self
            return out
        return self.compose([x**exponent]+[
            math.prod(exponent-i for i in range(j))*x**(exponent-j)
            for j in range(1, ORDER+1)])

    def __truediv__(self, other):
        return self*coerce(other)**-1

    def __rtruediv__(self, other):
        return coerce(other)*self**-1

    def diff(self, axis):
        out = np.zeros(len(INDICES))
        for i, a in enumerate(INDICES):
            b = list(a)
            b[axis] += 1
            if sum(b) <= ORDER:
                out[i] = b[axis]*self.c[INDEX[tuple(b)]]
        return Jet(coefficients=out)

    def exp(self):
        return self.compose([math.exp(self.value)]*4)

    def log(self):
        x = self.value
        return self.compose([math.log(x), 1/x, -1/x**2, 2/x**3])

    def sin(self):
        s, c = math.sin(self.value), math.cos(self.value)
        return self.compose([s, c, -s, -c])

    def cos(self):
        s, c = math.sin(self.value), math.cos(self.value)
        return self.compose([c, -s, -c, s])

    def atan(self):
        x = self.value
        return self.compose([math.atan(x), 1/(1+x*x), -2*x/(1+x*x)**2,
                             (6*x*x-2)/(1+x*x)**3])


def coerce(value):
    return value if isinstance(value, Jet) else Jet(value)


def variables(point):
    return [Jet.variable(x, i) for i, x in enumerate(point)]


def arrays(field):
    """Return values, gradients and Hessians with coordinate as the last index."""
    val = np.array([f.value for f in field])
    grad = np.array([[f.diff(i).value for i in range(DIM)] for f in field])
    hess = np.array([[[f.diff(i).diff(j).value for j in range(DIM)]
                     for i in range(DIM)] for f in field])
    return val, grad, hess


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def bump(s):
    """The standard smooth compact bump, extended by zero before derivatives."""
    if abs(s.value) >= 1:
        return Jet(0)
    return (-1/(1-s*s)).exp()
