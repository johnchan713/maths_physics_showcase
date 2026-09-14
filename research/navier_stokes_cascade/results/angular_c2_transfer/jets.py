"""Value, first derivative, and RAW second derivative bookkeeping.

These finite jets check differentiation identities, not global norm bounds.
Fractions stay exact. In particular dd is f'', not the Taylor coefficient f''/2.
"""
from dataclasses import dataclass
import mpmath as mp


@dataclass(frozen=True)
class Jet2:
    v: object
    d: object = 0
    dd: object = 0

    @staticmethod
    def cast(value):
        return value if isinstance(value, Jet2) else Jet2(value)

    def __add__(self, other):
        b = self.cast(other)
        return Jet2(self.v+b.v, self.d+b.d, self.dd+b.dd)

    __radd__ = __add__

    def __neg__(self):
        return Jet2(-self.v, -self.d, -self.dd)

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return self.cast(other)+-self

    def __mul__(self, other):
        b = self.cast(other)
        return Jet2(self.v*b.v, self.d*b.v+self.v*b.d,
                    self.dd*b.v+2*self.d*b.d+self.v*b.dd)

    __rmul__ = __mul__

    def inverse(self):
        if self.v == 0:
            raise ZeroDivisionError('A zero denominator cannot be normalized')
        return Jet2(1/self.v, -self.d/self.v**2,
                    (2*self.d**2-self.v*self.dd)/self.v**3)

    def __truediv__(self, other):
        return self*self.cast(other).inverse()

    def __rtruediv__(self, other):
        return self.cast(other)*self.inverse()

    def exp(self):
        value = mp.exp(self.v)
        # Differentiating exp(f) twice retains BOTH f'' and (f')^2.
        return Jet2(value, value*self.d, value*(self.dd+self.d*self.d))

    def log(self):
        if self.v <= 0:
            raise ValueError('The real logarithm needs a positive value')
        return Jet2(mp.log(self.v), self.d/self.v,
                    self.dd/self.v-(self.d/self.v)**2)


def physical_to_rows(moment_defects, eta, pressure, radius, freeze_normalization=False):
    """The original five rows with every angular normalization derivative.

    Input order is M,I,J,S,Cp. Bare numbers are rejected because their angular
    derivatives are unknown. The frozen-normalization switch is ONLY a failing
    diagnostic; its output cannot establish any actual bound.
    """
    if len(moment_defects) != 5 or not all(isinstance(v, Jet2) for v in moment_defects):
        raise ValueError('All five value/first/second derivative jets are required')
    if not mp.isfinite(pressure) or not mp.isfinite(radius) or pressure <= 0 or radius <= 0:
        raise ValueError('Finite positive physical scales required')
    if not mp.isfinite(eta) or abs(eta) > 1:
        raise ValueError('The angular coordinate must lie in [-1,1]')
    e = Jet2(eta, 1, 0)
    f = 1/(1+e*e)
    if freeze_normalization:
        f = Jet2(f.v)
    M, I, J, S, Cp = moment_defects
    xc = mp.exp(-6)
    angular = mp.sqrt(2)*radius**mp.mpf('1.5')*pressure*f*xc**mp.mpf('1.6')
    return [M/(radius*xc), (J-4*e*I)/angular, I/angular,
            -(S-8*e*M)/(radius*pressure**2*f*f*xc**mp.mpf('1.2')),
            Cp/(pressure**2*f*f*xc**mp.mpf('.2'))]
