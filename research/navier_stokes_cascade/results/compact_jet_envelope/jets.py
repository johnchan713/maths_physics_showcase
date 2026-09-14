"""Raw mixed jets. Stress outputs are FIRST order, not invented third derivatives."""
from dataclasses import dataclass
from fractions import Fraction
import mpmath as mp


def inverse_scalar(value):
    # Casting the integer divisor 5 and evaluating 1/5 in Python would silently
    # insert a binary float into otherwise high-precision or exact arithmetic.
    return Fraction(1)/value if isinstance(value,(int,Fraction)) else 1/value


@dataclass(frozen=True)
class FirstJet:
    v: object
    y: object = 0
    e: object = 0

    @staticmethod
    def cast(v):
        if isinstance(v, SecondJet):
            raise TypeError('Truncate a second jet explicitly with first()')
        return v if isinstance(v, FirstJet) else FirstJet(v)

    def __add__(self, other):
        b = self.cast(other)
        return FirstJet(self.v+b.v, self.y+b.y, self.e+b.e)

    __radd__ = __add__

    def __neg__(self):
        return FirstJet(-self.v, -self.y, -self.e)

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return self.cast(other)+-self

    def __mul__(self, other):
        b = self.cast(other)
        return FirstJet(self.v*b.v, self.y*b.v+self.v*b.y, self.e*b.v+self.v*b.e)

    __rmul__ = __mul__

    def inverse(self):
        if self.v == 0:
            raise ZeroDivisionError('Positive denominators cannot be replaced by zero')
        q = inverse_scalar(self.v)
        return FirstJet(q, -self.y*q*q, -self.e*q*q)

    def __truediv__(self, other):
        return self*self.cast(other).inverse()

    def __rtruediv__(self, other):
        return self.cast(other)*self.inverse()


@dataclass(frozen=True)
class SecondJet:
    v: object
    y: object = 0
    e: object = 0
    yy: object = 0
    ye: object = 0
    ee: object = 0

    @staticmethod
    def cast(v):
        if isinstance(v, FirstJet):
            raise TypeError('First-order data do not determine second derivatives')
        return v if isinstance(v, SecondJet) else SecondJet(v)

    def __add__(self, other):
        b = self.cast(other)
        return SecondJet(*(a+c for a, c in zip(self.values(), b.values())))

    __radd__ = __add__

    def values(self):
        return self.v, self.y, self.e, self.yy, self.ye, self.ee

    def __neg__(self):
        return SecondJet(*(-a for a in self.values()))

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return self.cast(other)+-self

    def __mul__(self, other):
        b = self.cast(other)
        return SecondJet(self.v*b.v, self.y*b.v+self.v*b.y, self.e*b.v+self.v*b.e,
                         self.yy*b.v+2*self.y*b.y+self.v*b.yy,
                         self.ye*b.v+self.y*b.e+self.e*b.y+self.v*b.ye,
                         self.ee*b.v+2*self.e*b.e+self.v*b.ee)

    __rmul__ = __mul__

    def compose(self, value, first, second):
        return SecondJet(value, first*self.y, first*self.e,
                         first*self.yy+second*self.y*self.y,
                         first*self.ye+second*self.y*self.e,
                         first*self.ee+second*self.e*self.e)

    def inverse(self):
        if self.v == 0:
            raise ZeroDivisionError('A zero denominator has no inverse jet')
        q = inverse_scalar(self.v)
        return self.compose(q, -q*q, 2*q*q*q)

    def __truediv__(self, other):
        return self*self.cast(other).inverse()

    def __rtruediv__(self, other):
        return self.cast(other)*self.inverse()

    def __pow__(self, n):
        if type(n) is not int:
            raise ValueError('Only exact integer powers are implemented')
        if n < 0:
            return self.inverse()**(-n)
        result = SecondJet(self.v*0+1)
        for _ in range(n):
            result = result*self
        return result

    def exp(self):
        value = mp.exp(self.v)
        return self.compose(value, value, value)

    def log(self):
        if self.v <= 0:
            raise ValueError('A real logarithm needs a positive argument')
        q = inverse_scalar(self.v)
        return self.compose(mp.log(self.v), q, -q*q)

    def first(self):
        return FirstJet(self.v, self.y, self.e)

    def partial_first(self, axis):
        if axis == 'y':
            return FirstJet(self.y, self.yy, self.ye)
        if axis == 'eta':
            return FirstJet(self.e, self.ye, self.ee)
        raise ValueError('Derivative axis must be y or eta')


def pressure_from_moments(y, eta, h, E, U, moments, axis_pressure):
    """(4.16), with all first output jets and the common axis datum retained."""
    if len(moments) != 5 or not all(isinstance(v, SecondJet) for v in
                                                    [E, U, axis_pressure, *moments]):
        raise ValueError('Fields, all five moments and axis pressure need second mixed jets')
    if not all(mp.isfinite(v) for v in (y, eta, h, E.v)):
        raise ValueError('Finite coordinates, h and E required')
    if abs(eta) > 1 or not 0 < h <= mp.mpf('.01') or E.v <= 0:
        raise ValueError('Outside the positive-field input domain')
    yj, ej = SecondJet(y, y=1), SecondJet(eta, e=1)
    X = yj.exp().first()
    e = ej.first()
    A, D = mp.mpf('.5')+h, mp.mpf('.5')-h
    d, L = 1-e*e, 1-2*h*e*e
    M, I, J, S, Cp = [v.first() for v in moments]
    Me, Ie, Je, Se, Cpe = [v.partial_first('eta') for v in moments]
    Ev, Uv = E.first(), U.first()
    H = mp.sqrt(2)*(yj/2).exp().first()*Ev
    Pi, Pie = axis_pressure.first()+Cp, axis_pressure.partial_first('eta')+Cpe
    W = 1-(2*D*e*M+d*Me)/X
    Q = -W+((1-h)*I-D*e*Ie-d*Je+2*(h-D)*e*J)/(X*H)
    N = -W*Uv+(D*(M-e*Me)+4*h*e*S-d*Se)/X+4*A*e*Pi-d*Pie
    return dict(W=W, Pi=Pi, Qs=Q, Ns=N, p1=X*Q/L, p2=X*N/(L*Ev))


def shear_from_fields(E, U):
    if not isinstance(E, SecondJet) or not isinstance(U, SecondJet) or E.v <= 0:
        raise ValueError('Positive E and second mixed field jets required')
    a = 1-2*E.partial_first('y')/E.first()
    bs = 2*U.partial_first('y')/E.first()
    if a.v <= 0:
        raise ValueError('The loop input requires positive a')
    ts = -bs/a
    return dict(a=a, bs=bs, ts=ts, vs=a*(1+ts*ts))
