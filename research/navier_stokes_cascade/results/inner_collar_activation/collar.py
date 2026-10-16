"""Flat Volterra quotients and the complete five-moment activation map.

Numerical seeds below are manufactured, NOT solutions of the analytic core.
All moment and pressure differences are formed without subtracting flat
numbers. The theorem and its exact-seed premise are in README.md.
"""
from dataclasses import dataclass
from functools import lru_cache
import mpmath as mp


@dataclass(frozen=True)
class Jet:
    """A value and its first eta derivative, with exact chain rules."""
    v: object
    d: object = 0

    def __post_init__(self):
        object.__setattr__(self, 'v', mp.mpf(self.v))
        object.__setattr__(self, 'd', mp.mpf(self.d))

    @staticmethod
    def cast(value):
        return value if isinstance(value, Jet) else Jet(value)

    def __add__(self, other):
        b = self.cast(other)
        return Jet(self.v+b.v, self.d+b.d)

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.v, -self.d)

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return self.cast(other)+-self

    def __mul__(self, other):
        b = self.cast(other)
        return Jet(self.v*b.v, self.d*b.v+self.v*b.d)

    __rmul__ = __mul__

    def __truediv__(self, other):
        b = self.cast(other)
        return Jet(self.v/b.v, (self.d*b.v-self.v*b.d)/(b.v*b.v))

    def __rtruediv__(self, other):
        return self.cast(other)/self

    def __pow__(self, n):
        return Jet(self.v**n, n*self.v**(n-1)*self.d) if n else Jet(1)


def exp(value):
    b = Jet.cast(value)
    return Jet(mp.exp(b.v), mp.exp(b.v)*b.d)


def exprel(value):
    """The entire function (exp(z)-1)/z, including its derivative at zero."""
    b = Jet.cast(value)
    if abs(b.v) < mp.mpf('.01'):
        term, total = Jet(1), Jet(1)
        for n in range(1, 45):
            term = term*b/(n+1)
            total += term
        return total
    return (exp(b)-1)/b


def flat_factor(s):
    s = mp.mpf(s)
    if not 0 <= s <= mp.mpf('.25'):
        raise ValueError('The inner quotient uses 0<=s<=1/4')
    return mp.exp(-1/s**2) if s else mp.mpf(0)


def ga(s):
    """sigma(s)=exp(-1/s^2)*ga(s), with ga(0)=e."""
    s = mp.mpf(s)
    return 1/(flat_factor(s)+mp.exp(-1/(1-s)**2))


def sigma(s):
    return flat_factor(s)*ga(s)


@lru_cache(maxsize=12)
def laguerre(order, digits):
    with mp.workdps(digits):
        nodes, weights = mp.gauss_quadrature(order, 'laguerre')
        return tuple(zip(nodes, weights))


def flat_quotient(s, coefficient, power=0, order=24):
    """Return int_0^s sigma(u)u^power*b(u)du / (sigma(s)s^(power+3)).

    A positive half-line integral supplies the smooth endpoint value b(0)/2.
    No raw division by a flat number is made, even below binary64 range.
    """
    s = mp.mpf(s)
    if power < 0 or int(power) != power:
        raise ValueError('A nonnegative integer power is required')
    flat_factor(s)
    if not s:
        return coefficient(mp.mpf(0))/2
    result = Jet(0) if isinstance(coefficient(s), Jet) else mp.mpf(0)
    for v, weight in laguerre(order, mp.mp.dps):
        ratio = 1+s*s*v
        u = s/mp.sqrt(ratio)
        result += weight*ratio**(-mp.mpf(power+3)/2)*ga(u)*coefficient(u)
    return result/(2*ga(s))


def seed(y, eta):
    """Closed-form manufactured fields and all five forward primitives.

    a=3, q=-2U_y/E=1+eta/5. Upstream moment values are manufactured.
    This seed deliberately has a nonzero leading-stress defect.
    """
    y, e = mp.mpf(y), Jet(eta, 1)
    x0, x = mp.mpf(2), 2*mp.exp(y)
    eb = 2+e/10+e**2/20
    q, ub = 1+e/5, e/3+e**2/30
    b, c = q*eb/2, ub-q*eb/2
    E, U = eb*mp.exp(-y), c+b*mp.exp(-y)
    H0, H = 2*eb, 2*eb*mp.exp(-y/2)
    M = (1+e/5)/10+x0*(c*mp.expm1(y)+b*y)
    I = (2+e/3)/10+2*x0*H0*mp.expm1(y/2)
    J = (1-e/5)/20+2*x0*H0*(c*mp.expm1(y/2)-b*mp.expm1(-y/2))
    S = (7+2*e**2)/100+x0*(c**2*mp.expm1(y)+2*c*b*y-(b**2-eb**2/2)*mp.expm1(-y))
    Cp = (3+e)/100-eb**2*mp.expm1(-2*y)/4
    Pi = -mp.mpf('.7')+e/5+3*e**2/20-eb**2*mp.expm1(-2*y)/4
    return dict(X=x, E=E, U=U, H=H, M=M, I=I, J=J, S=S, Cp=Cp, Pi=Pi, a=Jet(3), q=q)


def moment_map(state, h='.005', eta=0):
    """Original (4.16), including the pressure value AND eta derivative."""
    e, h = mp.mpf(eta), mp.mpf(h)
    A, D, d, L = mp.mpf('.5')+h, mp.mpf('.5')-h, 1-e*e, 1-2*h*e*e
    x, E, U, H = state['X'], state['E'].v, state['U'].v, state['H'].v
    M, I, J, S, Pi = [state[name] for name in ('M', 'I', 'J', 'S', 'Pi')]
    W = 1-(2*D*e*M.v+d*M.d)/x
    numerator = (1-h)*I.v-D*e*I.d-d*J.d+2*(h-D)*e*J.v
    Qs = -W+numerator/(x*H)
    Ns = -W*U+(D*(M.v-e*M.d)+4*h*e*S.v-d*S.d)/x+4*A*e*Pi.v-d*Pi.d
    return dict(W=W, numerator=numerator, Qs=Qs, Ns=Ns, p1=x*Qs/L, p2=x*Ns/(L*E))


class Activation:
    def __init__(self, width='.01', kappa0='.1', order=24):
        self.width, self.kappa0, self.order = mp.mpf(width), mp.mpf(kappa0), int(order)
        if not 0 < self.width <= mp.mpf('.01') or not 0 <= self.kappa0 <= mp.mpf('.5') or self.order < 8:
            raise ValueError('Invalid activation width, kappa0 or quadrature order')

    def field_coefficients(self, s, eta):
        """Coefficients after removing width*ea*s^3 from field differences."""
        s, t = mp.mpf(s), self.width
        base = seed(t*s, eta)
        ea = (1-self.kappa0)*sigma(s)
        lc = flat_quotient(s, lambda u: seed(t*u, eta)['a']/2, order=self.order)
        uc = flat_quotient(s, lambda u: seed(t*u, eta)['E']*seed(t*u, eta)['q']/2, order=self.order)
        w = t*s**3*lc
        ec = base['E']*lc*exprel(ea*w)
        e2c = 2*base['E']**2*lc*exprel(2*ea*w)
        return dict(base=base, ea=ea, lc=lc, uc=uc, ec=ec, e2c=e2c,
                    hc=mp.sqrt(2*base['X'])*ec, w=w,
                    E=base['E']*exp(ea*w), U=base['U']+ea*t*s**3*uc,
                    ratiohat=-w*exprel(-ea*w))

    @lru_cache(maxsize=512)
    def state(self, s, eta):
        s, eta, t = mp.mpf(s), mp.mpf(eta), self.width
        f = self.field_coefficients(s, eta)
        b, ea = f['base'], f['ea']

        def integrands(u):
            row = self.field_coefficients(u, eta)
            z, scale = row['base'], row['ea']*t*u**3
            uc, hc, e2c = row['uc'], row['hc'], row['e2c']
            return [z['X']*uc, z['X']*hc,
                    z['X']*(z['U']*hc+z['H']*uc+scale*uc*hc),
                    z['X']*(2*z['U']*uc+scale*uc**2-e2c/2), e2c/2]

        # Cache the five integrands together: all use the same positive rule.
        values = {}
        def entry(u, index):
            if u not in values:
                values[u] = integrands(u)
            return values[u][index]
        coefficients = [flat_quotient(s, lambda u, i=i: entry(u, i), power=3, order=self.order) for i in range(5)]
        hats = {name: t*t*s**6*c for name, c in zip(('M', 'I', 'J', 'S', 'Cp'), coefficients)}
        actual = dict(b)
        actual.update(E=f['E'], U=f['U'], H=mp.sqrt(2*b['X'])*f['E'])
        for name in ('M', 'I', 'J', 'S', 'Cp'):
            actual[name] = b[name]+ea*hats[name]
        actual['Pi'] = b['Pi']+ea*hats['Cp']
        return dict(fields=f, hats=hats, moment_coefficients=coefficients, actual=actual)

    def stress(self, s, eta):
        """The complete scaled ps difference and activated excess stress.

        The seed defect is NOT discarded: actual T = (F/Fr)Tr + ea*B.
        A cone test of B represents the exact-core theorem's conditional part.
        """
        s, e, t = mp.mpf(s), mp.mpf(eta), self.width
        state = self.state(s, e)
        f, dh = state['fields'], state['hats']
        b, ea, actual = f['base'], f['ea'], state['actual']
        ref = moment_map(b, eta=e)
        h = mp.mpf('.005')
        A, D, d, L = mp.mpf('.5')+h, mp.mpf('.5')-h, 1-e*e, 1-2*h*e*e
        x, E, Er = b['X'], f['E'].v, b['E'].v
        Eh, Uh, Hh = t*s**3*f['ec'], t*s**3*f['uc'], t*s**3*f['hc']
        Wh = -(2*D*e*dh['M'].v+d*dh['M'].d)/x
        Nh = (1-h)*dh['I'].v-D*e*dh['I'].d-d*dh['J'].d+2*(h-D)*e*dh['J'].v
        Qh = -Wh+Nh/(x*actual['H'].v)-ref['numerator']*Hh.v/(x*b['H'].v*actual['H'].v)
        pressure = 4*A*e*dh['Cp'].v-d*dh['Cp'].d
        Nsh = -(ref['W']*Uh.v+b['U'].v*Wh+ea*Wh*Uh.v) \
            +(D*(dh['M'].v-e*dh['M'].d)+4*h*e*dh['S'].v-d*dh['S'].d)/x+pressure
        p1h = x*Qh/L
        p2h = x/L*(Nsh/E-ref['Ns']*Eh.v/(Er*E))
        k = 1-ea
        a, q = b['a'].v, b['q'].v
        B = [a+p1h, q+p2h-k*q*f['ratiohat'].v]
        ts = q*mp.exp(-ea*f['w'].v)/a
        vs = k*a*(1+ts*ts)
        vr = a+q*q/a
        # Cancel first, then evaluate. Direct subtraction loses these terms
        # when ea lies far below the working precision relative to one.
        den_error = p1h+ts*p2h+q*q/a*ea*f['ratiohat'].v-k*q*f['ratiohat'].v*ts
        den = vr+den_error
        cross = p2h-ts*p1h-q*f['ratiohat'].v
        vs_error_over_ea = k*q*q/a*f['ratiohat'].v*(2+ea*f['ratiohat'].v)
        loss = max(vs-2, 0)*cross**2/(2*den**2)
        defect = [ref['p1']-a, ref['p2']-q]
        return dict(state=state, p1hat=p1h, p2hat=p2h, B_over_F=B,
                    ts=ts, vs=vs, denominator_over_ea=den, cross_over_ea=cross,
                    cone_loss=loss, seed_defect_over_Fr=defect,
                    denominator_error=den_error, vs_error_over_ea=vs_error_over_ea,
                    pressure_contribution_to_p2hat=x*pressure/(L*E),
                    angular_Jeta_contribution_to_p1hat=-d*dh['J'].d/(L*actual['H'].v))


def collar_threshold(vmin, vmax, error_bound, width_max='.01'):
    """An explicit sufficient dimensionless collar from the proved inequalities."""
    vmin, vmax, K, width = map(mp.mpf, (vmin, vmax, error_bound, width_max))
    if vmin <= 2 or vmax < vmin or K < 0 or width <= 0:
        raise ValueError('Positive seed growth and finite nonnegative bounds required')
    gamma = vmin-2
    allowance = min(1, vmin/2, vmin/(2*mp.sqrt(max(1, vmax-1))))
    algebra = (allowance/(K*width))**(mp.mpf(1)/3) if K else mp.inf
    activation = 1/mp.sqrt(max(16, mp.mpf(16)/9+mp.log(2*(vmax+1)/gamma)))
    return dict(smax=min(mp.mpf('.25'), algebra, activation), allowance=allowance,
                vmin=vmin, vmax=vmax, error_bound=K, width_max=width)
