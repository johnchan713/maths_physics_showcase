"""Exact heat-tail definitions and independent signed/positive stress formulas."""
from functools import lru_cache
import mpmath as mp


def step(x):
    x = mp.mpf(x)
    if x <= 0:
        return mp.mpf(0)
    if x >= 1:
        return mp.mpf(1)
    a, b = mp.exp(-1/x**2), mp.exp(-1/(1-x)**2)
    return a/(a+b)


def step_jets(x):
    x = mp.mpf(x)
    if not 0 < x < 1:
        return step(x), mp.mpf(0), mp.mpf(0)
    s, complement = step(x), step(1-x)
    first_log = 2/x**3+2/(1-x)**3
    second_log = -6/x**4+6/(1-x)**4
    product = s*complement
    return s, product*first_log, product*((1-2*s)*first_log**2+second_log)


def heat_derivative(h, Z, j=0):
    """Exact H^(j) via Tricomi U; the continuum bounds use the integral."""
    h, Z = map(mp.mpf, (h, Z))
    if not 0 < h < mp.mpf('.01') or Z < 0 or type(j) is not int or j < 0:
        raise ValueError('Construction h, nonnegative Z and integer derivative order required')
    coefficient = (-1)**j*mp.rf(h, j)*mp.rf(1+h, j)
    if Z == 0:
        return coefficient
    return coefficient*Z**(-1-h-j)*mp.hyperu(1+h+j, 2, 1/Z)


def axial_coefficients(h, eta):
    h, eta = map(mp.mpf, (h, eta))
    if not 0 < h < mp.mpf('.01') or abs(eta) > 1:
        raise ValueError('Construction h and closed eta interval required')
    L, d, D = 1-2*h*eta*eta, 1-eta*eta, mp.mpf('.5')-h
    alpha = 4*eta*eta/L**2
    beta = -2*(d-2*D*eta*eta)/L**2-8*h*eta*eta*d/L**3
    return L, alpha, beta


def positive_coefficient(h, eta, tH):
    h, eta, tH = map(mp.mpf, (h, eta, tH))
    L, alpha, _ = axial_coefficients(h, eta)
    return 2/L**2+8*h*eta*eta*(1-eta*eta)/L**3+alpha*tH


class HeatTail:
    """Manufacturable parameters of the SAME formula; not the giant global field."""
    def __init__(self, h='.005', xtail=1000, amplitude='1e-10', order=24, panels=8):
        self.h, self.xtail, self.amplitude = map(mp.mpf, (h, xtail, amplitude))
        if not 0 < self.h < mp.mpf('.01') or self.xtail < 100 or self.amplitude <= 0:
            raise ValueError('Positive fixed tail with Xtail>=100 and construction h required')
        if type(order) is not int or order < 8 or type(panels) is not int or panels < 1:
            raise ValueError('Positive quadrature resolution required')
        self.A, self.D, self.rho = mp.mpf('.5')+self.h, mp.mpf('.5')-self.h, self.h/1000
        self.xb = self.xtail*mp.exp(3)
        x, w = mp.gauss_quadrature(order, 'legendre')
        self.nodes = [((i+(1+x[j])/2)/panels, w[j]/(2*panels))
                      for i in range(panels) for j in range(order)]
        x, w = mp.gauss_quadrature(max(24, order), 'laguerre')
        self.laguerre = list(zip(x, w))

    @lru_cache(maxsize=32768)
    def state(self, y, eta):
        y, eta = map(mp.mpf, (y, eta))
        if not mp.mpf('.5') <= y <= 3 or abs(eta) > 1:
            raise ValueError('Use the closed terminal interval and closed eta interval')
        X = self.xtail*mp.exp(y)
        r, Z = mp.sqrt(2*X), 2*(1-eta*eta)/X
        H, Hp = heat_derivative(self.h, Z), heat_derivative(self.h, Z, 1)
        tH = -Z*Hp/H
        K = self.amplitude*mp.exp(-self.A*y)*H
        sigma, first, second = step_jets((3-y)/2)
        f, fp, fpp = 1-self.rho*sigma, self.rho*first/2, -self.rho*second/4
        a = 2+2*self.h-2*tH-2*fp/f
        return dict(X=X, r=r, Z=Z, H=H, Hp=Hp, tH=tH, K=K,
                    f=f, fp=fp, fpp=fpp, a=a, radial=2+2*self.h-2*tH)

    def integrate(self, fn, left):
        # All nonzero cutoff derivatives are supported in 1<v<3.
        left = max(mp.mpf(1), mp.mpf(left))
        return (3-left)*mp.fsum(w*fn(left+(3-left)*x) for x, w in self.nodes)

    def direct(self, y, eta):
        """Signed axial-viscous source, and the positive representation after parts."""
        y, eta = map(mp.mpf, (y, eta))
        if y >= 3:
            raise ValueError('Raw stress ratios are undefined at the zero edge')
        s = self.state(y, eta)
        L, alpha, beta = axial_coefficients(self.h, eta)
        def density(v, kind):
            t = self.state(v, eta)
            if kind == 'leading_radial':
                return t['r']*t['K']*t['radial']*t['fp']
            Hrad = t['r']**3*t['K']
            if kind == 'leading_time':
                return Hrad*t['fp']
            if kind == 'axial':
                return s['r']**2*mp.expm1(v-y)*t['K']**2*t['f']*t['fp']
            if kind == 'signed':
                return Hrad*(alpha*t['fpp']+beta*t['fp'])
            if kind == 'positive':
                return Hrad*positive_coefficient(self.h, eta, t['tH'])*t['fp']
            raise ValueError('Unknown integrand')
        values = {name: self.integrate(lambda v: density(v, name), y)
                  for name in ('leading_radial', 'leading_time', 'axial', 'signed', 'positive')}
        T0theta = (2*s['K']/s['r']*s['fp']+values['leading_radial']/s['r']**2
                   +values['leading_time']/(2*L*s['r']**2))
        T0z = eta*values['axial']/(L*s['r'])
        T1signed = -values['signed']/(2*s['r']**2)
        T1positive = alpha*s['r']*s['K']*s['fp']/2+values['positive']/(2*s['r']**2)
        return dict(T0theta=T0theta, T0z=T0z, T1signed=T1signed, T1positive=T1positive,
                    T1z=mp.mpf(0), a=s['a'], K=s['K'], correction_ratio=T1positive/T0theta,
                    leading_ratio=T0z/T0theta,
                    uniform_ratio_bound=2*self.xb/mp.mpf('.98')**2+7)

    @staticmethod
    def flat_kernel(u):
        u = mp.mpf(u)
        if u == 0:
            return 8*mp.e
        if not 0 < u < 2:
            raise ValueError('The flat coordinate must lie in [0,2)')
        c = 1-u/2
        a, b = mp.exp(-1/c**2), mp.exp(-4/u**2)
        return a*(8+u**3/c**3)/(a+b)**2

    def factored(self, delta, eta):
        """Integrate after removing the flat exponential; no tiny 0/0 ratio."""
        delta, eta = map(mp.mpf, (delta, eta))
        if not 0 < delta <= mp.mpf('.5') or abs(eta) > 1:
            raise ValueError('Positive collar delta<=.5 and closed eta interval required')
        y, L = 3-delta, axial_coefficients(self.h, eta)[0]
        alpha = axial_coefficients(self.h, eta)[1]
        s = self.state(y, eta)
        first, second, axial, viscous = [mp.mpf(0)]*4
        for w, weight in self.laguerre:
            root = mp.sqrt(1+delta**2*w/4)
            u = delta/root
            separation_scaled = (w/4)/(root*(1+root))
            separation = delta**3*separation_scaled
            exp_ratio = mp.expm1(separation)/separation if separation else 1
            t, Bflat = self.state(3-u, eta), self.flat_kernel(u)
            first += weight*Bflat*t['r']*t['K']*t['radial']
            second += weight*Bflat*t['r']**3*t['K']
            axial += weight*Bflat*exp_ratio*separation_scaled*t['K']**2*t['f']
            viscous += weight*Bflat*t['r']**3*t['K']*positive_coefficient(self.h, eta, t['tH'])
        b0 = (2*self.rho*s['K']/s['r']*self.flat_kernel(delta)
              +self.rho*delta**3/(8*s['r']**2)*(first+second/(2*L)))
        bz = eta*self.rho*s['r']*axial/(8*L)
        b1 = (alpha*self.rho*s['r']*s['K']*self.flat_kernel(delta)/2
              +self.rho*delta**3*viscous/(16*s['r']**2))
        return dict(b0=b0, bz=bz, b1=b1, correction_ratio=b1/b0,
                    leading_ratio=delta**6*bz/b0,
                    flat_log=-4/delta**2, float_flat_factor=float(mp.exp(-4/delta**2)))

    def endpoint(self, eta):
        eta = mp.mpf(eta)
        edge = self.state(mp.mpf(3), eta)
        L, alpha, _ = axial_coefficients(self.h, eta)
        return dict(b0=16*self.rho*mp.e*edge['K']/edge['r'],
                    bz=eta*self.rho*mp.e*edge['r']*edge['K']**2/(8*L),
                    b1=4*alpha*self.rho*mp.e*edge['r']*edge['K'],
                    correction_ratio=2*eta*eta*self.xb/L**2,
                    direction=[1, 0], stress=[0, 0])
