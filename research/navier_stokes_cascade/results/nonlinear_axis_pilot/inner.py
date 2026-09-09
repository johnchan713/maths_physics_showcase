"""Triangular Taylor construction of (B.15), with a labelled pressure seed.

Each list holds ordinary Taylor coefficients, not derivatives: a[k] is
f^(k)(eta0)/k!. Radial order n consumes one available eta derivative.
Arbitrary-exponent mpmath numbers preserve strictly positive, extremely small g.
This finite polynomial is an approximate nonlinear profile, not an existence proof.
"""
from dataclasses import dataclass
import mpmath as mp


def cut(a, m):
    return list(a[:m+1]) + [mp.mpf(0)]*max(0, m+1-len(a))


def add(*arrays, m):
    return [mp.fsum(a[k] if k < len(a) else 0 for a in arrays) for k in range(m+1)]


def scale(a, c, m):
    return [c*x for x in cut(a, m)]


def mul(a, b, m):
    return [mp.fsum(a[i]*b[k-i] for i in range(max(0, k-len(b)+1), min(k+1, len(a))))
            for k in range(m+1)]


def inv(a, m):
    if not a or not a[0]:
        raise ValueError("Nonzero constant coefficient required")
    b = [1/a[0]]
    for k in range(1, m+1):
        b.append(-mp.fsum(a[i]*b[k-i] for i in range(1, min(k+1, len(a))))/a[0])
    return b


def derivative(a):
    return [(i+1)*a[i+1] for i in range(len(a)-1)]


def product_n(a, b, n, m, right_derivative=False, right_radial_weight=False):
    terms = []
    for i in range(n+1):
        if i >= len(a) or n-i >= len(b):
            continue
        right = derivative(b[n-i]) if right_derivative else b[n-i]
        value = mul(a[i], right, m)
        terms.append(scale(value, n-i, m) if right_radial_weight else value)
    return add(*terms, m=m)


@dataclass(frozen=True)
class Parameters:
    h: object
    j: object
    sigma: object
    lam: object
    pressure: object

    @classmethod
    def from_dict(cls, p):
        values = [mp.mpf(p[k]) for k in ("h", "j", "sigma", "Lambda", "pressure_scale")]
        if not all(mp.isfinite(v) for v in values):
            raise ValueError("Finite parameters required")
        h, j, sigma, lam, pressure = values
        if not (0 < h < mp.mpf('.01') and 0 < j <= mp.mpf('.05') and 0 < sigma < 1 and lam >= 1 and pressure > 0):
            raise ValueError("Parameters outside pilot domain")
        return cls(*values)

    @property
    def A(self):
        return mp.mpf('.5')+self.h

    @property
    def D(self):
        return mp.mpf('.5')-self.h

    def H(self, eta):
        return self.D*eta+(1-eta*eta)*(4*eta+self.j)

    def zeta(self, eta):
        hh = self.H(eta)
        return -(1-2*self.h*eta*eta)*hh/(hh*hh+self.sigma*self.sigma)

    def hzero(self):
        return mp.findroot(self.H, (-self.j, mp.mpf(0)))

    def log_c(self):
        # On |Re eta|<=1+rho, |Im eta|<=rho, use |eta|<2.
        # rho=sigma/[4*(D+4+48+4j)] bounds the change in H by sigma/4.
        # Thus |H +/- i sigma|>=3sigma/4, and sigma^2/4 is conservative.
        bound = 4*(1+8*self.h)*(2*self.D+5*(8+self.j))/self.sigma**2
        return self.lam*(2*bound+1)

    def log_g(self, eta):
        # Split at the narrow rational feature; quadrature remains numerical.
        root = self.hzero()
        lo, hi = min(mp.mpf(0), eta), max(mp.mpf(0), eta)
        points = [lo]+[v for v in (root-10*self.sigma, root, root+10*self.sigma) if lo < v < hi]+[hi]
        points = sorted(set(points))
        integral = mp.quad(self.zeta, points) if lo != hi else mp.mpf(0)
        return self.lam*(integral if eta >= 0 else -integral)-self.log_c()


def axis_series(p, eta0, order):
    m = order
    one = cut([mp.mpf(1)], m)
    eta = cut([eta0, mp.mpf(1)], m)
    eta2 = mul(eta, eta, m)
    d = add(one, scale(eta2, -1, m), m=m)
    L = add(one, scale(eta2, -2*p.h, m), m=m)
    U = add(scale(eta, 4, m), [p.j], m=m)
    H = add(scale(eta, p.D, m), mul(d, U, m), m=m)
    W = add(one, scale(d, -4, m), scale(mul(eta, U, m), -2*p.D, m), m=m)
    denominator = add(mul(H, H, m), [p.sigma**2], m=m)
    zeta = scale(mul(mul(L, H, m), inv(denominator, m), m), -1, m)
    chi = mul(mul(H, H, m), inv(denominator, m), m)
    f = inv(add(one, eta2, m=m), m+1)
    pressure = scale(mul(f, f, m+1), -p.pressure, m+1)
    Z = add(scale(mul(add(one, scale(mul(eta, U, m), -2, m), m=m), U, m), -p.A, m),
            scale(H, -4, m), scale(mul(d, derivative(pressure), m), -1, m),
            scale(mul(eta, pressure, m), 4*p.A, m), m=m)
    return dict(eta=eta, d=d, L=L, U=U, H=H, W=W, zeta=zeta, chi=chi,
                pressure=pressure, Z=Z, inverse_L=inv(L, m))


@dataclass
class InnerProfile:
    parameters: Parameters
    eta0: object
    degree: int
    phi: list
    u: list
    pressure_increment: list
    log_g0: object

    def value(self, coefficients, Y, eta=None, radial_order=0, eta_order=0):
        delta = mp.mpf(0) if eta is None else eta-self.eta0
        terms = []
        for n in range(radial_order, len(coefficients)):
            row = coefficients[n]
            angular = mp.fsum(row[k]*mp.factorial(k)/mp.factorial(k-eta_order)*delta**(k-eta_order)
                             for k in range(eta_order, len(row)))
            terms.append(angular*mp.factorial(n)/mp.factorial(n-radial_order)*Y**(n-radial_order))
        return mp.fsum(terms)

    def point(self, Y, eta=None):
        eta = self.eta0 if eta is None else eta
        p = self.parameters
        phi = self.value(self.phi, Y, eta)
        uu = self.value(self.u, Y, eta)
        logg = self.log_g0 if eta == self.eta0 else p.log_g(eta)
        g = mp.exp(logg)
        return phi, 4*eta+p.j+uu/p.lam, -p.pressure/(1+eta*eta)**2+self.value(self.pressure_increment, Y, eta)/p.lam, g


def construct(p, eta0, degree, extra_eta=3):
    if not isinstance(degree, int) or not 2 <= degree <= 40 or not 2 <= extra_eta <= 6 or abs(eta0) > 1:
        raise ValueError("Invalid degree, derivative budget, or eta sample")
    total = degree+extra_eta
    base = axis_series(p, eta0, total)
    phi = [cut([mp.mpf(1)], total)]
    u = [cut([mp.mpf(0)], total)]
    pressure = [cut([mp.mpf(0)], total)]
    logg = p.log_g(eta0)
    g2 = [mp.exp(2*logg)]
    # The exponential coefficient recurrence uses (g^2)'=2 Lambda zeta g^2.
    for k in range(total):
        g2.append(2*p.lam*mp.fsum(base['zeta'][i]*g2[k-i] for i in range(k+1))/(k+1))
    for n in range(degree):
        m = total-n-1
        eta, d, U0, H, W, zeta, L1 = [cut(base[k], m) for k in ('eta','d','U','H','W','zeta','inverse_L')]
        B = [scale(add(scale(mul(eta, row, m), 2*p.D, m), mul(d, derivative(row), m), m=m), -1/mp.mpf(i+1), m)
             for i, row in enumerate(u)]
        uphi = product_n(u, phi, n, m)
        uphieta = product_n(u, phi, n, m, right_derivative=True)
        bphi = add(product_n(B, phi, n, m), product_n(B, phi, n, m, right_radial_weight=True), m=m)
        R1 = add(scale(mul(W, phi[n], m), n+1, m), scale(bphi, 1/p.lam, m),
                 scale(mul(add([mp.mpf(1)], scale(mul(eta,U0,m),-2,m),m=m),phi[n],m),p.h,m),
                 mul(add(mul(d,zeta,m),scale(eta,-2*p.h/p.lam,m),m=m),uphi,m),
                 mul(H,derivative(phi[n]),m),scale(mul(d,uphieta,m),1/p.lam,m),m=m)
        R1 = mul(L1, R1, m)
        uu = product_n(u, u, n, m)
        uueta = product_n(u, u, n, m, right_derivative=True)
        bu = product_n(B, u, n, m, right_radial_weight=True)
        K = add(scale(add([mp.mpf(1)], scale(mul(eta,U0,m),-4,m),m=m),p.A,m),scale(d,4,m),m=m)
        R2 = add(mul(K,u[n],m),scale(mul(eta,uu,m),-2*p.A/p.lam,m),
                 scale(mul(W,u[n],m),n,m),scale(bu,1/p.lam,m),mul(H,derivative(u[n]),m),
                 scale(mul(d,uueta,m),1/p.lam,m),scale(mul(eta,pressure[n],m),-4*p.A-2*n,m),
                 mul(d,derivative(pressure[n]),m),m=m)
        R2 = mul(L1,R2,m)
        phi.append(scale(add(scale(mul(base['chi'],phi[n],m),-1,m),scale(R1,1/p.lam,m),m=m),
                         1/mp.mpf(2*(n+1)*(n+2)),m))
        force0 = scale(mul(base['Z'],L1,m),-1,m) if n == 0 else []
        u.append(scale(add(force0,scale(R2,1/p.lam,m),m=m),1/mp.mpf(2*(n+1)**2),m))
        pressure.append(scale(mul(g2,product_n(phi,phi,n,m),m),1/mp.mpf(n+1),m))
    return InnerProfile(p, eta0, degree, phi, u, pressure, logg)


def diagnostics(profile, Y):
    """Evaluate original (4.9)/(B.15), not the recursion's remainder formulas."""
    p, eta = profile.parameters, profile.eta0
    A, D, lam = p.A, p.D, p.lam
    L, d = 1-2*p.h*eta**2, 1-eta**2
    phi, U, Pi, g = profile.point(Y)
    val = profile.value
    phiy, phiyy, phie = val(profile.phi,Y,radial_order=1), val(profile.phi,Y,radial_order=2), val(profile.phi,Y,eta_order=1)
    uy, uyy, ue = val(profile.u,Y,radial_order=1), val(profile.u,Y,radial_order=2), val(profile.u,Y,eta_order=1)
    avg = [scale(row,1/mp.mpf(n+1),len(row)-1) for n,row in enumerate(profile.u)]
    average_U = 4*eta+p.j+val(avg,Y)/lam
    average_U_eta = 4+val(avg,Y,eta_order=1)/lam
    W = 1-2*D*eta*average_U-d*average_U_eta
    Hc = D*eta+d*U
    ell = 1+Y*phiy/phi
    logEeta = lam*p.zeta(eta)+phie/phi
    Pi_eta = 4*p.pressure*eta/(1+eta**2)**3+val(profile.pressure_increment,Y,eta_order=1)/lam
    PiX = val(profile.pressure_increment,Y,radial_order=1)
    sq_terms = [-W*ell,-p.h*(1-2*eta*U),-Hc*logEeta]
    sn_terms = [-W*Y*uy/lam,-A*(1-2*eta*U)*U,-Hc*(4+ue/lam),-d*Pi_eta,4*A*eta*Pi,2*eta*Y*PiX/lam]
    theta_lhs = -2*L*lam*(Y*phiyy+2*phiy)/phi
    axial_lhs = -2*L*(Y*uyy+uy)
    theta_scale = max(mp.mpf(1),abs(theta_lhs),*(abs(t) for t in sq_terms))
    axial_scale = max(mp.mpf(1),abs(axial_lhs),*(abs(t) for t in sn_terms))
    # Pressure is tiny; normalize by its own nonzero magnitude, with no floor.
    pressure_rhs = (g*phi)**2
    pressure_error = abs(PiX-pressure_rhs)/max(abs(PiX),abs(pressure_rhs))
    X = Y/lam
    E = mp.sqrt(2*X)*g*phi if Y else mp.mpf(0)
    moments = inner_moments(profile,Y)
    return dict(Y=Y,Phi=phi,U=U,Pi=Pi,log_g=profile.log_g0,
                angular_residual=abs(theta_lhs-mp.fsum(sq_terms))/theta_scale,
                axial_residual=abs(axial_lhs-mp.fsum(sn_terms))/axial_scale,
                pressure_residual=pressure_error,
                omitted_eta_transport_gap=abs(Hc*logEeta)/theta_scale,
                omitted_axial_pressure_gap=abs(-d*Pi_eta+4*A*eta*Pi+2*eta*Y*PiX/lam)/axial_scale,
                angular_slope=(mp.mpf('.5')+Y*phiy/phi),
                a=-2*Y*phiy/phi, E=E, moments=moments,
                axis_U_error=abs(U-(4*eta+p.j)) if not Y else None)


def inner_moments(profile,Y):
    p = profile.parameters
    phi = [r[0] for r in profile.phi]
    U = [r[0]/p.lam for r in profile.u]
    U[0] += 4*profile.eta0+p.j
    g = mp.exp(profile.log_g0)
    def integral(coeff, extra=0):
        return mp.fsum(c*Y**(n+extra+1)/(n+extra+1) for n,c in enumerate(coeff))
    def convolution(a,b):
        return mul(a,b,len(a)+len(b)-2)
    return [integral(U)/p.lam,
            2*g*integral(phi,1)/p.lam**2,
            2*g*integral(convolution(U,phi),1)/p.lam**2,
            integral(convolution(U,U))/p.lam-g*g*integral(convolution(phi,phi),1)/p.lam**2,
            g*g*integral(convolution(phi,phi))/p.lam]


def direct_heat_join(profile,Y=mp.mpf(4)):
    """Necessary conditions only. Exact heat exterior has U=V0=0.

    Its logarithmic E slope is -A-Z H'/H in (-A,-1/2); at eta=+/-1
    it equals -A. Amplitude fitting cannot change this slope or eliminate U.
    """
    row = diagnostics(profile,Y)
    p = profile.parameters
    slope = row['angular_slope']
    gap = max(-p.A-slope, slope+mp.mpf('.5'), mp.mpf(0))
    return dict(U_gap=abs(row['U'])/max(1,abs(4*profile.eta0+p.j)),
                M_gap=abs(row['moments'][0])/(Y/p.lam*max(1,abs(4*profile.eta0+p.j))),
                heat_slope_necessary_gap=gap,
                matched=False, reason="Direct joining requires zero axial flow, zero radial-average axial flow, and a compatible swirl slope; a collar has not been constructed.")
