"""Coefficient bounds and exact finite-series controls for the core theorem.

The all-index proof and existence argument are in README.md. Finite formal
jets test algebra, not convergence or the actual enormous-parameter solution.
"""
from fractions import Fraction as Q
from math import comb, factorial
import mpmath as mp


def weight(n, k, rho=Q(1)):
    if any(type(v) is not int or v < 0 for v in (n, k)) or rho <= 0:
        raise ValueError('Nonnegative integer indices and positive rho required')
    return Q(factorial(k)*comb(n+k, k), 20**n*(n+1)**2*(k+1)**2)/rho**k


def mixed_allocation(i, j, nu):
    if min(i, j) < 0 or nu not in (1, 2):
        raise ValueError('Invalid mixed derivative indices')
    return Q((i+1)*j, (i+j+1)*(i+j+nu))


def source_sums():
    """Termwise fixed-coefficient suprema, already including 1/L<=2."""
    return dict(R1_linear=2*(64+Q(37, 100))+128+128,
                R1_quadratic_without_zeta=4+10+Q(8, 100)+4+10+10,
                R1_zeta_multiplier=Q(10),
                R2_linear=2*(1+8*9+5*4)+128+128,
                R2_quadratic=8+4+10+10,
                R2_pressure=16+10+8)


def core_budget(scale='Md64', digits=80):
    """Independent outward arithmetic for ONLY the core, in logarithms."""
    if scale not in ('P16', 'Md4', 'Md64') or type(digits) is not int or digits < 60:
        raise ValueError('Unknown pressure scale or inadequate interval precision')
    mp.iv.dps = digits
    p, log, exp = mp.iv.mpf, mp.iv.log, mp.iv.exp
    upper = lambda v: mp.make_mpf(v._mpi_[1])
    lower = lambda v: mp.make_mpf(v._mpi_[0])

    def add(*values):
        shift = p(max(upper(v) for v in values))
        return shift+log(sum((exp(v-shift) for v in values), p(0)))

    lp = log(p(16)) if scale == 'P16' else exp(p(int(scale[2:])))+11
    lr0, lrho = log(p('1e-20'))-log(p(256)), log(p('1e-20'))-log(p(2048))
    la, ld, lm = log(p(256)), log(p(80))-lrho, log(p(200))-2*log(p('1e-20'))
    lz, lv = log(p(1024))+2*lp-lr0, p(30720)
    ls = add(log(p(2)), lv, log(p(40))+lz)
    q1 = add(log(p(40)), log(p(10))+lm)+2*la
    q2 = add(log(p(32))+2*la, log(p(2720))+4*la)
    m1, m2 = ld+add(log(p(386))+la+ls, q1+2*ls), ld+add(log(p(442))+la+ls, q2+2*ls)
    l1, l2 = ld+add(log(p(386))+la, log(p(2))+q1+ls), ld+add(log(p(442))+la, log(p(2))+q2+ls)
    lf, llip = lv+add(m1, m2), lv+add(l1, l2)
    leval = log(p(10000))-6*lrho
    lk = 16*(log(p('1e6'))+leval+add(p(0), ls, lf, llip, lm, lz, 2*lp))
    llam = log(p('1e12'))+8*lk-2*log(p('1e-16'))
    logs = dict(P=lp, rho=lrho, zeta_bound=lm, Z_over_L_bound=lz,
                angular_inverse_bound=lv, ball_norm_bound=ls, J_R1_bound=m1, J_R2_bound=m2,
                map_displacement_numerator=lf, map_Lipschitz_numerator=llip,
                real_evaluation_bound=leval, K=lk, Lambda=llam,
                map_displacement=lf-llam, map_Lipschitz=llip-llam,
                six_parameter_two_radial_error=leval+lf-llam,
                axial_entry_error=leval+ls-llam,
                minimum_log_C=llam+log(2*exp(lm)+1))
    ceiling = log(p('1e-50'))
    checks = dict(map_invariant=upper(logs['map_displacement']) < lower(ceiling),
                  strict_contraction=upper(logs['map_Lipschitz']) < lower(ceiling),
                  normalized_jet_error=upper(logs['six_parameter_two_radial_error']) < lower(ceiling),
                  axial_entry_error=upper(logs['axial_entry_error']) < lower(ceiling),
                  same_selected_Lambda_formula=True)
    return dict(scale=scale, digits=digits, logs=logs, checks=checks)


class Series:
    """Finite rational Taylor algebra in (Y,eta), explicitly a diagnostic.

    Caps are shared within a calculation. Eta padding is retained until the
    final low-order comparison; differentiation can consume this padding.
    """
    ny, nk = 5, 10

    def __init__(self, terms=None):
        self.terms = {tuple(key): Q(value) for key, value in (terms or {}).items()
                      if value and 0 <= key[0] <= self.ny and 0 <= key[1] <= self.nk}

    @classmethod
    def cast(cls, value):
        return value if isinstance(value, cls) else cls({(0, 0): value})

    def __add__(self, other):
        out = dict(self.terms)
        for key, value in self.cast(other).terms.items():
            out[key] = out.get(key, 0)+value
        return Series(out)

    __radd__ = __add__

    def __neg__(self):
        return Series({key: -value for key, value in self.terms.items()})

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return self.cast(other)+-self

    def __mul__(self, other):
        out = {}
        for (n, k), a in self.terms.items():
            for (m, j), b in self.cast(other).terms.items():
                key = n+m, k+j
                if key[0] <= self.ny and key[1] <= self.nk:
                    out[key] = out.get(key, 0)+a*b
        return Series(out)

    __rmul__ = __mul__

    def inverse(self):
        c = self.terms.get((0, 0), 0)
        if not c:
            raise ValueError('A nonzero constant is required for a formal inverse')
        z, term, total = -(self-c)/c, Series.cast(1), Series.cast(1)
        for _ in range(self.ny+self.nk):
            term = term*z
            total += term
        return total/c

    def __truediv__(self, other):
        if not isinstance(other, Series):
            return Series({key: value/Q(other) for key, value in self.terms.items()})
        return self*other.inverse()

    def __pow__(self, degree):
        if type(degree) is not int or degree < 0:
            raise ValueError('Nonnegative integer series power required')
        out = Series.cast(1)
        for _ in range(degree):
            out = out*self
        return out

    def eta(self):
        return Series({(n, k-1): k*v for (n, k), v in self.terms.items() if k})

    def radial(self):
        return Series({(n-1, k): n*v for (n, k), v in self.terms.items() if n})

    def euler(self):
        return Series({(n, k): n*v for (n, k), v in self.terms.items() if n})

    def average(self):
        return Series({(n, k): v/(n+1) for (n, k), v in self.terms.items()})

    def integral(self):
        return Series({(n+1, k): v/(n+1) for (n, k), v in self.terms.items()})

    def J(self, nu):
        if nu not in (1, 2):
            raise ValueError('Only regular radial inverses J1,J2 are used')
        return Series({(n+1, k): v/((n+1)*(n+nu)) for (n, k), v in self.terms.items()})

    def shift(self, degree):
        if degree < 0 and any(n < -degree for n, _ in self.terms):
            raise ValueError('Cannot silently discard a singular integration constant')
        return Series({(n+degree, k): v for (n, k), v in self.terms.items()})

    def limit(self, ny, nk=None):
        return Series({(n, k): v for (n, k), v in self.terms.items()
                       if n <= ny and (nk is None or k <= nk)})

    def zero(self, ny, nk):
        return not self.limit(ny, nk).terms


def finite_angular_inverse(value, chi):
    """Exact inverse on finite radial jets; convergence needs the norm proof."""
    term, total = value, value
    for _ in range(Series.ny):
        term = -(chi*term).J(2)/2
        total += term
    return total


def formal_core(degree=3, output_eta=2, pressure_mode='full'):
    """Solve finite triangular coefficient equations with rational arithmetic.

    Diagnostic h=1/200, P=16, Lambda=256 and small toy j,sigma are not the
    selected giant schedule. Neither their finite jet nor their g amplitude
    is offered as the analytic solution certified by the norm argument.
    """
    if degree < 1 or output_eta < 0 or pressure_mode not in ('full', 'zero', 'drop_eta', 'drop_radial'):
        raise ValueError('Invalid formal diagnostic request')
    Series.ny, Series.nk = degree+3, output_eta+degree+4
    one, Y, e = Series.cast(1), Series({(1, 0): 1}), Series({(0, 1): 1})
    h, lam, j, sig = Q(1, 200), Q(256), Q(1, 100), Q(1, 5)
    t, A, D = 1/lam, Q(1, 2)+h, Q(1, 2)-h
    d, L, Ustar = 1-e*e, 1-2*h*e*e, 4*e+j
    invL = L.inverse()
    Hstar = D*e+d*Ustar
    Wstar = 1-d*Ustar.eta()-2*D*e*Ustar
    den = (Hstar*Hstar+sig*sig).inverse()
    chi, zeta = Hstar*Hstar*den, -L*Hstar*den
    Pi0 = -Q(5, 2)*16**2*((1+e*e)**2).inverse()
    Z = -A*(1-2*e*Ustar)*Ustar-Hstar*Ustar.eta()-d*Pi0.eta()+4*A*e*Pi0
    # g_eta=Lambda*zeta*g enforces the azimuthal datum, not an unrelated g.
    g = {(0, 0): Q(1, 100)}
    for k in range(Series.nk):
        g[0, k+1] = lam*sum(zeta.terms.get((0, j), 0)*g.get((0, k-j), 0) for j in range(k+1))/(k+1)
    g = Series(g)

    def inverse_angular(value):
        return finite_angular_inverse(value, chi)

    phi0, u0 = inverse_angular(one), -(Y*Z*invL)/2

    def sources(phi, u):
        U = Ustar+t*u
        W = Wstar+t*(-2*D*e*u.average()-d*u.average().eta())
        H = Hstar+t*d*u
        p = (g*g*phi*phi).integral()
        R1 = ((W+h*(1-2*e*U)+d*u*zeta)*phi+W*phi.euler()+H*phi.eta())*invL
        pressure = -4*A*e*p+d*p.eta()-2*e*p.euler()
        if pressure_mode == 'zero':
            pressure = Series()
        elif pressure_mode == 'drop_eta':
            pressure -= d*p.eta()
        elif pressure_mode == 'drop_radial':
            pressure += 2*e*p.euler()
        R2 = ((A*(1-4*e*Ustar)+d*Ustar.eta())*u-2*A*e*t*u*u
              +W*u.euler()+Hstar*u.eta()+d*t*u*u.eta()+pressure)*invL
        return R1, R2, p, U, W, H

    phi, u = one, Series()
    for _ in range(degree+1):
        R1, R2, *_ = sources(phi, u)
        phi = (phi0+t*inverse_angular(R1.J(2))/2).limit(degree)
        u = (u0+t*R2.J(1)/2).limit(degree)
    R1, R2, p, U, W, H = sources(phi, u)
    Pi = Pi0+t*p
    # Original sources are computed separately and always keep full pressure.
    Sq_phi = -W*(phi+phi.euler())-h*(1-2*e*U)*phi-H*(lam*zeta*phi+phi.eta())
    Sn = -W*U.euler()-A*(1-2*e*U)*U-H*U.eta()-d*Pi.eta()+4*A*e*Pi+2*e*Pi.euler()
    angular = -2*L*lam*(Y*phi.radial().radial()+2*phi.radial())-Sq_phi
    axial = -2*L*(Y*u.radial().radial()+u.radial())-Sn
    pressure_residual = p.radial()-g*g*phi*phi
    # Forward primitives, with zero constants, independently check stress.
    F, Sq = g*phi, Sq_phi/phi
    Qs = (Y*F*Sq).integral().shift(-2)/F
    Ns = Sn.integral().shift(-1)
    angular_stress = Y*Qs*invL/lam+2*phi.euler()/phi
    axial_stress_times_E = Y*Ns*invL/lam+2*Y*u.radial()/lam
    return dict(angular=angular, axial=axial, pressure=pressure_residual,
                angular_stress=angular_stress, axial_stress_times_E=axial_stress_times_E,
                phi=phi, u=u, p=p, g=g, zeta=zeta,
                degree=degree, output_eta=output_eta,
                compared_radial_degree=degree-1, eta_padding=Series.nk-output_eta,
                diagnostic_is_actual_analytic_core=False)
