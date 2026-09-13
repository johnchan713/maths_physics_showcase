"""Small exact polynomial oracle for the ORIGINAL inner source equations.

This file does not use the positive majorants. It expands the source equations
before the B.15 normalization, so missing pressure or angular terms survive as
nonzero polynomials. All coefficients are exact rational numbers.
"""
from fractions import Fraction


class Polynomial:
    def __init__(self, terms=None):
        self.terms = {tuple(k): Fraction(v) for k, v in (terms or {}).items() if v}

    @staticmethod
    def cast(value):
        return value if isinstance(value, Polynomial) else Polynomial({(): value})

    @classmethod
    def variable(cls, name):
        return cls({(name,): 1})

    def __add__(self, other):
        terms = dict(self.terms)
        for monomial, value in self.cast(other).terms.items():
            terms[monomial] = terms.get(monomial, 0)+value
        return Polynomial(terms)

    __radd__ = __add__

    def __neg__(self):
        return Polynomial({k: -v for k, v in self.terms.items()})

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return self.cast(other)+-self

    def __mul__(self, other):
        terms = {}
        for a, x in self.terms.items():
            for b, y in self.cast(other).terms.items():
                monomial = tuple(sorted(a+b))
                terms[monomial] = terms.get(monomial, 0)+x*y
        return Polynomial(terms)

    __rmul__ = __mul__

    def __bool__(self):
        return bool(self.terms)


def source_identities():
    names = ('A D eta d h t Ustar Ustar_eta Hstar Wstar zeta '
             'Phi Phi_eta DPhi u u_eta Du average_u average_u_eta '
             'Pi0 Pi0_eta p p_eta Dp')
    v = {n: Polynomial.variable(n) for n in names.split()}
    (A, D, eta, d, h, t, U0, Ue, H0, W0, zeta, phi, phie, Dphi,
     u, ue, Du, au, aue, Pi0, Pi0e, p, pe, Dp) = v.values()
    U = U0+t*u
    H = H0+t*d*u
    B = -2*D*eta*au-d*aue
    W = W0+t*B
    Pi = Pi0+t*p

    # Multiply the original angular equation by t*Phi. The order-one
    # contribution is Hstar*zeta*Phi=-L*chi*Phi.
    original_angular = t*(W*(phi+Dphi)+h*(1-2*eta*U)*phi+H*phie)+H*zeta*phi
    leading_angular = H0*zeta*phi
    LR1 = (W+h*(1-2*eta*U)+d*u*zeta)*phi+W*Dphi+H*phie

    original_axial = (A*(1-2*eta*U)*U+W*t*Du+H*(Ue+t*ue)
                      -4*A*eta*Pi+d*(Pi0e+t*pe)-2*eta*t*Dp)
    leading_axial = A*(1-2*eta*U0)*U0+H0*Ue-4*A*eta*Pi0+d*Pi0e
    LR2 = ((A*(1-4*eta*U0)+d*Ue)*u-2*A*eta*t*u*u
           +W*Du+H0*ue+d*t*u*ue-4*A*eta*p+d*pe-2*eta*Dp)
    return {
        'angular_remainder_exact': not (original_angular-leading_angular-t*LR1),
        'axial_remainder_exact': not (original_axial-leading_axial-t*LR2),
        'angular_transport_is_nonzero': bool(t*H*phie),
        'axial_pressure_derivative_is_nonzero': bool(t*d*pe),
        'radial_pressure_term_is_nonzero': bool(2*eta*t*Dp),
        'mixed_average_derivative_is_nonzero': bool(t*t*d*aue*Dphi),
    }


def stress_values(p1r, p2r, kappa, ratio, dp1=0, dp2=0):
    """Exact finite rational evaluation of (4.11), including shear signs."""
    p1r, p2r, kappa, ratio, dp1, dp2 = map(Fraction, (p1r, p2r, kappa, ratio, dp1, dp2))
    if p1r <= 0 or kappa <= 0 or ratio <= 0:
        raise ValueError('Positive radial component, shear factor and E ratio required')
    a = kappa*p1r
    bs = -kappa*p2r*ratio
    ts = -bs/a
    p1, p2 = p1r+dp1, p2r+dp2
    return dict(a=a, bs=bs, ts=ts, vs=a+bs*bs/a,
                Pc=p1+p2*ts, Jc=p2-p1*ts,
                vr=p1r+p2r*p2r/p1r)


def stress_identities():
    """Check cancellation uniformly in symbolic perturbations after clearing p1r."""
    p, q, k, e, x, y = [Polynomial.variable(n) for n in ('p', 'q', 'k', 'e', 'x', 'y')]
    ratio = 1+e
    Pc_p = (p+x)*p+(q+y)*q*ratio
    vr_p = p*p+q*q
    J_p = (q+y)*p-(p+x)*q*ratio
    vs_p = k*(p*p+q*q*ratio*ratio)
    return {
        'Pc_difference_exact': not (Pc_p-vr_p-(x*p+q*q*e+y*q*ratio)),
        'J_difference_exact': not (J_p-(y*p-p*q*e-x*q*ratio)),
        'vs_difference_exact': not (vs_p-k*vr_p-k*q*q*(2*e+e*e)),
        'zero_kappa_division_needed': True,
    }
