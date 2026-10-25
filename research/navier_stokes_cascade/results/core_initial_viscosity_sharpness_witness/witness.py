"""Local Fourier-jet witness for sharp first-real graph-viscosity rates.

This is an exact algebraic obstruction inside the inherited norm class.  It
does not assert that the actual pulse equals the witness or provide missing
continuum profile bounds.
"""
from dataclasses import dataclass
from fractions import Fraction as Q


def exact(x):
    return x if isinstance(x,Q) else Q(x)


@dataclass(frozen=True)
class Mode:
    """Coefficient of cos(phi), sin(phi), with phi=y+tau."""
    cosine: Q=Q(0)
    sine: Q=Q(0)

    def __post_init__(self):
        object.__setattr__(self,'cosine',exact(self.cosine))
        object.__setattr__(self,'sine',exact(self.sine))

    def __add__(self,other):
        if not isinstance(other,Mode):other=Mode(other,0)
        return Mode(self.cosine+other.cosine,self.sine+other.sine)

    def __neg__(self):return Mode(-self.cosine,-self.sine)
    def __sub__(self,other):return self+(-other)
    def __mul__(self,a):
        a=exact(a);return Mode(a*self.cosine,a*self.sine)
    __rmul__=__mul__

    def derivative(self):
        """Both active graph and temporal derivatives act as d/dphi."""
        return Mode(self.sine,-self.cosine)

    def inverse_derivative(self):
        if self.haar()!=0:raise ValueError('Only zero-Haar modes are invertible')
        return -self.derivative()

    def haar(self):return Q(0)
    def nonzero(self):return self.cosine!=0 or self.sine!=0


def exponents(kappa=Q(1,100000)):
    k=exact(kappa)
    if not 0<k<=Q(1,100):raise ValueError('Retained positive small-kappa regime required')
    alpha=2-3*k;D=alpha-k;nu=3-6*k
    return dict(kappa=k,alpha=alpha,D=D,nu=nu,
                old_principal=Q(1,2),signed_principal=alpha-Q(1,2))


def local_witness(kappa=Q(1,100000)):
    """Exact active-interior principal-cross -> fast-viscosity chain.

    On a collar where all cutoffs are one and the graph coefficient is
    constant, take old cosine coefficient A=eps^(1/2)e_r and signed cosine
    coefficient L=2 eps^(alpha-1/2)sin(phi)e_theta.  Then
    Lambda_rtheta=1/2(A_r L_theta+L_r A_theta)=eps^alpha sin(phi).
    We record only coefficients and epsilon exponents, so no huge powers are
    allocated.
    """
    o=exponents(kappa);k,alpha,D,nu=(o[n] for n in ('kappa','alpha','D','nu'))
    W=Mode(0,1)                         # eps^alpha Lambda_rtheta
    graph_W=W.derivative()              # eps^(alpha-k) cos(phi)
    connection_W=2*W                    # eps^alpha 2W/R at R=1
    source=graph_W                       # unique order-D coefficient
    fast=source.derivative()             # -Nt^{-1}source = N source
    time=fast.derivative()               # Nt fast = -source
    graph2_fast=fast.derivative().derivative()
    viscosity=-graph2_fast               # -eps Dr^2 fast, coefficient only
    return dict(orders=o,
      covariance=dict(mode=W,order=alpha,connection_order=alpha,
        realization='A=eps^(1/2)e_r; L=2 eps^(alpha-1/2) sin(phi)e_theta'),
      source=dict(mode=source,order=D,connection_mode=connection_W,
        connection_order=alpha),
      fast=dict(mode=fast,order=D,time_mode=time),
      viscosity=dict(mode=viscosity,order=D+1-2*k,
        one_graph_order=D+1-k,zero_graph_order=D+1),
      exact_fast_time_cancellation=time+source,
      strict_margin=D+1-2*k-nu,
      all_modes_zero_Haar=all(m.haar()==0 for m in (W,source,fast,viscosity)))


def logical_consequence(kappa=Q(1,100000)):
    """What the witness does and does not prove."""
    w=local_witness(kappa);v=w['viscosity'];o=w['orders']
    return dict(source_rate_attained=w['source']['order']==o['D'] and w['source']['mode'].nonzero(),
      radial_viscosity_rate_attained=v['order']==o['nu'] and v['mode'].nonzero(),
      leading_viscosity_unique_within_vector_laplacian=v['one_graph_order']>o['nu'] and v['zero_graph_order']>o['nu'],
      inherited_upper_bounds_force_positive_gain=False,
      actual_pulse_special_cancellation_excluded=False,
      required_actual_identity=('vanishing active principal-cross graph coefficient, graph-frequency orthogonality, '
        'or an exact same-order cancellation proved for the actual current pulse'))
