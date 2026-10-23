"""Full cylindrical/graph curl coefficients and finite radial-chart transport."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('curl_operator_columns',HERE.parent/'core_canonical_haar_column_jet_closure/bounds.py')
columns=importlib.util.module_from_spec(sp);sp.loader.exec_module(columns)
entrance=columns.entrance
P,exact,order,product=entrance.P,entrance.exact,entrance.order,entrance.product

class CurlBudget:
    """Bounds for eta curl_*(h cross A), eta=epsilon^kappa, h=n/|n|^2.

    The retained native-scale identities give eta<=1, eta*epsilon<=1,
    eta*M<=1. No graph or curl coefficient callback is supplied.
    """
    def __init__(self,homogeneous,Rmin,Rmax,dimension=5):
        self.homogeneous=homogeneous;self.Rmin,self.Rmax=map(exact,(Rmin,Rmax))
        if not 0<self.Rmin<=self.Rmax:raise ValueError('Same closed positive radial hull required')
        if dimension!=5 or type(dimension) is not int:raise ValueError('Original five coefficient variables required')
        if not callable(getattr(homogeneous,'lift',None)):raise TypeError('Full-normal reciprocal/lift recipe required')
        self.dimension=dimension
        self.power_value=max(Q(1),self.Rmax*self.Rmax,1/self.Rmin)
        self.Ufloor=min(Q(1),self.Rmin**3)

    def lift(self,j):
        order(j);v=self.homogeneous.lift(j)
        if not isinstance(v,P):raise TypeError('Positive full-normal lift jet polynomials required')
        return v

    def radial_power(self,j):
        """Raw R jets of d_r R^(d_r-1), using 1/3<d_r<5/2."""
        order(j);return P(3*self.power_value*factorial(j+1)/self.Rmin**j)

    def reciprocal_radius(self,j):
        order(j);return P(Q(factorial(j),1)/self.Rmin**(j+1))

    def directions_R(self,j):return P(2*int(j==0))+2*self.radial_power(j)
    def directions_U(self,j):return P(int(j==0))+3*self.radial_power(j)

    @lru_cache(None)
    def raw(self,j,conjugated=False):
        """SUM of induced infinity norms of D0 and every first-order matrix.

        conjugated=True first replaces D_R by (d_r R^(d_r-1))*D_R,
        as required for partial_R=(d_r R^(d_r-1))*partial_U.
        D0 is unchanged before coefficient composition with R(U).
        """
        order(j)
        V=self.directions_U if conjugated else self.directions_R
        return 2*(product(V,self.lift,j)+
            product(self.directions_R,lambda i:self.lift(i+1),j)+
            product(self.reciprocal_radius,self.lift,j))

    def inverse_map(self,j):
        """Raw derivative of R(U)=U^(1/d_r), j>=1; no analytic profile needed."""
        order(j)
        if not j:raise ValueError('Composition uses increments of positive order')
        return Q(self.Rmax*factorial(j+2),2)/self.Ufloor**j

    def increment(self,j):
        order(j)
        return Q(0) if not j else Q(int(j==1))+self.inverse_map(j)/factorial(j)

    @lru_cache(None)
    def increment_power(self,j,k):
        order(j);order(k)
        if not k:return Q(int(j==0))
        if j<k:return Q(0)
        return sum((self.increment(i)*self.increment_power(j-i,k-1) for i in range(1,j+1)),Q(0))

    def compose(self,normalized,j):
        """Normalized finite Taylor jets under R(U); each outer order <=j."""
        order(j)
        return sum((self.increment_power(j,k)*normalized(k) for k in range(j+1)),P())

    def normalized(self,j,chart='R'):
        order(j)
        if chart not in ('R','U'):raise ValueError('Choose the original R chart or its explicit U=R^d_r re-expression')
        raw=lambda k:Q(5**k,factorial(k))*self.raw(k,conjugated=chart=='U')
        return raw(j) if chart=='R' else self.compose(raw,j)

    def transport(self,normalized,chart='R'):
        """Also use for separately supplied signed-target/source jet inputs."""
        if chart not in ('R','U'):raise ValueError('Consistent R or U coefficient chart required')
        return normalized if chart=='R' else lambda j:self.compose(normalized,j)

def join_root(template,operator,chart='R'):
    """Replace the independent D_j slot, consistently transporting ALL jets for U."""
    roots=entrance.inverse.roots
    if not isinstance(template,roots.RootBudget) or not isinstance(operator,CurlBudget):raise TypeError('Original root recipe and complete curl budget required')
    if chart not in ('R','U'):raise ValueError('Consistent R or U coefficient chart required')
    convert=lambda f:operator.transport(f,chart)
    return roots.RootBudget(H=convert(template.H),target=convert(template.target),wave=convert(template.wave),
        curl=lambda j:operator.normalized(j,chart),loghalf=convert(template.loghalf),
        K0=template.K0,beta=template.beta,R=template.R,mass=template.mass,width=template.width,trace=template.trace)

def from_original(p,selector_inputs,ell,*,F,G,btilde,path,tilt,template,cutoff,chart='R'):
    p=entrance.reference.normal.frame_inputs(p)
    phase,column,c,tc,g=columns.from_original(p,selector_inputs,ell,F=F,G=G,btilde=btilde,path=path,tilt=tilt)
    operator=CurlBudget(column.homogeneous,p['Rmin'],p['Rmax'])
    base=columns.root_join(template,column,cutoff)
    return phase,column,operator,join_root(base,operator,chart),c,tc,g

def certificate():
    return dict(parent_commit='95da5280b58140f76fd390b07e27940d6a455d5a',
        result='conditional complete cylindrical/graph curl coefficient recipe and explicit radial-chart conjugation',
        exact_curl='B=er cross Dr(h cross A)+ez cross epsilon partial_Z(h cross A)+ez*(h cross A)_theta/R',
        coefficients='Dmu=cross(a_mu)cross(h); D0=sum_mu cross(a_mu)cross(partial_mu h)+(eta/R)ez etheta^T cross(h)',
        raw='2[(Vchart*h)_j+(V_R*h_shift)_j+(r*h)_j]; V_R=2 delta_j0+2 g_j, V_U=delta_j0+3 g_j',
        radial_chart='first conjugate D_R by d_r R^(d_r-1), then compose every coefficient and root input jet with R(U); no chart jets reused silently',
        input_order='D_j uses normal through j+1; full phase normal uses background through j+2; whole axial-moment chain n needs full background through n+13 via this route',
        full_pressure_memory_and_all_covariance_mean_residual_obligations_retained=True,
        actual_background_suprema_weighted_targets_full_residual_or_blowup_certified=False)
