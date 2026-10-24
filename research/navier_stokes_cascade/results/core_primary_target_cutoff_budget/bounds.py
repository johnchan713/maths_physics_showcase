"""Original leading-target, flat-log and fixed-cutoff norm recipes."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb,factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
    sp=importlib.util.spec_from_file_location(name,HERE.parent/path)
    m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
curl=load('target_cutoff_curl','core_complete_curl_operator_budget/bounds.py')
P,exact,order=curl.P,curl.exact,curl.order
RADIUS=Q(1,6057072)
def conv(a,b,n):return sum((a(i)*b(n-i) for i in range(n+1)),P())
def seminorm(callback,j):
    order(j);v=exact(callback(j))
    if v<0:raise ValueError('Nonnegative fixed-profile normalized seminorm bound required')
    return v

class ChartBudget:
    """Same physical branch 1/4<=s<=4, |eta|<=1 and R>=Rmin>0."""
    def __init__(self,Rmin,width):
        self.Rmin,self.width=map(exact,(Rmin,width))
        if self.Rmin<=0 or not 0<self.width<=1:raise ValueError('Original positive radial hull and unchanged flat width in (0,1] required')

    def log_s(self,j):
        order(j);return Q(0) if not j else sum((Q(32**i*64**(j-i),3) for i in range(j+1)),Q(0))

    @lru_cache(None)
    def exponential(self,kind,j):
        order(j)
        if kind not in ('eta','scale'):raise ValueError('Original s^-D or s^-(1+h) factor required')
        value,beta=(Q(2),Q(1,2)) if kind=='eta' else (Q(16),Q(2))
        if not j:return value
        return sum((Q(i,j)*beta*self.log_s(i)*self.exponential(kind,j-i) for i in range(1,j+1)),Q(0))

    def increment(self,kind,j):
        order(j)
        if not j:return Q(0)
        if kind=='log_s':return self.log_s(j)
        if kind=='log2_s':return 2*self.log_s(j) # log 2 >1/2
        if kind!='profile':raise ValueError('Original profile/log-scale coordinate map required')
        xjet=Q(2,j)/self.Rmin**j+self.log_s(j)
        etajet=2*self.exponential('eta',j)+self.exponential('eta',j-1)
        return xjet+etajet

    @lru_cache(None)
    def power(self,kind,n,k):
        order(n);order(k)
        if not k:return Q(int(n==0))
        if n<k:return Q(0)
        return sum((self.increment(kind,j)*self.power(kind,n-j,k-1) for j in range(1,n+1)),Q(0))

class TargetBudget:
    """profile(j) bounds J_j(F), F=yb^3*T0/zeta, on the closed native shell.

    These are fixed smooth profile seminorms from the existing edge
    factorizations; no current weighted target callback is supplied.
    """
    def __init__(self,chart,profile):
        if not isinstance(chart,ChartBudget) or not callable(profile):raise TypeError('Original chart and fixed factored-profile seminorms required')
        self.chart,self.profile=chart,profile

    def quotient(self,j):
        order(j)
        return sum((seminorm(self.profile,i)*comb(j-i+2,2) for i in range(j+1)),Q(0))

    def composed(self,j):
        return sum((self.quotient(k)*self.chart.power('profile',j,k) for k in range(j+1)),Q(0))

    def target(self,j):
        """J_j(s^-(1+h)*T0/zeta) <= this * delta^(-j-3)."""
        order(j)
        return P(sum((self.chart.exponential('scale',i)*self.composed(j-i) for i in range(j+1)),Q(0)))

    def loghalf(self,j):
        """Positive-order jets of log sqrt(zeta), with root pole delta^-3j."""
        order(j)
        if not j:return P() # unused by RootBudget.E, NOT a bound for log sqrt(zeta) itself
        return P(sum(((self.chart.width**2+4)*Q(k+1,2)*self.chart.power('profile',j,k)
            for k in range(1,j+1)),Q(0)))

class CutoffBudget:
    """Fixed normalized ONE-DIMENSIONAL cutoff seminorms, not new profiles."""
    def __init__(self,chart,*,dyadic,mesh,transverse,temporal,radius=RADIUS):
        self.chart=chart;self.dyadic,self.mesh,self.transverse,self.temporal=dyadic,mesh,transverse,temporal
        self.radius=exact(radius)
        if not isinstance(chart,ChartBudget) or self.radius<=0 or not all(callable(f) for f in (dyadic,mesh,transverse,temporal)):
            raise ValueError('Same chart, positive radius and all four fixed profile seminorm callbacks required')

    def dyadic_jets(self,j):
        return P(sum((seminorm(self.dyadic,k)*self.chart.power('log2_s',j,k) for k in range(j+1)),Q(0)))

    def mesh_jets(self,j):
        order(j)
        value=sum((seminorm(self.mesh,a)*seminorm(self.mesh,b)*seminorm(self.mesh,j-a-b)
            for a in range(j+1) for b in range(j-a+1)),Q(0))
        return P({3*j:value}) # original mesh width S^-3, not a replacement width

    def auxiliary_jets(self,j):
        order(j)
        return P(sum((seminorm(self.transverse,i)*(2/self.radius)**i*
            seminorm(self.temporal,j-i)*(1/self.radius)**(j-i) for i in range(j+1)),Q(0)))

    def normalized(self,j,include_slow=True):
        order(j)
        if type(include_slow) is not bool:raise TypeError('Choose explicitly whether the external slow partition is included once')
        if not include_slow:return self.auxiliary_jets(j)
        return conv(lambda i:conv(self.dyadic_jets,self.mesh_jets,i),self.auxiliary_jets,j)

    def raw(self,j):return factorial(j)*self.normalized(j)

def from_original(p,selector_inputs,ell,*,F,G,btilde,path,tilt,template,profile,dyadic,mesh,transverse,temporal,chart='R'):
    """Join all primary norm slots on the same original R chart, then optional U transport."""
    p=curl.entrance.reference.normal.frame_inputs(p)
    coordinates=ChartBudget(p['Rmin'],template.width);target=TargetBudget(coordinates,profile)
    cut=CutoffBudget(coordinates,dyadic=dyadic,mesh=mesh,transverse=transverse,temporal=temporal)
    roots=curl.entrance.inverse.roots
    base=roots.RootBudget(H=template.H,target=target.target,wave=template.wave,curl=template.curl,
        loghalf=target.loghalf,K0=template.K0,beta=template.beta,R=template.R,mass=template.mass,
        width=template.width,trace=template.trace)
    bundle=curl.from_original(p,selector_inputs,ell,F=F,G=G,btilde=btilde,path=path,tilt=tilt,
        template=base,cutoff=cut.raw,chart=chart)
    return (*bundle,target,cut)

def certificate():
    return dict(parent_commit='dd73978319d6f96c410964552b1fdddfa8e1437a',
        result='conditional fixed-profile leading-target, flat-log and complete-cutoff norm input closure',
        factor='F=yb^3 T0/zeta is smooth on the whole closed native shell by the exact inner/outer factorizations; T0/zeta itself retains its cubic outer pole',
        target='J_n(s^-(1+h) T0/zeta)<=T_n delta^(-n-3), which supplies the original delta^(-3n-3) root slot without changing its map',
        loghalf='positive-order jets only; exact width^2 and outer constant4 retained, all nonlinear profile-chart derivatives included',
        cutoff='fixed profile seminorms, mesh S^(3n), dyadic log2(s), dual auxiliary covectors, ci L=2r0; slow partition included once in wave norm and remains outside H',
        primitive_profile_validity_value_cone_and_full_pressure_mean_residual_obligations_retained=True,
        numerical_actual_profile_seminorms_complete_residual_or_blowup_certified=False)
