"""Current scalar gain from the existing coarse mean constructor, not new inputs."""
from fractions import Fraction as Q
from functools import lru_cache
from math import ceil
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


fresh=load('scalar_gain_fresh',HERE.parent/'core_reweighted_fresh_mean_closure/budgets.py')
prefix=fresh.prefix
P,order,exact=prefix.P,prefix.order,prefix.exact
SCALARS=('P','Jtheta','Jz')


def slow_reset(n,old,inc,b,G,M):
    """All scalar remainders AFTER the exact reserved-hull linear solve.

    Inputs through n+2: tangential order B and radial B+1, 0<B<=2;
    b has its extra epsilon and G is bounded. Output order 2B.
    G occurs in the radial-axial flux even in the more general formula.
    This does not cancel the separate G*u and 2G*z linear scalar leaks
    if the prescribed G=0 reserved-hull contract is violated.
    """
    order(n)
    if set(old)!=set(prefix.STATE) or set(inc)!=set(prefix.STATE):raise ValueError('Three complete mean components required')
    x,u,z=(inc[k] for k in prefix.STATE)
    ox,ou,oz=(old[k] for k in prefix.STATE)
    S=2*(b+ox)*x+x*x
    H=(b+ox)*z+(G+oz)*x+x*z
    T=2*ou*u+u*u
    K=oz*u+ou*z+u*z
    Zf=2*oz*z+z*z
    tau=(n+1)+(n+1)*(n+2)
    p=M*(tau*x+S+T+(n+1)*H)
    q=M*(tau*x+S+T+(n+1)*H)
    return dict(P=p,Qrem=q,Jtheta=M*K,Jz=M*Zf+q*Q(1,2))


def fast_change(n,old,inc,M):
    """Actual zero-Haar fast field, including A1; output order 2B."""
    order(n)
    if set(old)!=set(prefix.STATE) or set(inc)!=set(prefix.STATE):raise ValueError('Three complete mean components required')
    x,u,z=(inc[k] for k in prefix.STATE)
    ox,ou,oz=(old[k] for k in prefix.STATE)
    S=2*ox*x+x*x
    H=ox*z+oz*x+x*z
    T=2*ou*u+u*u
    Zf=2*oz*z+z*z
    K=oz*u+ou*z+u*z
    p=M*(S+T+(n+1)*H)
    q=M*(S+T+(n+1)*H)
    return dict(P=p,Qrem=q,Jtheta=M*K,Jz=M*Zf+q*Q(1,2))


def wave_change(n,D,Dnext,M):
    """Complete covariance change, all radial weights included in M."""
    order(n)
    return dict(P=M*(2*D+(n+1)*Dnext),Jtheta=M*D,
                Jz=M*(2*D+Q(n+1,2)*Dnext))


def rates(kappa=Q(1,100000)):
    k=exact(kappa)
    if not 0<k<Q(1,4):raise ValueError('Retained 0<kappa<1/4 required')
    B=1-k;C=B+Q(1,2);D=C-k
    return dict(B=B,C=C,D=D,old_postslow=2*B,wave_change=C,fast_change=2*B,
                postfast=min(2*B,C),fresh_slow=C,fast_feedback=B,
                slow_feedback=C,scalar_margin=2*B-C,
                complete_feedback=min(B,C),fast_row_gain_proved=False)


class ScalarGainBudget:
    """No outgoing scalar callback: all polynomials use the existing prefix.

    Conditional on the SAME exact five-map, G=0 on its entire potential
    hull, complete zero-Haar fast fields and original actual input/recipe
    contracts. This is not a verifier of those analytical premises.
    """
    def __init__(self,means):
        self.means=means;self.before=means.before;self.generated=means.generated
        p=Q(means.reweighting['original_pole']);L=Q(means.reweighting['L'])
        # Means/covariances are first put at pole 6p, products at 12p;
        # this constant also embeds lower-pole moment inputs into 12p.
        # Do not assume operator constants are monotone on a long shell.
        self.embedding=max(Q(1),L)**ceil(12*p)

    def embed(self,d):return {k:self.embedding*v for k,v in d.items()}

    def moment(self,n):return self.embedding*self.means.data('fixed',n)['Kmom']

    @lru_cache(None)
    def old_postslow(self,n):
        order(n);N=n+2;b=self.before
        return slow_reset(n,self.embed(b.postfast(N)),self.embed(b.slow_increment(N)),
                          b.data('fixed',N)['b'],b.data('fixed',N)['G'],self.moment(N))

    @lru_cache(None)
    def sourced_change(self,n):
        order(n)
        D=lambda j:self.embedding*P(self.generated.covariance(j).c)
        return wave_change(n,D(n),D(n+1),self.moment(n+1))

    @lru_cache(None)
    def fresh_fast_change(self,n):
        order(n);N=n+1
        return fast_change(n,self.embed(self.before.final(N)),self.embed(self.means.fast(N)),self.moment(N))

    @lru_cache(None)
    def postfast(self,n):
        order(n);a,b,c=self.old_postslow(n),self.sourced_change(n),self.fresh_fast_change(n)
        return {k:fresh.monomial_upper(a[k]+b[k]+c[k]) for k in SCALARS}

    @lru_cache(None)
    def slow_increment(self,n):
        order(n);a=self.postfast(n);an=self.postfast(n+1)
        k=self.means.data('slow',n);kn=self.means.data('slow',n+1)
        return dict(u=k['H']*(k['Kp']*a['P']+k['Kz']*a['Jz']),
                    z=k['H']*k['Kg']*a['Jtheta'],x=(n+1)*kn['H']*kn['Kpsi']*an['Jtheta'])

    def whole_axial_moment(self,n,c_rho):
        """c_rho is the whole through-n product norm; do not freeze its jets."""
        order(n);a=self.postfast(n)
        return a['Jz']+c_rho*a['P']

    def unsupported_gain(self,*args,**kwargs):
        raise NotImplementedError('Actual zero-Haar post-wave fast rows at D remain unproved')
