"""Quantitative multiharmonic regular ledger; coarse mean feedback blocks gain."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


fresh=load('multi_residual_fresh',HERE.parent/'core_reweighted_fresh_mean_closure/budgets.py')
join,prefix=fresh.join,fresh.prefix
P,order,exact=prefix.P,prefix.order,prefix.exact
convert,hull=fresh.convert,fresh.monomial_upper


def harmonic_constants(harmonics):
    h=tuple(harmonics)
    if not h or any(type(m) is not int or m==0 for m in h) or len(set(h))!=len(h):
        raise ValueError('A nonempty finite set of distinct nonzero integer harmonics required')
    if set(h)!={-m for m in h}:raise ValueError('Retain conjugate harmonics of a real wave')
    H=max(map(abs,h));ratio=max(Q(abs(n),abs(m)) for m in h for n in h)
    return dict(harmonics=sorted(h),max_abs=H,max_input_ratio=ratio,
        quadratic_outputs=sorted({m+n for m in h for n in h}),
        nonzero_outputs=sorted((set(h)|{m+n for m in h for n in h})-{0}))


def increment_recipe(n,data,harmonics=(-2,-1,1,2)):
    """Inputs through n+3; T,R,Pi are UNCUT sourced increment bounds.

    Full source-normal pressure is retained in Pi. oldU,newU are complete
    CUT velocity bounds. dM bounds the actual two mean increments, not a
    subtraction of unrelated norm bounds. The cutoff slot bounds only the
    retained temporal pulse cutoff; spatial localizers remain in T,R,Pi.
    """
    order(n);hc=harmonic_constants(harmonics);H=hc['max_abs'];ratio=hc['max_input_ratio']
    keys=('b','G','graph','c','r')+prefix.LINEAR+('T','R','Pi','oldU','newU','oldM','dM')
    d=prefix.require(data,keys);D=Q(n+4)
    b,G,g,c,r=(d[k] for k in ('b','G','graph','c','r'))
    e,bR,Bz,N,K,psi=(d[k] for k in prefix.LINEAR)
    T,R,pi,old,new,M,dm=(d[k] for k in ('T','R','Pi','oldU','newU','oldM','dM'))
    dr=g*D;A=dr+D+r
    groups=dict(slow=D+b*dr+G*D,phase=2*H*e,base=bR+Bz+b*r,
        viscosity=D*D*(g*g+1)+r*dr+r*r,
        mixed=4*H*N*(dr+D)+2*H*N*(dr+D+r),connection=4*H*N*r)
    E=sum(groups.values(),P() if isinstance(T,P) else Q(0))
    pressure=psi*(dr+D)*pi
    principal_curl=psi*(c*D+K+4*H*H*N*N)*R
    linear=psi*E*(T+R)+pressure+principal_curl
    Cww=(1+ratio)*A
    Cmw=2*A+2*H*N
    both_wave_crosses=2*Cww*old*new
    square=Cww*new*new
    old_mean_cross=Cmw*M*new
    new_mean_feedback=Cmw*dm*(old+new)
    wave_stage=linear+both_wave_crosses+square+old_mean_cross
    return dict(F=wave_stage+new_mean_feedback,wave_stage=wave_stage,mean_feedback=new_mean_feedback,
        linear=linear,pressure_gradient=pressure,principal_on_uncut_curl=principal_curl,
        both_old_new_wave_crosses=both_wave_crosses,new_wave_square=square,
        old_mean_new_wave=old_mean_cross,groups=groups,Cww=Cww,Cmw=Cmw,input_order=n+3)


def rates(kappa=Q(1,100000),fast_order=None,slow_order=None):
    k=exact(kappa)
    if not 0<k<Q(1,4):raise ValueError('Original positive small kappa required')
    B=1-k;C=B+Q(1,2);D=C-k
    fast=B if fast_order is None else exact(fast_order)
    slow=B if slow_order is None else exact(slow_order)
    groups=dict(complete_increment_linear=D,both_old_new_wave_crosses=D,
        new_wave_square=2*B-k,old_mean_new_wave=B+B-Q(1,2),
        fresh_fast_mean_current_wave=fast,fresh_slow_mean_current_wave=slow)
    return dict(B=B,C=C,desired=D,groups=groups,guaranteed=min(groups.values()),
        sufficient_fast_margin=fast-D,sufficient_slow_margin=slow-D,
        improvement_from_coarse_size_only=False,
        note='Stronger mean orders are candidate sufficient hypotheses, not certified actual map bounds')


def source_pole(p):
    p=exact(p)
    if p<0:raise ValueError('Finite nonnegative original pole required')
    return dict(old_wave=p,new_wave=3*p,old_means=2*p,fresh_mean_increment=6*p,
        wave_cross=4*p,wave_square=6*p,old_mean_new_wave=5*p,mean_feedback=9*p,source=9*p)


class NextResidualBudget:
    """SAME original prefix and generated inverse; no new input estimate slot.

    All polynomials use the parent's estimate argument S+K, not a changed
    physical band. This regular source has only the old coarse order B.
    Exact additive tails stay outside every source/inverse callback.
    """
    def __init__(self,means):
        self.means=means;self.generated=means.generated;self.before=means.before
        self.harmonics=harmonic_constants(self.generated.harmonics)
        self.guaranteed_order=rates()['B']
        self.pole_ledger=source_pole(Q(means.reweighting['original_pole']))
        if self.harmonics['harmonics']!=[-2,-1,1,2]:raise ValueError('This composed stage is the first sourced +/-1,+/-2 wave')

    @lru_cache(None)
    def uncut_increment(self,n):
        order(n);out={k:P() for k in ('T','R','Pi')}
        for m in self.generated.harmonics:
            f=self.generated.modes[abs(m)]
            A=join.through(f.potential,n);An=join.through(f.potential,n+1)
            fixed=self.before.fixed(n)
            c=join.curl_budget(n,join.through(f.velocity,n),A,An,join.convert(fixed['graph']),join.convert(fixed['r']))
            for key in ('T','R'):out[key]+=convert(c[key])
            out['Pi']+=convert(join.through(f.pressure,n))
        return {k:hull(v) for k,v in out.items()}

    @lru_cache(None)
    def regular_source(self,n):
        order(n);N=n+3;b=self.before;g=self.generated;f=self.means
        fixed=b.fixed(N);linear=dict(b.linear(N))
        # Raw cutoff jets -> complete factorial through-order norm on the
        # same five coefficient variables, before the product rule.
        linear['cutoff']=hull(convert(join.through(g.cutoff_raw,N)))
        oldM=hull(sum(b.final(N).values(),P()))
        fast=f.fast(N);slow=f.slow_increment(N)
        dm=hull(sum(fast.values(),P())+sum(slow.values(),P()))
        data={k:fixed[k] for k in ('b','G','graph','c','r')}
        data.update(linear);data.update(self.uncut_increment(N))
        data.update(oldU=b.wave(N)['U'],newU=hull(convert(g.addition(N)['U'])),oldM=oldM,dM=dm)
        result=increment_recipe(n,data,self.generated.harmonics)
        # This call retains the current angular-zero pressure and rows as
        # well. Their size is not an improved zero-mode residual estimate.
        result['full_current_mean_rows']=f.final_rows(n)
        result['order_certificate']=dict(guaranteed=self.guaranteed_order,desired=rates()['desired'],
            improved_order_certified=False,source_pole=self.pole_ledger['source'])
        return result

    def raw_source(self,n):
        order(n);return factorial(n)*self.regular_source(n)['F']

    def inverse_source(self,module):return lambda n:module.P(self.raw_source(n).c)
