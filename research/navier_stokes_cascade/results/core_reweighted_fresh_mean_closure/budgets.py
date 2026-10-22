"""Explicit higher-pole reweighting and fresh post-sourced mean SIZE bounds."""
from functools import lru_cache
from fractions import Fraction as Q
from math import ceil,comb,factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


join=load('fresh_reweight_join',HERE.parent/'core_sourced_covariance_moment_join/budgets.py')
prefix=join.prefix
P,order,exact=prefix.P,prefix.order,prefix.exact
X=P({1:1})


def convert(p):return P(p.c)


def monomial_upper(poly):
    """For x>=1, sum c_j*x^j <= (sum c_j)*x^max(j), exactly."""
    if not isinstance(poly,P):raise TypeError('Positive estimate polynomial required')
    return P({max(poly.c):sum(poly.c.values(),Q(0))}) if poly.c else P()


def cut_recipe(pole,L,width,gap):
    """SAME two-edge proof, with the actual cutoff plateau gap retained."""
    pole,L,width,gap=map(exact,(pole,L,width,gap))
    if pole<0 or L<=0 or not 0<width<=1 or not 0<gap<=L/2:
        raise ValueError('Nonnegative pole and original positive shell/width/gap required')
    recipe=join.edge.primitive_recipe(pole,L,width)
    r=min(Q(recipe['collar']),gap);q=recipe['kernel_degree'];a,b=width**2,Q(4)
    D=lambda c:sum((Q(comb(q,j)*factorial(j))*(r*r/c)**j for j in range(q+1)),Q(0))/(2*c)
    recipe['collar']=str(r)
    recipe['branches']=dict(left=(D(a),b/(L-r)**2),right=(D(b),a/(L-r)**2),
                            middle=(L*join.edge.atlas.edge_constant(pole,width),(a+b)/(r*r)))
    return recipe


def shift_recipe(p,L,width,gap):
    """K=3+E6+E12+sum(branch coefficient * 3^ceil(exponent)).

    Exact finite POWER recipe only: the giant selected integer need not be
    allocated to establish polynomial-in-S upper bounds.
    """
    p=exact(p)
    if p<0:raise ValueError('Nonnegative original pole required')
    recipes={pole:cut_recipe(pole,L,width,gap) for pole in (6*p,12*p)}
    terms=[]
    # Keep both calls even when p=0; a larger finite constant is harmless.
    for pole in (6*p,12*p):
        for name,(coefficient,exponent) in recipes[pole]['branches'].items():
            terms.append(dict(pole=str(pole),branch=name,coefficient=str(coefficient),power=ceil(exponent)))
    E6=join.edge.atlas.edge_constant(6*p,width);E12=join.edge.atlas.edge_constant(12*p,width)
    return dict(original_pole=str(p),L=str(exact(L)),width=str(exact(width)),gap=str(exact(gap)),
                constant=str(3+E6+E12),terms=terms,
                covered_operator_poles={k:list(map(str,v)) for k,v in join.next_pole_ledger(p)['required_operator_poles'].items()},
                bound_argument='x=S+K, an estimate argument only; no physical scale or band changes')


def materialize_shift(recipe,max_power=10000):
    """Moderate diagnostics only; the exponent guard is not a proof gate."""
    order(max_power)
    if any(t['power']>max_power for t in recipe['terms']):raise ValueError('Retain the giant fixed shift as its exact power recipe')
    return Q(recipe['constant'])+sum((Q(t['coefficient'])*3**t['power'] for t in recipe['terms']),Q(0))


def translate(poly,K):
    """Exact P(S+K), without changing the original physical S."""
    if not isinstance(poly,P):raise TypeError('Positive estimate polynomial required')
    K=exact(K)
    if K<0:raise ValueError('Nonnegative estimate shift required')
    out=P()
    for n,c in poly.c.items():
        out+=P({j:c*comb(n,j)*K**(n-j) for j in range(n+1)})
    return out


class FreshMeanBudget(prefix.PrefixBudget):
    """Old fields + complete sourced wave, then actual fresh fast/slow SIZE.

    All output polynomials are evaluated at x=S+K. Original operator
    envelopes MUST be the retained positive kernel/profile recipes, so
    they dominate the same pole-independent profile factors. A generic
    operator norm, or an arbitrary callback, does not prove that fact.
    Actual field enclosures and recipe provenance are still retained
    original input premises, not certified by this numerical constructor.
    """
    def __init__(self,generated,p,L,width,gap,trace=None):
        self.generated=generated;self.before=generated.before
        self.reweighting=shift_recipe(p,L,width,gap)
        domains={k:[Q(v) for v in vals] for k,vals in self.reweighting['covered_operator_poles'].items()}
        if not join.operator_coverage(p,domains)['pass_all']:raise ValueError('Incomplete next-pole recipe')
        before=self.before
        def fixed(n):
            old=before.fixed(n);return dict(old,Kp=old['Kp']*X,Kmom=old['Kmom']*X)
        def wave(n):return dict(U=monomial_upper(before.wave(n)['U']+convert(generated.addition(n)['U'])))
        super().__init__(fixed=fixed,initial=before.final,wave=wave,linear=before.linear,
            KT1=lambda n:before.KT1(n)*X,KA1=lambda n:before.KA1(n)*X,
            time_inverse=before.time_inverse,slow=before.slow,trace=trace)

    @lru_cache(None)
    def postfast(self,n):
        order(n)
        old=self.data('initial',n);inc=self.fast(n)
        return {k:monomial_upper(old[k]+inc[k]) for k in prefix.STATE}

    @lru_cache(None)
    def final(self,n):
        order(n)
        old=self.postfast(n);inc=self.slow_increment(n)
        return {k:monomial_upper(old[k]+inc[k]) for k in prefix.STATE}

    @lru_cache(None)
    def final_rows(self,n):
        order(n)
        wave=self.data('wave',n+3)['U']
        return prefix.full_mean(n,dict(self.data('fixed',n+3),**self.final(n+3),w=wave*wave))

    def source(self,*args,**kwargs):
        raise NotImplementedError('The next nonzero source requires the +/-1,+/-2 harmonic grammar; no residual gain is certified here')

    def raw_source(self,*args,**kwargs):return self.source(*args,**kwargs)
    def inverse_source(self,*args,**kwargs):return self.source(*args,**kwargs)
