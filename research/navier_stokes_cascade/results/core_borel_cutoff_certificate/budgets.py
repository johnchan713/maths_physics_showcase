"""Euler jets of the implemented Borel cutoff and unchanged-scale certificates."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb,factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m);return m


radial=load('borel_radial_join',HERE.parent/'core_background_radial_reconstruction/budgets.py')
head=radial.head
exact=head.exact
order=head.phase.order


@lru_cache(None)
def gaussian_polynomial(j):
    """g^(j)(x)=P_j(1/x)*exp(-1/x^2), P_next=-u^2 P'+2u^3 P."""
    order(j)
    if not j:return {0:1}
    out={}
    for d,c in gaussian_polynomial(j-1).items():
        if d:out[d+1]=out.get(d+1,0)-d*c
        out[d+3]=out.get(d+3,0)+2*c
    return {d:c for d,c in out.items() if c}


def gaussian_raw(j):return sum(abs(c)*factorial((d+1)//2) for d,c in gaussian_polynomial(j).items())


@lru_cache(None)
def step_bounds(m):
    order(m);A=[gaussian_raw(j) for j in range(m+1)];B=[Q(81)]
    for j in range(1,m+1):B.append(81*sum(Q(comb(j,k))*2*A[k]*B[j-k] for k in range(1,j+1)))
    coarse=[sum(Q(comb(j,k))*A[k]*B[j-k] for k in range(j+1)) for j in range(m+1)]
    raw=[Q(1)]+coarse[1:]
    # Same squared-flat step; independent global bounds already proved in
    # cutoff_repair_derivative_budget. The all-order quotient route stays recorded.
    if m>=1:raw[1]=min(raw[1],Q(32))
    if m>=2:raw[2]=min(raw[2],Q(1216))
    return dict(gaussian=A,reciprocal=B,coarse=coarse,raw=raw)


@lru_cache(None)
def stirling(j):
    order(j)
    if not j:return (1,)
    prev=stirling(j-1)+(0,)
    return tuple((prev[k-1] if k else 0)+k*prev[k] for k in range(j+1))


def cutoff_bounds(m):
    order(m);raw=step_bounds(m)['raw'];euler=[Q(1)]
    for j in range(1,m+1):euler.append(sum(Q(stirling(j)[k])*2**k*raw[k] for k in range(1,j+1)))
    return dict(step_raw=raw,Euler=euler,K=[max(euler[:j+1]) for j in range(m+1)],
                identity='chi(s)=1-sigma(2s-1), sigma squared-flat logistic; same implemented function',
                universal_retained_schedule_verified=False)


def ceil_log2(value):
    value=exact(value)
    if value<1:raise ValueError('Pad each full normalized constant to at least one')
    p,q=value.numerator,value.denominator;k=max(0,p.bit_length()-q.bit_length())
    if p>q*2**k:k+=1
    return k


def verify_schedule(stages,previous_log2_upper=0,boundary_log2_upper=0):
    """Check PROVIDED intervals for the existing log2(a_n); never select scales.

    constants must include all normalized m=0..n, and gain_lower is a
    certified lower bound for n*h. Input interval/constant provenance remains
    a premise. A finite prefix cannot certify the infinite retained sequence.
    """
    prev,boundary=map(exact,(previous_log2_upper,boundary_log2_upper))
    if min(prev,boundary)<0 or not stages:raise ValueError('Nonempty supplied prefix and nonnegative domain/previous logs required')
    rows=[]
    for n,s in enumerate(stages,1):
        if set(s)!={'n','gain_lower','log2_lower','log2_upper','constants'} or s['n']!=n:raise ValueError('Consecutive complete fixed-scale certificate records required')
        gain,lo,hi=map(exact,(s['gain_lower'],s['log2_lower'],s['log2_upper']))
        constants=tuple(map(exact,s['constants']))
        if gain<=0 or lo<0 or hi<lo or len(constants)!=n+1 or min(constants)<1:raise ValueError('Valid fixed interval, positive gain and EVERY normalized row through n required')
        powers=[ceil_log2(c) for c in constants]
        gates=dict(separation_from_previous=lo>=prev+1,common_q_domain=lo>=boundary+1,
                   all_normalized_rows_on_entire_support=all(gain*lo>=n+c for c in powers))
        rows.append(dict(n=n,gates=gates,binary_constant_upper=powers,normalized_log_slack=[gain*lo-n-c for c in powers]))
        prev=hi
    return dict(rows=rows,all_supplied_prefix_gates=all(all(r['gates'].values()) for r in rows),
                infinite_schedule_certified=False,scales_selected_or_changed=False,
                scope='conditional verification of supplied existing binary-log intervals only; applies to all smaller q, not only an endpoint sample')


def comparison_dependence(h,Lead1,Lead2,M11,geometry):
    """Cphi(K2)=D+E*max(1,a*K2) from the SAME parent seven-term recipe."""
    h,Lead1,Lead2,M11=map(exact,(h,Lead1,Lead2,M11))
    if min(Lead1,Lead2,M11)<0:raise ValueError('Nonnegative full native tuple bounds required')
    g=dict(geometry);Rmin,Rmax=g.pop('Rmin'),g.pop('Rmax')
    def leading(j):
        if j not in (1,2):raise ValueError('Only proved low-order input supplied')
        return Lead1 if j==1 else Lead2
    def profile(n,j):
        if (n,j)!=(1,1):raise ValueError('Only proved coefficient-one first-jet input supplied')
        return M11
    b=head.BackgroundBudget(h,leading,profile,lambda j:Q(1),Rmin,Rmax,'euler')
    base=b.comparison(**g);chart=head.chart_constants(1,Rmin,Rmax)
    carriers=base['carriers'];T=exact(g['Lhigh'])*(carriers['theta']+carriers['axial'])*(chart['F']+chart['G'])
    a=2*(1+2*h)*M11;E=5*T;D=base['Cphi']-E*max(Q(1),a)
    if D<0:raise ValueError('Positive seven-term decomposition failed')
    return dict(D=D,E=E,a=a,lambda_max=exact(g['lambda_max']),tilt_upper=exact(g['U']),
                scope='all other full background, carrier and geometry bounds fixed; no accepted selector changed')


def comparison_value(dep,K):
    K=exact(K)
    D,E,a=map(exact,(dep['D'],dep['E'],dep['a']))
    if K<1 or min(D,a)<0 or E<=0:raise ValueError('Padded Euler cutoff upper norm and nonnegative recipe data required')
    return D+E*max(Q(1),a*K)


def selector_gate(dep,lambda_min,U,radius,Lreq):
    lam,U,r,L=map(exact,(lambda_min,U,radius,Lreq))
    if min(lam,r,L)<=0 or U<0:raise ValueError('Existing positive growth radius and length data required')
    if lam>dep['lambda_max'] or U<dep['tilt_upper']:raise ValueError('Selector growth and tilt data are inconsistent with the phase recipe')
    beta=lam/(4*(1+U*U)**2);K=cutoff_bounds(2)['K'][2]
    lower=comparison_value(dep,12);upper=comparison_value(dep,K)
    return dict(K2_certified_upper=K,K2_midpoint_lower=Q(12),Cphi_lower_for_this_recipe=lower,Cphi_certified_upper=upper,
                beta=beta,necessary_recipe_margin=beta*L*L-144*r*r*lower*lower,
                sufficient_recipe_margin=beta*L*L-144*r*r*upper*upper,
                necessary_recipe_gate=beta*L*L>=144*r*r*lower*lower,
                sufficient_recipe_gate=beta*L*L>=144*r*r*upper*upper,
                scope='only this normal-floor sufficient recipe, relative to its inputs; no conclusion of actual cone failure or full selector acceptance')


def fixture():
    f=radial.fixture();old=f['result'];dep=comparison_dependence(f['h'],old['leading']['through'][1],old['leading']['through'][2],old['first']['through'][1],f['geometry'])
    actual_cutoff=radial.low_comparison(f['h'],8,Q(1,2),4,(7,11,13),(5,17,19),f['state'],cutoff_bounds(2)['K'][2],f['geometry'])
    return dict(radial=f,dependence=dep,joined=actual_cutoff,
                gate=selector_gate(dep,1,2,Q(1,6057072),100000))
