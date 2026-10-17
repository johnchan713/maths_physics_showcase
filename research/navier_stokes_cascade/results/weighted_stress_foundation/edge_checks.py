"""Independent finite exact kernel/ODE jets and labelled heat-tail controls."""
from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
from math import factorial
from pathlib import Path
import sys
import mpmath as mp
from edge_bounds import binomial,volterra_coefficients,flat_power_polynomial,cutoff_integer_log

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m


@lru_cache(None)
def maps():return load('weighted_foundation_rational_maps',HERE.parent/'modulation_repair_all_order/map_checks.py')


def truncate(v,n):return maps().jet(v.terms,v.n,n)


def volterra_jet(j,c=4,method='kernel'):
    m=maps();terms={};c=Q(c)
    for b in range(j.order+1):
        values=volterra_coefficients([j.terms.get((a,b),Q(0)) for a in range(j.order-b+1)],c,method)
        terms.update({(a,b):v for a,v in enumerate(values) if v})
    return m.jet(terms,2,j.order)


def divide_delta(v,k):
    if type(k) is not int or k<0 or k>v.order:raise ValueError('Available nonnegative division order required')
    if any(a[0]<k for a in v.terms):raise ValueError('Jet is not divisible by delta**k')
    return maps().jet({(a-k,b):c for (a,b),c in v.terms.items()},2,v.order-k)


def rising(h,j):
    out=Q(1)
    for k in range(j):out*=h+k
    return out


@lru_cache(None)
def outer_jets(center=Q(1),order=8,method='kernel',drop_beta=False):
    """Actual finite Gamma pole jets, or a labelled power-tail eta=0 control.

    b0,bz,b1 are divided by their common positive rho*e factor. The pole
    heat amplitude is normalized to one and rb=16, Xb=128; not giant data.
    """
    if center not in (Q(-1),Q(0),Q(1)):raise ValueError('Use heat poles or the labelled power control')
    m=maps();M=order+3;s=m.variable(0,order=M);eta=center+m.variable(1,order=M)
    h=Q(1,200);A=Q(1,2)+h;L=1-2*h*eta*eta;d=1-eta*eta
    r=16*m.exponential(-s/2);expA=m.exponential(A*s)
    if center:
        Z=2*d*m.exponential(s)/128
        H=sum(((-1)**j*rising(h,j)*rising(1+h,j)*Z**j/factorial(j) for j in range(M+1)),m.jet(0,2,M))
        Hp=sum(((-1)**(j+1)*rising(h,j+1)*rising(1+h,j+1)*Z**j/factorial(j) for j in range(M+1)),m.jet(0,2,M))
        K=expA*H;tH=-Z*Hp/H
    else:K=expA;tH=m.jet(0,2,M) # power fixture; actual eta=0 Gamma checked separately
    Bflat=m.exponential((1-s/2)**-2-1)*(8+s**3*(1-s/2)**-3)
    Bheat=2+2*h-2*tH;alpha=4*eta*eta/L**2
    C=alpha*(1-h+tH) if drop_beta else 2*L**-2+8*h*eta*eta*d/L**3+alpha*tH
    V=lambda j:volterra_jet(j,4,method)
    b0=2*K/r*Bflat+s**3/(8*r**2)*(V(r*K*Bheat*Bflat)+V(r**3*K*Bflat)/(2*L))
    b1=alpha*r*K*Bflat/2+s**3/(16*r**2)*V(r**3*K*C*Bflat)
    numerator=V(r**2*K**2*Bflat)-r**2*V(K**2*Bflat) # all f endpoint jets equal one
    bz=eta*divide_delta(numerator,3)/(8*L*r)
    b0,b1,bz=[truncate(v,order) for v in (b0,b1,bz)]
    sd=m.variable(0,order=order)
    norm=(b0*b0+sd**12*bz*bz)**Q(1,2)
    return dict(b0=b0,bz=bz,b1=b1,W=b1/norm,nz=sd**6*bz/norm)


@lru_cache(None)
def kernel_stats():
    rows=[];count=0
    for c in (Q(1,10000),Q(1),Q(4)):
        for degree in range(16):
            j=[Q(0)]*16;j[degree]=Q(1)
            a=volterra_coefficients(j,c,'kernel');b=volterra_coefficients(j,c,'ode')
            count+=len(a);rows.append(a==b)
    # The omitted Jacobian first changes the coefficient of delta^3 for J=delta.
    j=[Q(0),Q(1),Q(0),Q(0)];control=volterra_coefficients(j,4)
    return dict(count=count,match=all(rows),missing_jacobian_coefficient=control[3],
        mixed_coefficients_equal=all(not (volterra_jet(maps().jet({(a,b):Q((-1)**a,a+b+1) for a in range(9) for b in range(9-a)},2,8),4,'kernel')-
            volterra_jet(maps().jet({(a,b):Q((-1)**a,a+b+1) for a in range(9) for b in range(9-a)},2,8),4,'ode')).terms for _ in [0]))


@lru_cache(None)
def outer_stats():
    count=0;equal=True;endpoints=True;traces=[]
    for center in (Q(-1),Q(0),Q(1)):
        a=outer_jets(center);b=outer_jets(center,method='ode')
        equal &= all(not (a[k]-b[k]).terms for k in a)
        count+=len(a)*45
        L=1-2*Q(1,200)*center**2;expected=2*center**2*128/L**2
        endpoints &= a['W'].value()==expected and a['b0'].value()>0
        traces.append(str(a['W'].value()))
    zero=outer_jets(Q(0));wrong=outer_jets(Q(0),drop_beta=True)
    return dict(count=count,equal=equal,endpoints=endpoints,traces=traces,
        zero_first=zero['W'].partial((1,0)),zero_second=zero['W'].partial((2,0)),
        cubic=zero['W'].partial((3,0)),eta_second=zero['W'].partial((0,2)),
        omitted_beta_cubic=wrong['W'].partial((3,0)),
        axial_factors_nonzero=all(outer_jets(e)['bz'].value()!=0 for e in (Q(-1),Q(1))),
        outer_pole_heat_jets_not_analytic_claim=True)


@lru_cache(None)
def heat_stats():
    tail=load('weighted_foundation_actual_heat_tail',HERE.parent/'outer_edge_stability/tail.py')
    with mp.workdps(70):
        heat=tail.HeatTail(order=32,panels=2);rows=[]
        target=heat.xb/8
        for delta in map(mp.mpf,('.02','.002','.0002','.00002')):
            a=heat.factored(delta,mp.mpf(0));ratio=a['b1']/mp.sqrt(a['b0']**2+delta**12*a['bz']**2)
            scaled=ratio/(delta**3*target)
            rows.append(dict(delta=str(delta),scaled=float(scaled),relative_error=float(abs(scaled-1)),positive=bool(ratio>0)))
        # Independent positive Gamma integral versus Tricomi U at the heat endpoint.
        Z=mp.mpf(2)/heat.xb
        direct=mp.quad(lambda v:mp.exp(-v)*v**heat.h*(1+Z*v)**(-heat.h),[0,1,mp.inf])/mp.gamma(1+heat.h)
        gap=abs(direct-tail.heat_derivative(heat.h,Z))
        return dict(rows=rows,checks=dict(eta_zero_actual_Gamma_cubic_limit_refines=all(b['relative_error']<a['relative_error']/5 for a,b in zip(rows,rows[1:])),
            actual_Gamma_cubic_limit_agrees_with_Xb_over_eight=rows[-1]['relative_error']<1e-4,
            beta_term_gives_positive_eta_zero_interior_correction=all(r['positive'] for r in rows),
            independent_positive_Gamma_integral_matches_Tricomi_at_eta_zero=gap<mp.mpf('1e-55'),
            finite_diagnostic_tail_is_not_the_actual_giant_schedule=True))


@lru_cache(None)
def budget_stats():
    previous=0;rows=[];all_checks=[]
    h=Q(1,200)
    for n in range(2,15):
        logs=[(n+1)**2*(m+1) for m in range(n+1)]
        ell=cutoff_integer_log(h,n,logs,previous,boundary=2)
        margins=[n*h*(ell+extra)-v-n for v in logs for extra in (0,1,10,1000000)]
        all_checks.extend(x>0 for x in margins)
        rows.append(dict(n=n,ell=ell,minimum_margin=str(min(margins))))
        assert ell>=previous+1;previous=ell
    return dict(rows=rows,all_smaller_q=all(all_checks),count=len(all_checks))


@lru_cache(None)
def diagnostic_review():
    k,o,h,b=kernel_stats(),outer_stats(),heat_stats(),budget_stats()
    return dict(exact_kernel_components=k['count'],exact_factored_mixed_components=o['count'],
        exact_budget_inequalities=b['count'],outer_endpoint_traces=o['traces'],heat_rows=h['rows'],
        checks=dict(kernel_moments_match_independent_endpoint_ODE_through_order_15=k['match'],
            mixed_parameter_kernel_jets_match=k['mixed_coefficients_equal'],
            missing_radial_Jacobian_changes_delta_cubed=k['missing_jacobian_coefficient']==-Q(1,8),
            all_factored_outer_components_match_independent_kernel_ODE=o['equal'],
            eta_zero_and_both_heat_pole_endpoint_quotients_agree=o['endpoints'],
            eta_zero_first_two_radial_jets_vanish_but_third_is_3Xb_over_4=o['zero_first']==o['zero_second']==0 and o['cubic']==96,
            eta_zero_second_parameter_jet_is_4Xb=o['eta_second']==512,
            omitting_beta_erases_required_nonzero_cubic=o['omitted_beta_cubic']==0 and o['cubic']!=0,
            axial_factor_not_lost_at_heat_poles=o['axial_factors_nonzero'],
            exact_single_schedule_controls_whole_smaller_q_domain=b['all_smaller_q'],
            generic_weighted_error_can_reverse_sign=1-Q(1,100)/(Q(1,200))<0,
            fixed_flat_constants_and_order_loss=flat_power_polynomial(4,-3,1)=={3:Q(8),1:Q(-3)},
            **h['checks']))
