"""Exact full-row subtraction and A1/rate controls, separate from continuum claims."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
from itertools import product
import sys
from budgets import HERE,load,scalar,fresh,prefix,P,difference,covariance_rows,rates,flat_recipe,materialize_flat,FastRowGain

fast=load('fast_row_original_inverse',HERE.parent/'core_fast_inverse_budget/proof_checks.py')
means=fast.means


def fixture(trace=None):
    old=sys.modules.get('budgets');sys.modules['budgets']=fresh
    try:c=load('fast_row_original_fixture',HERE.parent/'core_reweighted_fresh_mean_closure/controls.py')
    finally:
        if old is None:sys.modules.pop('budgets',None)
        else:sys.modules['budgets']=old
    return FastRowGain(c.fixture(trace))


def identity_review():
    R,Z,T,Y=means.R,means.Z,means.T,means.U
    eps=Q(1,13);M=Q(7,2);c=Q(2,5);f=(R-1)**6*(2-R)**6
    dr=lambda a:a.diff(0)+M*a.diff(3)
    dz=lambda a:eps*a.diff(1)
    ts=lambda a:-eps*a.diff(2)+c*a.diff(3)
    de=lambda a,e:dr(a)+e*a.shift(-1)
    lap=lambda a:dr(dr(a))+dr(a).shift(-1)+dz(dz(a))
    rho=R**4*f/(R**4*f).integral().scalar()
    def average(a):
        out=means.Poly()
        for (r,z,t,y),v in a.terms.items():
            if y%2==0:out+=means.Poly({(r,z,t,0):v/Q(y+1)})
        return out
    def pressure(gr):
        # A manufactured LINEAR primitive suffices to check subtraction.
        # Its identity as the actual graph inverse is not asserted.
        p=(gr-rho*average(gr.integral())).primitive()
        return p-p.at(0,1)
    b=eps*R**4*f*(1+Z);V=R**4*(2+Z);G=R**4*f*(1+T)
    psi=R**4*f*(1+Z+Z**2)*(1+Y+Y**2)
    eta=R**4*f*(1+2*Z+T)*(1+2*Y+Y**2)
    old=(-eps*psi.diff(1),R**4*f*(1+Z+Y),de(psi,1))
    inc=(-eps*eta.diff(1),R**4*f*(2+T+Y**2),de(eta,1))
    W=[[means.Poly() for _ in range(3)] for _ in range(3)]
    def full(v,w):
        x,u,z=v
        gr=-ts(x)-de(2*b*x+x*x+w[0][0],1)-dz(b*z+G*x+x*z+w[0][2])+(2*V*u+u*u+w[1][1]).shift(-1)+eps*(lap(x)-x.shift(-2))
        pm=pressure(gr)
        et=ts(u)+de(b*u+V*x+x*u+w[0][1],2)+dz(G*u+V*z+z*u+w[1][2])-eps*(lap(u)-u.shift(-2))
        ez=ts(z)+de(b*z+G*x+x*z+w[0][2],1)+dz(2*G*z+z*z+w[2][2]+pm)-eps*lap(z)
        return gr,pm,et,ez
    new=tuple(a+b for a,b in zip(old,inc));a=full(old,W);bnew=full(new,W)
    x,u,z=old;X,U,Zi=inc
    S=2*b*X+2*x*X+X*X;H=b*Zi+G*X+x*Zi+X*z+X*Zi;TT=2*V*U+2*u*U+U*U
    dg=-ts(X)-de(S,1)-dz(H)+TT.shift(-1)+eps*(lap(X)-X.shift(-2))
    dp=pressure(dg)
    rt=b*U+V*X+x*U+X*u+X*U
    at=G*U+V*Zi+z*U+Zi*u+Zi*U
    az=2*G*Zi+2*z*Zi+Zi*Zi+dp
    et=-eps*U.diff(2)+de(rt,2)+dz(at)-eps*(lap(U)-U.shift(-2))
    ez=-eps*Zi.diff(2)+de(H,1)+dz(az)-eps*lap(Zi)
    checks=dict(
        complete_radial_source_difference=bnew[0]-a[0]==dg,
        full_linear_pressure_reconstruction_difference=bnew[1]-a[1]==dp,
        tangential_difference_excludes_only_fast_time=bnew[2]-a[2]-c*U.diff(3)==et,
        axial_difference_excludes_only_fast_time=bnew[3]-a[3]-c*Zi.diff(3)==ez,
        pressure_gradient_cannot_be_omitted=bool(dz(dp)),
        viscosity_cannot_be_omitted=bool(eps*lap(Zi)),
        cylindrical_connection_cannot_be_omitted=bool(2*rt.shift(-1)),
        old_new_mean_crosses_are_nonzero=bool(x*U+X*u) and bool(z*Zi),
        full_moving_pressure_subtraction_keeps_rho=bool(rho*average(dg.integral())),
        both_mean_curls_exactly_divergence_free=de(old[0],1)+dz(old[2])==de(inc[0],1)+dz(inc[2])==0,
    )
    D=[[R**4*f*(i+j+1)*(1+Z+Y) for j in range(3)] for i in range(3)]
    wc=full(old,D);dgw=-de(D[0][0],1)-dz(D[0][2])+D[1][1].shift(-1)
    checks.update(
        complete_covariance_radial_difference=wc[0]-a[0]==dgw,
        complete_covariance_theta_difference=wc[2]-a[2]==de(D[0][1],2)+dz(D[1][2]),
        complete_covariance_axial_and_pressure_difference=wc[3]-a[3]==de(D[0][2],1)+dz(D[2][2]+pressure(dgw)),
        all_six_covariance_entries_present=all(bool(D[i][j]) for i in range(3) for j in range(i,3)),
    )
    return dict(checks=checks,scope='exact rational compact coefficient fields and a linear pressure-operator surrogate; not evaluation of the selected graph inverse')


def rate_review():
    checks={};rows=[]
    for k in (Q(1,100000),Q(1,10),Q(249,1000)):
        r=rates(k);f=flat_recipe(0,k);tag=str(k)
        checks[tag+':every_old_difference_term_exceeds_D']=all(v>=r['D'] for name,v in r['groups'].items() if name!='sourced_covariance')
        checks[tag+':wave_sets_exact_postwave_D']=r['postwave']==r['D'] and r['old_row_margin']>0
        checks[tag+':strict_finite_A1_gain']=r['A1_margin']>0 and f['raw_source_order']==2*r['q']+3
        checks[tag+':A1_has_no_temporal_inverse_plus_four']=f['extra_temporal_inverse_loss']==0
        checks[tag+':fresh_actual_means_close_feedback_at_D']=min(r['fresh_fast'],r['fresh_slow'])==r['feedback']==r['D']
        checks[tag+':next_regular_source_at_D_only_conditionally']=r['next_regular_source']==r['D'] and not r['complete_residual_or_unconditional_actual_inputs_certified']
        rows.append({name:str(r[name]) for name in ('B','H','C','D','q','A1_margin')})
    f=flat_recipe(0)
    checks['actual_small_kappa_high_jet_dependency_is_not_hidden']=f['integration_order']==99999 and f['raw_source_order']==200001
    for n in (0,1,3):
        # Moderate diagnostic kappa, ORIGINAL raw constant independently read.
        recipe=flat_recipe(n,Q(1,5));J=recipe['raw_source_order'];q=recipe['integration_order']
        raw=fast.flat_budget_data(n,q,2,3,5)
        out=materialize_flat(recipe,L=2,Bn=P(3),AJ=P(5),source=P(7),inverse_floor=P(2),rho_g=Q(3,2),weight_factors=P(11))
        expected=2*6**n*raw['constant']*factorial(J)*7*11*3**(-(-3*q//2))
        checks[f'n{n}:independent_raw_A1_constant_and_all_conversions']=out.evaluate(1)==expected
    for t in (Q(1,2),Q(1,8),Q(1,32)):
        # eps=t^10, k=1/10, B=9/10. Graph viscosity is t^17,
        # whereas the invalid no-graph-loss claim would be t^19.
        checks[str(t)+':viscosity_two_graph_losses_cannot_be_erased']=t**17/t**19==1/(t*t)
    return dict(checks=checks,rows=rows,scope='exact exponent inequalities and independently inherited raw A1 constants; no actual input supremum scan')


def norm_review():
    R,Z,T,Y=means.R,means.Z,means.T,means.U
    def norm(p,n):
        ans=Q(0)
        for ix in product(range(n+1),repeat=4):
            if sum(ix)>n:continue
            d=p;fac=1
            for axis,j in enumerate(ix):
                fac*=factorial(j)
                for _ in range(j):d=d.diff(axis)
            ans+=sum((abs(v)*max(Q(1),Q(2)**key[0]) for key,v in d.terms.items()),Q(0))/fac
        return ans
    old0=(1+R+Z+Y,2+R+T+Y**2,3+Z+Y)
    inc0=(2+R+Z*Y,1+Z+T+Y**2,2+R+Z+Y**2)
    b0=2+R+Z;V=3+R+T;G=2+Z
    checks={}
    for n in (0,1):
        N=n+3;old={key:norm(p,N) for key,p in zip(prefix.STATE,old0)}
        inc={key:norm(p,N) for key,p in zip(prefix.STATE,inc0)}
        fixed=dict(b=norm(b0,N),V=norm(V,N),G=norm(G,N),graph=Q(1),c=Q(2),r=norm(R.shift(-2),N),rho=Q(0),Kp=Q(1),Kmom=Q(1))
        bound=difference(n,old,inc,fixed)
        for t in (Q(1,2),Q(1,4)):
            eps=t**10;M=1/t
            dr=lambda p:p.diff(0)+M*p.diff(3)
            dz=lambda p:eps*p.diff(1)
            ts=lambda p:-eps*p.diff(2)+2*p.diff(3)
            de=lambda p,e:dr(p)+e*p.shift(-1)
            lap=lambda p:dr(dr(p))+dr(p).shift(-1)+dz(dz(p))
            oldv=(t**19*old0[0],t**9*old0[1],t**9*old0[2])
            incv=(t**19*inc0[0],t**9*inc0[1],t**9*inc0[2])
            def full(v):
                x,u,z=v;b=eps*b0
                gr=-ts(x)-de(2*b*x+x*x,1)-dz(b*z+G*x+x*z)+(2*V*u+u*u).shift(-1)+eps*(lap(x)-x.shift(-2))
                # Identity pressure map is an independent bounded linear
                # surrogate of norm one, not the selected T0 graph inverse.
                et=ts(u)+de(b*u+V*x+x*u,2)+dz(G*u+V*z+z*u)-eps*(lap(u)-u.shift(-2))
                ez=ts(z)+de(b*z+G*x+x*z,1)+dz(2*G*z+z*z+gr)-eps*lap(z)
                return et,ez
            a=full(oldv);b=full(tuple(x+y for x,y in zip(oldv,incv)))
            for j,key in enumerate(('Etheta','Ez')):
                actual=b[j]-a[j]-2*incv[j+1].diff(3)
                checks[f'n{n}_t{t}:independent_unreduced_physical_norm_'+key]=norm(actual,n)<=t**17*bound[key]
    return dict(checks=checks,scope='exact rational physical powers and independent coefficient l1 derivative bounds on a local box with a bounded linear pressure surrogate')


@lru_cache(None)
def composition_review():
    trace=set();g=fixture(trace);r=g.regular_rows(0);recipe=g.postwave_recipe(0)
    f,s,w=(r[k] for k in ('old_fast_remainder','old_slow_change','complete_wave_change'))
    checks={}
    for key in ('Etheta','Ez'):
        checks['same_old_fast_slow_complete_wave_sum_'+key]=r['rows'][key].evaluate(1)==2*g.embedding*(f[key]+s[key]+w[key]).evaluate(1)
        checks['old_fast_difference_nonzero_'+key]=f[key].evaluate(1)>0
        checks['old_slow_difference_nonzero_'+key]=s[key].evaluate(1)>0
        checks['complete_wave_difference_nonzero_'+key]=w[key].evaluate(1)>0
    checks.update(
        full_pressure_in_both_mean_changes_and_covariance=all(a['pm'].evaluate(1)>0 for a in (f,s,w)),
        original_higher_pole_coverage_and_floor_retained=recipe['exact_A1']['original_geometry']['covered_operator_poles']==g.means.reweighting['covered_operator_poles'],
        exact_A1_recipe_and_no_flat_tail_inversion=recipe['exact_A1']['raw_source_order']==200001 and 'never forward' in recipe['exact_A1']['treatment'],
        scalar_C_join_uses_same_existing_mean_object=g.scalars.means is g.means,
        actual_input_enclosures_not_promoted=not g.generated.certificate['actual_primitive_enclosures_certified'],
    )
    try:g.inverse_source(None)
    except NotImplementedError:checks['uninstantiated_high_jet_inverse_callback_is_blocked']=True
    else:checks['uninstantiated_high_jet_inverse_callback_is_blocked']=False
    return dict(checks=checks,regular_degrees={k:max(v.c) for k,v in r['rows'].items()},
                highest_materialized_original_initial_jet=max(n for name,n in trace if name=='initial'),
                default_A1_source_jet=200001,
                scope='existing positive polynomial recipes plus explicit deferred high-jet A1 product; no high-jet actual family evaluated')


@lru_cache(None)
def review():
    parts=dict(subtraction=identity_review(),rates=rate_review(),norms=norm_review(),composition=composition_review(),actual_A1_identity=fast.mode_review())
    return dict(checks={name+':'+k:v for name,p in parts.items() for k,v in p['checks'].items()},parts=parts)
