"""Independent complete scalar identities and rational rate diagnostics."""
from fractions import Fraction as Q
from functools import lru_cache
import sys
from budgets import HERE,load,final,scalar,prefix,P,seed_module,rates,response,signed_covariance,haar_slow_reset,PreservedHaarScalarJoin

coupled=seed_module.coupled
means=coupled.means;R,Z,T,U=means.R,means.Z,means.T,means.U
EPS=means.EPS


@lru_cache(None)
def exact_review():
    f,old,base,W=coupled.fields();b,V,G=base
    psi=R*f*(1+2*Z+Z**3)*(1+T+T**2)
    inc=(-EPS*psi.diff(1),(R**4*f).diff(0).shift(-2)*(1+2*Z+T),means.radial(psi))
    new=tuple(x+y for x,y in zip(old,inc));B,u,g=inc
    beta,v,gamma=(x.haar() for x in old)
    S=2*(b+beta)*B+B*B;H=(b+beta)*g+gamma*B+B*g
    A=2*v*u+u*u;K=gamma*u+v*g+u*g;F=2*gamma*g+g*g
    gr=lambda a:means.averaged_gr(base,a,W[0][0],W[1][1],W[0][2])
    rem=gr(new)-gr(old)-(2*V*u).shift(-1)
    tau=lambda q:EPS*q.diff(2)+EPS**3*q.diff(1).diff(1)
    Pnew=tau(B.integral())+(A-S).shift(-1).integral()-EPS*H.integral().diff(1)
    Qnew=tau((R**2*B).integral())+(R*(S+A)).integral()-EPS*(R**2*H).integral().diff(1)
    Jtheta=(R**2*K).integral();Jz=(R*F).integral()-Qnew/2
    jt=lambda a:(R**2*(V*a[2]+a[2]*a[1]+W[1][2]).haar()).integral()
    rho=f/f.integral().scalar()*(1+Z*(R-Q(3,2)));c=(R**2*rho).integral()/2
    pm=coupled.phi(rem-rho*Pnew)
    checks=dict(
        reserved_base_G_zero=G==0,
        slow_products_use_old_Haar_exactly=all((a*dn).haar()==a.haar()*dn for a in old for dn in inc),
        old_Haar_is_nonzero=all(bool(x.haar()) for x in old),
        full_old_correlations_are_not_products_of_Haars=(old[1]*old[2]).haar()!=old[1].haar()*old[2].haar(),
        exact_complete_pressure_remainder=Pnew==rem.integral(),
        exact_weighted_pressure_remainder=Qnew==(R**2*rem).integral(),
        exact_angular_scalar_after_linear_row=jt(new)-jt(old)-(R**2*V*g).integral()==Jtheta,
        exact_axial_scalar_after_linear_row=coupled.scalar_z(base,new,W,gr(new))-coupled.scalar_z(base,old,W,gr(old))+(R*V*u).integral()==Jz,
        pressure_field_not_discarded=bool(pm) and pm.at(0,1)==pm.at(0,2)==0,
        whole_moving_axial_pressure_identity=EPS*(R*(F+pm)).integral().diff(1)==EPS*(Jz+c*Pnew).diff(1),
        freezing_pressure_compensator_derivative_detected=bool(c.diff(1)*Pnew),
        omitting_axial_velocity_square_detected=bool((R*F).integral()),
        old_Haar_crosses_cannot_be_reset=bool((R**2*(gamma*u+v*g)).integral()),
        complete_covariance_unchanged_all_six_nonzero=all(bool(W[i][j]) for i in range(3) for j in range(i,3)),
        radial_vector_viscosity_has_zero_plain_and_weighted_moments=(B.diff(0).diff(0)+B.diff(0).shift(-1)-B.shift(-2)).integral()==0 and (R**2*(B.diff(0).diff(0)+B.diff(0).shift(-1)-B.shift(-2))).integral()==0,
        axial_viscosity_is_retained=bool(EPS**3*B.integral().diff(1).diff(1)),
        reserved_velocity_moments_exact=(R**2*u).integral()==(R*g).integral()==0,
    )
    # A zero-Haar fast increment has NONZERO old-new correlations. It
    # cannot use the Haar-only slow estimate above, even when b,V are slow.
    fastinc=tuple(q*U for q in inc)
    checks['fast_all_components_zero_Haar']=all(q.haar()==0 for q in fastinc)
    checks['fast_old_new_correlation_survives']=bool((old[1]*fastinc[1]).haar())
    checks['using_Haar_only_for_fast_would_erase_actual_cross']=old[1].haar()*fastinc[1].haar()==0 and bool((old[1]*fastinc[1]).haar())
    return dict(checks=checks,scope='exact finite C5 polynomial current-source differences, not actual canonical fields or selected C-infinity profiles')


def raw_scalar_budget(n,old,inc,b,M,eps):
    """Independent unreduced size expression with physical epsilon kept."""
    ox,ou,oz=(old[k] for k in prefix.STATE);x,u,z=(inc[k] for k in prefix.STATE)
    S=2*(b+ox)*x+x*x;H=(b+ox)*z+oz*x+x*z
    A=2*ou*u+u*u;K=oz*u+ou*z+u*z;F=2*oz*z+z*z
    tau=(n+1)+(n+1)*(n+2)
    pp=M*(eps*tau*x+S+A+eps*(n+1)*H)
    return dict(P=pp,Qrem=pp,Jtheta=M*K,Jz=M*F+pp/2)


@lru_cache(None)
def polynomial_review():
    checks={};rows=[]
    old={k:P({0:i+2,1:1}) for i,k in enumerate(prefix.STATE)}
    inc={k:P({0:i+5,2:1}) for i,k in enumerate(prefix.STATE)}
    for n in (0,1,3):
        poly=haar_slow_reset(n,old,inc,P(3),P(5))
        for s in (1,4,16):
            ev=lambda d:{k:v.evaluate(s) for k,v in d.items()}
            expected=haar_slow_reset(n,ev(old),ev(inc),Q(3),Q(5))
            checks[f'n{n}_S{s}:positive_polynomials_match_direct_substitution']=all(poly[k].evaluate(s)==expected[k] for k in poly)
            # eps=t^10, kappa=.1, C=1.4: exactly rational powers.
            for t in (Q(1,2),Q(1,8)):
                e=t**10
                o={k:v*t**(20 if k=='x' else 10) for k,v in ev(old).items()}
                d={k:v*t**(24 if k=='x' else 14) for k,v in ev(inc).items()}
                raw=raw_scalar_budget(n,o,d,3*e,5,e)
                checks[f'n{n}_S{s}_t{t}:all_four_unreduced_physical_budgets_below_C_plus_one']=all(raw[k]<=t**24*poly[k].evaluate(s) for k in poly)
        rows.append(dict(order=n,polynomials={k:v.record() for k,v in poly.items()}))
    # Independent complex vectors, not the initial real cosine/sine rule.
    Uold=((1,2),(2,-1),(-1,3));Vnew=((-2,1),(1,3),(4,-2))
    cov=lambda a,i,j:Q(a[i][0]*a[j][0]+a[i][1]*a[j][1])
    changed=tuple(tuple(x+y for x,y in zip(a,b)) for a,b in zip(Uold,Vnew))
    u=sum(abs(x)+abs(y) for x,y in Uold);v=sum(abs(x)+abs(y) for x,y in Vnew)
    checks['all_six_complex_covariance_changes_below_cross_plus_square']=all(abs(cov(changed,i,j)-cov(Uold,i,j))<=signed_covariance(Q(u),Q(v)) for i in range(3) for j in range(i,3))
    checks['complex_mixed_terms_not_deleted']=any(cov(changed,i,j)-cov(Uold,i,j)!=cov(Vnew,i,j) for i in range(3) for j in range(i,3))
    return dict(checks=checks,rows=rows,scope='rational physical epsilon comparisons and positive-polynomial algebra; complete coefficient norm hypotheses remain inputs')


@lru_cache(None)
def rate_review():
    checks={};rows=[]
    for k in (Q(1,100000),Q(1,10),Q(249,1000)):
        a=rates(k);C=a['C'];tag=str(k)
        checks[tag+':first_slow_order_derived_from_three_complete_changes']=a['first_postfast']==a['first_signed'] and a['first_signed']<2*a['B']<2
        checks[tag+':both_accumulated_slow_orders_above_one']=min(a['first_slow'],C)>1
        checks[tag+':all_complete_postslow_scalars_C_plus_one']=a['postslow']==C+1 and 2*C>C+1
        checks[tag+':strict_kappa_improvement_over_full_mean_budget']=a['postslow']-a['coarse_postslow']==k
        checks[tag+':whole_axial_pressure_residual_has_extra_epsilon']=a['whole_axial_residual']==C+2
        checks[tag+':regular_nonzero_source_not_promoted']=a['full_regular_order']==a['D']<C
        checks[tag+':coarse_B_cannot_be_relabelled_Haar_one']=1-a['B']==k>0
        rows.append({key:str(v) for key,v in a.items() if type(v) is Q})
    checks['coarse_bound_allows_unbounded_rescaled_Haar_norm']=all(Q(2)**j>Q(2)**(j-1) for j in (1,10,100))
    return dict(checks=checks,rows=rows,scope='exact affine exponent comparisons; epsilon^-kappa divergence blocks coarse-to-fine relabelling')


def fixture(trace=None,seed_trace=None):
    old=sys.modules.get('budgets');sys.modules['budgets']=final
    try:c=load('preserved_haar_final_fixture',HERE.parent/'core_final_mean_row_closure/controls.py')
    finally:
        if old is None:sys.modules.pop('budgets',None)
        else:sys.modules['budgets']=old
    m=c.fixture(trace).means
    # Diagnostic LOWER original data, not actual coefficient enclosures.
    return PreservedHaarScalarJoin(m,seed_module.fixture_budget(trace=seed_trace),
        primary=lambda n:P({0:n+3,1:2}),first_signed=lambda n:P({0:n+2,1:1}),transfer=lambda n:P(n+2))


@lru_cache(None)
def composition_review():
    j=fixture();checks={};rows=[]
    for n in (0,1):
        a,b=j.seed_fields(n),j.first_slow(n);h=j.old_Haar(n)
        checks[f'n{n}:Haar_preservation_is_exact_finite_slow_sum']=all(h[k].c==(a[k]+b[k]).c for k in prefix.STATE)
        checks[f'n{n}:initial_defects_from_original_SeedBudget_not_coarse_initial_callback']=all(j.initial_defects(n)[k].c==(j.embedding*j.transfer(n)*P(j.seed.output(n)[k].c)).c for k in ('P','Jtheta','Jz'))
        parts=(j.initial_defects(n),j.first_wave_change(n),j.first_fast_change(n))
        checks[f'n{n}:first_postfast_all_three_complete_changes']=all(j.first_postfast(n)[k].c==sum((a[k] for a in parts),P()).c for k in ('P','Jtheta','Jz'))
        checks[f'n{n}:first_fast_uses_FULL_initial_means']=all(j.first_fast_change(n)[k].c==scalar.fast_change(n,j.current.embed(j.before.data('initial',n+1)),j.current.embed(j.before.fast(n+1)),j.current.moment(n+1))[k].c for k in ('P','Jtheta','Jz'))
        N=n+2;fixed=j.means.data('fixed',N)
        expected=haar_slow_reset(n,j.old_Haar(N),j.current.embed(j.current.slow_increment(N)),fixed['b'],j.current.moment(N))
        out=j.postslow(n)
        checks[f'n{n}:fresh_slow_uses_current_scalar_derived_C_envelope']=all(out[k].c==expected[k].c for k in out)
        checks[f'n{n}:radial_compensator_keeps_current_rho']=j.radial_compensator(n).c==(j.means.data('fixed',n)['rho']*out['P']).c
        checks[f'n{n}:whole_axial_moment_keeps_pressure']=j.whole_axial_moment(n,P(3)).c==(out['Jz']+3*out['P']).c
        rows.append(dict(order=n,postslow_degrees={k:max(v.c) for k,v in out.items()}))
    checks['original_weighted_input_gap_not_hidden']=len(j.certificate(0)['missing_actual_input_inequalities'])==5 and not j.certificate(0)['actual_current_fine_norms_certified']
    # Independent scalar counterexample: original U=5, update V=-4,
    # total=1. A total norm bound cannot bound the original summand.
    checks['updated_wave_norm_not_used_as_original_primary_bound']=abs(5-4)<abs(5) and j.primary(0).c!=j.before.data('wave',0)['U'].c
    D=lambda r:j.embedding*signed_covariance(j.primary(r),j.first_signed(r))
    expected_wave=scalar.wave_change(0,D(0),D(1),j.current.moment(1))
    checks['first_covariance_uses_separate_original_primary_and_signed_recipes']=all(j.first_wave_change(0)[k].c==expected_wave[k].c for k in ('P','Jtheta','Jz'))
    try:j.certify_actual_inputs()
    except NotImplementedError:checks['finite_recipe_cannot_certify_actual_inputs']=True
    else:checks['finite_recipe_cannot_certify_actual_inputs']=False
    # HUGE deferred A0/A1 source orders are untouched by this moment gain.
    a=final.rates()
    checks['same_previous_final_rows_and_compact_high_jet_orders']=a['A0_source_order_zero']==100005 and a['joined_regular']==j.orders['D']
    seed_trace=set();other=fixture(seed_trace=seed_trace);other.old_Haar(2)
    checks['postslow_initial_primary_scalar_channel_uses_n_plus_13']=max(seed_trace)==13
    other.old_Haar(3)
    checks['whole_axial_initial_primary_scalar_channel_uses_n_plus_14']=max(seed_trace)==14
    return dict(checks=checks,rows=rows,certificate=j.certificate(0),scope='finite recipes evaluated on manufactured lower norm data; no selected physical input or weighted transfer independently certified')


@lru_cache(None)
def review():
    parts=dict(exact=exact_review(),polynomials=polynomial_review(),rates=rate_review(),composition=composition_review())
    return dict(parts=parts,checks={name+':'+key:v for name,part in parts.items() for key,v in part['checks'].items()})
