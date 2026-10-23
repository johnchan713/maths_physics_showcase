"""Independent full row differences, raw epsilon budgets and dependency traces."""
from fractions import Fraction as Q
from functools import lru_cache
import sys
from budgets import HERE,load,signed,final,preserved,scalar,prefix,fast,fresh,P,rates,fast_envelope,postfast_scalar_envelope,old_Haar_envelope,final_row_envelope,mean_feedback,operator_ledger,mean_shift,dependency_trace,PostSignedMeanRecipes

a=load('post_signed_exact_algebra',HERE.parent/'current_target_moment_audit/target_algebra.py')
s,R,Z,T,u=a.s,a.R,a.Z,a.T,a.u


@lru_cache(None)
def exact_review():
    ep=a.eps;V=a.V;b=a.b;G=s.S(0);q=a.state()
    def engine(q,p):
        B,v,g,W=(q[k] for k in ('beta','v','gamma','W'))
        gr=-a.ts(B)-a.de(2*b*B+B*B+W['rr'],1)-a.dz(b*g+G*B+B*g+W['rz'])+(2*V*v+v*v+W['tt'])/R+ep*(a.lap(B)-B/R**2)
        Et=a.ts(v)+a.de(b*v+V*B+B*v+W['rt']-a.st,2)+a.dz(G*v+V*g+g*v+W['tz'])-ep*(a.lap(v)-v/R**2)
        Ez=a.ts(g)+a.de(b*g+G*B+B*g+W['rz']-a.sz,1)+a.dz(2*G*g+g*g+W['zz']+p)-ep*a.lap(g)
        return dict(gr=a.ex(gr),Et=a.ex(Et),Ez=a.ex(Ez))
    def increment(scale,isfast):
        c=(u+1/u)/2 if isfast else 1
        psi=scale*R**3*a.F*(1+Z+Z*T)*c
        return dict(beta=-a.dz(psi),gamma=a.de(psi,1),v=scale*(5*R**2*a.F+R**3*s.diff(a.F,R))*(1+Z+T)*c)
    def plus(q,d):return dict(W=q['W'],**{k:a.ex(q[k]+d[k]) for k in ('beta','v','gamma')})
    def difference(old,d,dp):
        B,v,g=(old[k] for k in ('beta','v','gamma'));X,Uv,Zv=(d[k] for k in ('beta','v','gamma'))
        S=2*b*X+2*B*X+X*X;H=b*Zv+G*X+B*Zv+X*g+X*Zv;St=2*V*Uv+2*v*Uv+Uv*Uv
        gr=-a.ts(X)-a.de(S,1)-a.dz(H)+St/R+ep*(a.lap(X)-X/R**2)
        Rt=b*Uv+V*X+B*Uv+X*v+X*Uv;At=G*Uv+V*Zv+g*Uv+Zv*v+Zv*Uv
        Az=2*G*Zv+2*g*Zv+Zv*Zv+dp
        # Exact remainder AFTER c Nt tangential term is separated.
        et=-ep*s.diff(Uv,T)+a.de(Rt,2)+a.dz(At)-ep*(a.lap(Uv)-Uv/R**2)
        ez=-ep*s.diff(Zv,T)+a.de(H,1)+a.dz(Az)-ep*a.lap(Zv)
        return dict(gr=a.ex(gr),Et=a.ex(et),Ez=a.ex(ez))
    # Pressure outputs are arbitrary COMPLETE smooth fields here. The
    # identity holds for every pair, hence for the actual linear T0
    # reconstruction. Its compact/time/moving identities are tested
    # separately below; these fields are not claimed to solve T0.
    p0=R**2*a.F*(1+Z+T+(u+1/u)/2)
    pf=R**2*a.F*(1+Z**2+T)*(2+u+1/u)
    ps=R**2*a.F*(1+Z*T+Z**3)
    f=increment(ep**2,True);sl=increment(ep**3,False)
    qf=plus(q,f);qe=plus(qf,sl)
    old,mid,end=engine(q,p0),engine(qf,p0+pf),engine(qe,p0+pf+ps)
    rf,rs=difference(q,f,pf),difference(qf,sl,ps)
    checks={}
    for name,oldrow,newrow,rem,d in (('fast',old,mid,rf,f),('slow',mid,end,rs,sl)):
        checks[name+':complete_current_radial_source_difference']=a.zero(newrow['gr']-oldrow['gr']-rem['gr'])
        for key,field in (('Et','v'),('Ez','gamma')):
            checks[name+':whole_tangential_difference_'+key]=a.zero(newrow[key]-oldrow[key]-a.ci*a.L(d[field])-rem[key])
        checks[name+':actual_three_component_divergence']=a.zero(a.de(d['beta'],1)+a.dz(d['gamma']))
        checks[name+':both_weighted_Haar_velocity_moments']=a.zero(a.integ(R**2*a.haar(d['v']))) and a.zero(a.integ(R*a.haar(d['gamma'])))
    checks['all_new_fast_components_zero_Haar']=all(a.zero(a.haar(f[k])) for k in f)
    checks['new_slow_has_zero_temporal_principal_term']=all(a.zero(a.L(sl[k])) for k in sl)
    checks['old_new_fast_correlations_survive']=not a.zero(a.haar(q['v']*f['v']))
    checks['new_slow_Haar_old_field_is_old_Haar_plus_previous_slow']=all(a.zero(a.haar(qf[k])-a.haar(q[k])) for k in ('beta','v','gamma'))
    checks['nonzero_reserved_axial_background_can_leave_linear_scalar_leaks']=not a.zero(a.integ(R**3*sl['v'])) and not a.zero(a.integ(2*R**2*sl['gamma']))
    checks['full_recomputed_pressure_increment_derivatives_not_erased']=not a.zero(a.dz(pf)) and not a.zero(a.dz(ps))
    checks['all_six_complete_covariances_remain_fixed_through_means']=all(q['W'][k]==qf['W'][k]==qe['W'][k] and not a.zero(q['W'][k]) for k in q['W'])
    # Full row identity follows from the two independently expanded
    # differences and the original exact fast inverse, not a Q1-only
    # estimate accidentally labelled full.
    for key,field in (('Et','v'),('Ez','gamma')):
        checks['complete_two_step_identity_'+key]=a.zero(end[key]-old[key]-a.ci*a.L(f[field])-rf[key]-rs[key])
    ops=load('post_signed_raw_compact',HERE.parent/'core_fast_inverse_budget/proof_checks.py')
    original=ops.mode_review()
    checks['original_exact_A1_time_error_has_positive_residual_sign']=original['checks']['actual_fast_time_error_is_A1_of_current_residual'] and original['checks']['wrong_fast_time_error_sign_detected']
    checks['original_actual_A1_velocity_and_moving_cutoff_are_retained']=original['checks']['actual_A1_velocity_nonzero'] and original['checks']['moving_cutoff_Z_term_is_required']
    return dict(checks=checks,scope='exact graph/Laurent conservative row differences with arbitrary complete pressure outputs, plus separate original formal-phase compact inverse controls; not an actual continuum five-map realization')


@lru_cache(None)
def budget_review():
    old=sys.modules.get('budgets');sys.modules['budgets']=final
    try:fc=load('post_signed_raw_physical',HERE.parent/'core_final_mean_row_closure/controls.py')
    finally:
        if old is None:sys.modules.pop('budgets',None)
        else:sys.modules['budgets']=old
    checks={};rows=[]
    oldfields={k:P(i+2) for i,k in enumerate(prefix.STATE)}
    fi={k:P(i+5) for i,k in enumerate(prefix.STATE)};sl={k:P(i+8) for i,k in enumerate(prefix.STATE)}
    fixed={k:P(2) for k in ('b','V','G','graph','c','r','rho','Kp','Kmom')}
    for n in (0,1,2):
        out=final_row_envelope(n,dict(Etheta=P(3),Ez=P(4)),P(5),oldfields,fi,sl,fixed)
        for t in (Q(1,2),Q(1,8)):
            ev=lambda d:{k:v.evaluate(1) for k,v in d.items()}
            base=ev(fixed);base.update(b=2*t**10,graph=2/t)
            o={k:v*t**(19 if k=='x' else 9) for k,v in ev(oldfields).items()}
            f={k:v*t**(23 if k=='x' else 13) for k,v in ev(fi).items()}
            ss={k:v*t**(24 if k=='x' else 14) for k,v in ev(sl).items()}
            rawf=fc.physical_difference(n,o,f,base,t**10)
            raws=fc.physical_difference(n,{k:o[k]+f[k] for k in o},ss,base,t**10)
            for key in ('Etheta','Ez'):
                haar=3 if key=='Etheta' else 4;tail=0 if key=='Etheta' else 5*t**19
                raw=t**18*haar+tail+rawf[key]+raws[key]
                checks[f'n{n}_t{t}:full_regular_mean_row_bound_'+key]=raw<=t**18*out[key].evaluate(1)
            checks[f'n{n}_t{t}:full_pressure_not_promoted_with_row_gain']=rawf['pm']<=t**13*out['fast_difference']['pm'].evaluate(1)
        # All new scalar norms are derived sums of old, covariance and
        # FULL old/fast correlations. Their functionals agree exactly.
        oldS={k:P(i+2) for i,k in enumerate(('P','Jtheta','Jz'))}
        post=postfast_scalar_envelope(n,oldS,P(4),P(5),oldfields,fi,P(6))
        direct=postfast_scalar_envelope(n,{k:v.evaluate(2) for k,v in oldS.items()},Q(4),Q(5),{k:v.evaluate(2) for k,v in oldfields.items()},{k:v.evaluate(2) for k,v in fi.items()},Q(6))
        checks[f'n{n}:complete_postfast_scalar_polynomial_substitution']=all(v.evaluate(2)==direct[k] for k,v in post.items())
        before=old_Haar_envelope(oldfields,sl)
        checks[f'n{n}:preserved_Haar_recipe_keeps_previous_slow']=all(before[k].c==(oldfields[k]+sl[k]).c for k in before)
        fb=mean_feedback(n,P(3),fi,graph=P(2),reciprocal=P(2),normal=P(2))
        for t in (Q(1,2),Q(1,8)):
            d=n+4;dm=sum(v.evaluate(1) for v in fi.values());Uwave=3*t**5
            raw=(2*(2*d/t+d+2)+8/t**5)*dm*t**13*Uwave
            checks[f'n{n}_t{t}:derived_new_fast_mean_feedback_order_D']=raw<=t**13*fb.evaluate(1)
    for k in (Q(1,100000),Q(1,10),Q(249,1000)):
        r=rates(k);C0=r['C0'];tag=str(k)
        checks[tag+':fast_and_slow_full_difference_strictly_above_C0']=r['tangential_terms']['fast_difference']>C0 and r['tangential_terms']['slow_difference']>C0
        checks[tag+':all_regular_mean_rows_close_at_C0']=r['final_tangential']==r['final_radial']==r['final_regular_mean']==C0
        checks[tag+':fresh_current_postfast_scalars_stay_C']=r['postfast_scalar']==r['C']<r['B']+r['D'] and r['C']<2*r['D']
        checks[tag+':fresh_Haar_scalar_reset_C_plus_one']=r['postslow_scalar']==r['C']+1>C0
        checks[tag+':both_compact_counts_have_strict_margin']=r['tangential_terms']['A1']>C0 and r['A0']>C0
        checks[tag+':preceding_compact_integers_do_not_meet_this_strict_rule']=r['D']+(r['A1_integrations']-1)*k<=C0 and r['B']+(r['A0_integrations']-1)*k<=C0
        rows.append({key:str(v) for key,v in r.items() if type(v) in (int,Q)})
    return dict(checks=checks,rates=rows,scope='positive-polynomial functionals and independent rational physical epsilon term budgets; lower actual weighted bounds are not evaluated')


@lru_cache(None)
def recipe_review():
    checks={};rows=[];geometry=dict(L='4',width='1/2',gap='1')
    shift=mean_shift(30,1,geometry);led=shift['ledger'];cover=shift['recipe']['covered_operator_poles']
    checks['new_postfast_4q_domain_exceeds_signed_pressure_2q']=led['postfast_source']=='120' and led['post_signed_rows']=='60'
    checks['original_shift_proves_all_2q_4q_domains']=all(set(v)<=set(cover[k]) for k,v in led['required_domains'].items())
    checks['wave_only_pole_does_not_forget_old_mean_domain']=operator_ledger(3,1)['q']=='6'
    checks['same_original_width_hull_and_gap_in_operator_recipe']=all(shift['recipe'][k]==geometry[k] for k in geometry)
    K=fresh.materialize_shift(shift['recipe'])
    checks['moderate_new_kernel_constants_are_dominated']=all(K>=fresh.join.edge.rational_upper(fresh.cut_recipe(p,4,Q(1,2),1)) for p in (60,120))
    calls=[]
    def rows_at(n):calls.append(n);return dict(Etheta=P(n+2),Ez=P(n+3))
    ff=fast_envelope(1,rows_at,P(2),P(3),P(5))
    factor=6*24*2**5*5*120
    checks['fast_envelope_reads_tangential_n4_and_radial_n5']=calls==[5,6]
    checks['fast_envelope_keeps_A1_and_raw_factorial_conversion']=ff['z'].evaluate(1)==factor*4*8 and ff['u'].evaluate(1)==factor*7
    for k in (Q(1,100000),Q(1,10)):
        r=rates(k);t=dependency_trace(0,r['A1_integrations'],r['A0_integrations'])
        checks[str(k)+':new_next_target_input_is_A0_source_plus_14']=t['maxima']['next_signed_target_input']==t['A0_source_order']+14
        checks[str(k)+':new_prefast_row_is_A0_source_plus_10']=t['maxima']['post_signed_row']==t['A0_source_order']+10
        checks[str(k)+':older_DAGs_and_actual_profiles_not_evaluated']=not t['old_dependency_DAGs_expanded'] and not t['actual_jet_or_profile_values_evaluated']
        rows.append(t)
    checks['order_zero_A1_A0_and_next_target_counts']=rows[0]['A1_source_order']==100005 and rows[0]['A0_source_order']==200003 and rows[0]['maxima']['next_signed_target_input']==200017
    # Construct same-prefix metadata using a minimal stub: these recipes
    # never materialize the prior huge coefficient/norm DAG.
    class Common:
        reweighting=dict(geometry,original_pole='1')
    common=Common()
    class Final:
        means=common
    class Haar:
        means=common
    class Current:
        final=Final();haar=Haar()
        def target_inputs(self,n):return dict(stage='after-latest-fast-and-slow-before-next-generic-signed')
    j=PostSignedMeanRecipes(Current(),wave_pole=30,original_pole=1,geometry=geometry)
    c=j.certificate(0);r=j.post_signed_rows(0)
    checks['new_rows_derived_from_old_final_recipe_and_all_six_covariance_entries']=r['old_rows']['method']=='FinalMeanRowClosure.tangential_recipe' and 'all six' in r['covariance']['entries']
    checks['A1_actual_velocity_kept_in_fresh_fast_recipe']='g_d-A1 g_d' in c['fast_recipe']['actual_axial_velocity']
    q=Q(c['common_pole']['q'])
    checks['new_A1_source_is_post_signed_row_recipe']=c['compact']['A1']['source_recipe']['method']=='PostSignedMeanRecipes.post_signed_rows' and Q(c['compact']['A1']['source_pole'])==2*q
    checks['new_A0_source_is_FULL_final_gr_not_stale_gr']='NEW fast + NEW slow' in c['compact']['A0']['source_recipe']['method'] and Q(c['compact']['A0']['source_pole'])==4*q
    checks['one_common_pole_uses_whole_huge_request_graph']=q==1200114 and c['common_pole']['maximum_current_target_input']==200017 and q>30
    checks['huge_operator_constants_deferred_instead_of_allocated']=c['operator_shift']['recipe']['method']=='fresh.shift_recipe' and not c['operator_shift']['recipe']['materialized']
    checks['fast_difference_reads_increment_jets_through_n_plus_three']=c['fast_recipe']['order']==3 and c['fast_recipe']['radial_rows_order']==8
    try:PostSignedMeanRecipes(Current(),wave_pole=30,original_pole=1,geometry=dict(geometry,width='1/3'))
    except ValueError:checks['different_actual_width_or_hull_is_rejected']=True
    else:checks['different_actual_width_or_hull_is_rejected']=False
    checks['whole_new_axial_pressure_product_stays_current']='derivative of c_rho' in c['scalar_recipe']['whole_axial_pressure_moment']
    checks['all_high_jet_source_profiles_floor_and_raw_kernel_recipes_retained']=all('N_' in r['factors'][1] and 'N_' in r['factors'][2] and r['source_factorial']==r['raw_source_order'] and not r['materialized_actual_input_bound'] for r in c['compact'].values())
    try:j.certify_actual_inputs()
    except NotImplementedError:checks['no_actual_input_or_full_cycle_certificate_from_finite_recipes']=True
    else:checks['no_actual_input_or_full_cycle_certificate_from_finite_recipes']=False
    return dict(checks=checks,request_traces=rows,operator_ledger=led,certificate=c,scope='finite exact request graph, moderate kernel constants and deferred same-stage recipes; not actual high-jet input verification')


@lru_cache(None)
def review():
    parts=dict(exact=exact_review(),budgets=budget_review(),recipes=recipe_review())
    return dict(parts=parts,checks={name+':'+key:v for name,part in parts.items() for key,v in part['checks'].items()})
