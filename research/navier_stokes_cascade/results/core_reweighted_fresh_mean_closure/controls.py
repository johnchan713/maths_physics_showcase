"""Exact kernel domination and fresh full-pressure mean composition controls."""
from fractions import Fraction as Q
from functools import lru_cache
from budgets import HERE,load,join,prefix,P,X,convert,monomial_upper,cut_recipe,shift_recipe,materialize_shift,translate,FreshMeanBudget


def fixture(trace=None):
    tangent=join.tangent
    p=tangent.normal.frame.inputs(Rmin=1,Rmax=2,bmin=Q(1,16),bmax=Q(1,8),Fmax=Q(1,64),
        gmax=Q(1,32),gmin=Q(1,64),u=Q(1,16384),lam_min=Q(1,1024),lam_max=Q(1,128),
        cmin=1,cmax=1,Dbox=6,Db=Q(1,1048576),Mz=Q(1,1024))
    sp=tangent.normal.selector.fixture();L=tangent.normal.normal_length_certificate(p)['length_rule_floor']
    sp['length']=max(Q(25),Q(-(-L.numerator//L.denominator)))
    ell=tangent.normal.selector.select_tail(sp)['ell_tail'] # manufactured ONLY
    def operators(n):return prefix.haar.weighted_operator_budget(2,4,Q(1,2),P(n+2),P(n+3),1)
    C2=join.edge.rational_upper(cut_recipe(2,4,Q(1,2),1))
    C4=join.edge.rational_upper(cut_recipe(4,4,Q(1,2),1))
    E2=join.edge.atlas.edge_constant(2,Q(1,2));E4=join.edge.atlas.edge_constant(4,Q(1,2))
    def fixed(n):
        d={k:P(i+n+1) for i,k in enumerate(prefix.FIXED)}
        return dict(d,Kp=P(n+2)*(2+C2+C4),Kmom=P(n+3)*(E2+E4))
    b=prefix.PrefixBudget(fixed=fixed,initial=lambda n:{k:P(i+n+1) for i,k in enumerate(prefix.STATE)},
        wave=lambda n:{k:P(i+n+1) for i,k in enumerate(('T','R','Pi','U'))},
        linear=lambda n:{k:P(i+n+1) for i,k in enumerate(prefix.LINEAR)},
        KT1=lambda n:operators(n)['KT'],KA1=lambda n:operators(n)['KA'],time_inverse=lambda n:6*X,
        slow=lambda n:{k:P(i+n+1) for i,k in enumerate(('H','Kp','Kz','Kg','Kpsi'))},trace=trace)
    T=join.P
    w=join.GeneratedWave(p,sp,ell,b,F=lambda n:T(n+2),G=lambda n:T(n+3),btilde=lambda n:T(n+1),
                        path=lambda n:T({1:1}) if not n else T(int(n==1)))
    return FreshMeanBudget(w,1,4,Q(1,2),1)


def kernel_review():
    checks={};rows=[]
    for p,L,width,gap in ((Q(0),4,Q(1,2),1),(Q(1),4,Q(1,2),1),(Q(3,2),6,Q(1,3),Q(1,2))):
        recipe=shift_recipe(p,L,width,gap);K=materialize_shift(recipe)
        C={q:join.edge.rational_upper(cut_recipe(q,L,width,gap)) for q in (6*p,12*p)}
        E={q:join.edge.atlas.edge_constant(q,width) for q in (6*p,12*p)}
        tag=str(p)
        checks[tag+':one_fixed_shift_covers_all_four_new_domains']=all(K>=v for v in
            (1+C[6*p],E[6*p],2+C[6*p]+C[12*p],E[6*p]+E[12*p]))
        for q in (6*p,12*p):
            old=prefix.haar.weighted_operator_budget(q,L,width,P(1),P(1),gap)
            checks[tag+':pole'+str(q)+':same_actual_gap_kernel_recipe']=old['KT'].evaluate(1)==1+C[q] and old['KA'].evaluate(1)==old['weight_floor_inverse']*E[q]
        checks[tag+':same_original_width_and_gap']=recipe['width']==str(width) and recipe['gap']==str(Q(gap))
        rows.append(dict(pole=tag,shift_numerator_bits=K.numerator.bit_length(),shift_denominator_bits=K.denominator.bit_length()))
    K=materialize_shift(shift_recipe(1,4,Q(1,2),1));f=P({0:2,1:3,3:5})
    checks['exact_binomial_translation_no_new_S_ge_K_gate']=K>1 and translate(f,K).evaluate(1)==f.evaluate(1+K)
    checks['positive_polynomial_inputs_monotone_under_estimate_shift']=all(f.evaluate(s)<=f.evaluate(s+K) for s in (1,2,10))
    hull=monomial_upper(f)
    checks['sparse_polynomial_hull_dominates_every_nonnegative_monomial']=hull.c=={3:Q(10)} and all(hull.evaluate(s)>=f.evaluate(s) for s in (1,2,10))
    checks['zero_sparse_hull_stays_zero']=monomial_upper(P()).c=={}
    huge=shift_recipe(1,4,Q(1,100),Q(1,1000))
    checks['giant_constant_retained_exactly_without_materialization']=max(t['power'] for t in huge['terms'])>10000
    try:materialize_shift(huge)
    except ValueError:checks['materializer_guard_is_explicit_diagnostic_limit']=True
    else:checks['materializer_guard_is_explicit_diagnostic_limit']=False
    return dict(checks=checks,rows=rows,scope='exact rational upper recipes, not sampled continuum operator norms')


@lru_cache(None)
def fresh_review():
    trace=set();f=fixture(trace);g=f.generated;checks={}
    pre=f.rows(0);post=f.rows(0,True);final=f.final(0);end=f.final_rows(0)
    old=f.before.final(0);fast=f.fast(0);slow=f.slow_increment(0)
    checks['old_accumulated_means_survive_fresh_fast_and_slow']=all(final[k].evaluate(1)==(old[k]+fast[k]+slow[k]).evaluate(1) and final[k].evaluate(2)>=(old[k]+fast[k]+slow[k]).evaluate(2) for k in prefix.STATE)
    complete=f.before.wave(0)['U']+convert(g.addition(0)['U'])
    checks['complete_changed_wave_keeps_old_cross_and_sourced_square']=f.wave(0)['U'].c==monomial_upper(complete).c and f.wave(0)['U'].evaluate(1)==complete.evaluate(1)
    checks['all_three_scalars_and_pressure_recomputed_after_fast']=all(post[k].evaluate(1)>pre[k].evaluate(1) for k in ('gr','P','pm','Jtheta','Jz'))
    fresh=f.rows(0,True);nxt=f.rows(1,True);k=f.slow(0);kn=f.slow(1)
    expected=dict(u=k['H']*(k['Kp']*fresh['P']+k['Kz']*fresh['Jz']),z=k['H']*k['Kg']*fresh['Jtheta'],x=kn['H']*kn['Kpsi']*nxt['Jtheta'])
    checks['slow_map_uses_fresh_postfast_pressure_and_all_three_scalars']=all(slow[v].c==expected[v].c for v in prefix.STATE)
    wave=f.wave(3)['U'];independent=prefix.full_mean(0,dict(f.fixed(3),**f.final(3),w=wave*wave))
    checks['final_pressure_and_scalars_recomputed_after_slow']=all(end[v].c==independent[v].c for v in ('gr','P','pm','Etheta','Ez','Jtheta','Jz'))
    checks['postslow_pressure_contains_slow_mean_products']=end['pm'].evaluate(1)>post['pm'].evaluate(1)
    checks['full_pressure_Kp_and_global_moment_Kmom_reweighted']=all(f.fixed(0)[v].c==(f.before.fixed(0)[v]*X).c for v in ('Kp','Kmom'))
    # A1 is the KA1 axial mean operation, not a disposable zero-mean reset.
    time=f.time_inverse(0);rows=f.rows(4)
    without_A1=24*24*2**4*time*(rows['Etheta']+rows['Ez'])
    checks['axial_A1_operation_retained_in_fast_tangential_bound']=f.fast(0)['z'].evaluate(1)>without_A1.evaluate(1)
    high=max(n for name,n in trace if name=='initial')
    checks['postslow_order_zero_uses_finite_original_mean_jets_through_31']=high==31
    checks['full_cut_sourced_pressure_retained']=g.addition(0)['Pi'].evaluate(1)>0
    checks['normal_and_tangent_from_original_frame_inputs']=g.certificate['tangent_and_length_inputs_derived_from_retained_frame_bounds']
    checks['no_actual_input_enclosure_or_residual_promotion']=not g.certificate['actual_primitive_enclosures_certified'] and not g.certificate['complete_covariance_means_or_residual_certified']
    for name in ('source','raw_source','inverse_source'):
        try:getattr(f,name)(0)
        except NotImplementedError:checks['next_harmonic_grammar_blocks_'+name]=True
        else:checks['next_harmonic_grammar_blocks_'+name]=False
    return dict(checks=checks,highest_initial_jet=high,final_mean_degrees={v:max(final[v].c) for v in prefix.STATE},
        final_row_degrees={v:max(end[v].c) for v in ('gr','P','pm','Etheta','Ez','Jtheta','Jz')},
        scope='exact positive-polynomial composition on manufactured original recipe inputs; size only')


@lru_cache(None)
def review():
    parts=dict(kernel=kernel_review(),fresh=fresh_review())
    return dict(checks={name+':'+k:v for name,r in parts.items() for k,v in r['checks'].items()},parts=parts)
