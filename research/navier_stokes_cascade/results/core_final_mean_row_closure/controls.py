"""Independent exact pressure/compact-mode and mixed-rate diagnostics."""
from fractions import Fraction as Q
from functools import lru_cache
import sys
from budgets import HERE,load,haar,fast,prefix,P,rates,oscillatory_recipe,remaining_row_envelope,FinalMeanRowClosure

ops=load('final_mean_original_operators',HERE.parent/'core_fast_inverse_budget/proof_checks.py')
means=ops.means;U,Z,T=ops.U,ops.Z,ops.T


@lru_cache(None)
def compact_review():
    shape=(U-1)**6*(2-U)**6;mass=shape.integral().scalar()
    center=(U*shape).integral().scalar()/mass
    rho=shape/mass*(1+Z*(U-center));chi=(shape.primitive()-shape.primitive().at(0,1))/mass+Z*shape/100
    grbar=shape*(1+T+Z+(U-center)**2);Pnew=grbar.integral()
    hbar=grbar-rho*Pnew
    pm=hbar.primitive()-hbar.primitive().at(0,1)-chi*hbar.integral()
    f=shape*(1+Z+Z**2+T);omega=Q(5,2)
    K,J,I=ops.compact(ops.CP(f),chi,omega)
    A=ops.mul(chi.diff(0),J);F=(ops.CP(f),ops.CP(),ops.CP())
    residual=ops.add(ops.graph(K,omega),ops.neg(F))
    checks=dict(
        full_moving_compensator_has_exact_unit_mass=rho.integral()==1 and bool(rho.diff(1)),
        compensated_Haar_source_has_zero_whole_radial_moment=hbar.integral()==0 and bool(hbar),
        Haar_pressure_field_not_zero_or_order_promoted=bool(pm),
        Haar_radial_residual_is_exactly_minus_rho_P=pm.diff(0)-grbar==-rho*Pnew,
        compensated_A0_zero_mode_exactly_zero=chi.diff(0)*hbar.integral()==0,
        nonzero_mode_compact_identity_exact=ops.zero(ops.add(ops.graph(K,omega),ops.neg(ops.add(F,ops.neg(A))))),
        pointwise_A0_nonzero_despite_zero_Haar=not ops.zero(A),
        full_pointwise_radial_nonzero_residual_is_minus_A0=residual==ops.neg(A),
        compact_pressure_both_zero_mode_endpoints=pm.at(0,1)==pm.at(0,2)==0,
        compact_pressure_nonzero_mode_endpoints=K[0].at(0,1)+K[1].at(0,1)==0 and K[2].at(0,1)==0 and K[0].at(0,2)+K[2].at(0,2)==0 and K[1].at(0,2)==0,
        moving_cutoff_derivatives_commute_with_graph=ops.graph(ops.diff(K,1),omega)==ops.diff(ops.graph(K,omega),1),
        full_pressure_compensation_stays_zero_after_Z_derivative=hbar.diff(1).integral()==0,
        freezing_rho_derivative_is_detected=(grbar.diff(1)-rho*Pnew.diff(1)).integral()==rho.diff(1).integral()*Pnew and bool(rho.diff(1)*Pnew),
    )
    badrho=(1+Z)*rho;badh=grbar-badrho*Pnew
    checks['nonunit_compensator_leaves_nonflat_A0_zero_mode']=badh.integral()==-Z*Pnew and bool(chi.diff(0)*badh.integral())
    checks['omitting_compensation_leaves_nonflat_A0_zero_mode']=bool(chi.diff(0)*grbar.integral())
    checks['zero_Haar_is_not_pointwise_zero']=not ops.zero(A) and bool(rho*Pnew)
    # The slow pressure derivative includes the missing rho_Z P term
    # pointwise even though its integral alone is zero.
    frozen=grbar.diff(1)-rho*Pnew.diff(1)
    checks['pointwise_pressure_source_needs_rho_Z_P']=frozen-hbar.diff(1)==rho.diff(1)*Pnew and bool(frozen-hbar.diff(1))
    checks['e0_U_to_R_jacobian_restores_exact_radial_moment']=all(
        Q(rv)**(1-d)/d*(d*Q(rv)**(d-1))==1 for d in (1,2,3) for rv in (1,2,3))
    return dict(checks=checks,scope='formal nonzero endpoint phases at rational frequency and finite C5 compact polynomial data, not the selected C-infinity profiles or their huge jets')


def physical_difference(n,old,inc,fixed,eps):
    """Independent raw physical-size calculation with epsilon factors kept."""
    d=Q(n+4);b,V,G,g,c,r,rho,kp,km=(fixed[k] for k in ('b','V','G','graph','c','r','rho','Kp','Kmom'))
    x,u,z=(old[k] for k in prefix.STATE);X,Uv,Zv=(inc[k] for k in prefix.STATE)
    dr=d*g;vis=d*d*(g*g+1)+r*dr+r*r
    srr=2*b*X+2*x*X+X*X;srz=b*Zv+G*X+x*Zv+X*z+X*Zv
    stt=2*V*Uv+2*u*Uv+Uv*Uv
    gr=d*(eps+c)*X+(dr+r)*srr+eps*d*srz+r*stt+eps*vis*X
    pressure=kp*(gr+rho*km*gr)
    rt=b*Uv+V*X+x*Uv+X*u+X*Uv
    at=G*Uv+V*Zv+z*Uv+Zv*u+Zv*Uv
    az=2*G*Zv+2*z*Zv+Zv*Zv+pressure
    return dict(Etheta=eps*d*Uv+(dr+2*r)*rt+eps*d*at+eps*vis*Uv,
                Ez=eps*d*Zv+(dr+r)*srz+eps*d*az+eps*vis*Zv,pm=pressure,gr=gr)


@lru_cache(None)
def rate_review():
    checks={};rows=[]
    for k in (Q(1,100000),Q(1,10),Q(249,1000)):
        a=rates(k);B,C,D=a['B'],a['C'],a['D'];tag=str(k)
        checks[tag+':both_fresh_difference_orders_stronger_than_C']=a['Q1_terms']['fresh_fast_difference']>C and a['Q1_terms']['fresh_slow_difference']>C
        checks[tag+':only_two_fresh_A1_integrations_needed']=a['Q1_terms']['fresh_A1']==Q(3,2)>C and a['final_Q1']==Q(3,2)
        checks[tag+':strict_A0_margin_and_minimal_integer']=B+a['A0_integration_order']*k>C and B+(a['A0_integration_order']-1)*k<=C
        checks[tag+':all_three_final_regular_mean_rows_C']=a['final_tangential']==a['final_radial']==a['full_angular_zero_regular']==C
        checks[tag+':radial_P_and_A0_both_stronger_than_C']=a['radial_compensator']>C and a['A0_order']>C
        checks[tag+':full_radial_source_B_derived_with_covariance_graph_and_centripetal_terms']=a['radial_source']==a['radial_source_groups']['covariance_graph']==a['radial_source_groups']['centripetal']==B
        checks[tag+':combined_regular_still_D_not_full_physical']=a['joined_regular']==D<C and not a['complete_physical_tail_or_blowup_certified']
        rows.append({key:str(a[key]) for key in ('B','C','D','final_Q1','A0_integration_order','A0_source_order_zero')})
    # eps=t^10, kappa=1/10, B=.9, C=1.4, D=1.3. Independent
    # physical expression retains epsilon on slow time, viscosity, axial
    # pressure and fluxes; the majorant drops only nonnegative gains.
    old={k:P(i+2) for i,k in enumerate(prefix.STATE)}
    finc={k:P(i+5) for i,k in enumerate(prefix.STATE)}
    sinc={k:P(i+8) for i,k in enumerate(prefix.STATE)}
    fixed={k:P(2) for k in ('b','V','G','graph','c','r','rho','Kp','Kmom')}
    for n in (0,1,2):
        bound=remaining_row_envelope(n,old,finc,sinc,fixed)
        for t in (Q(1,2),Q(1,8)):
            ev=lambda data:{k:v.evaluate(1) for k,v in data.items()}
            base=ev(fixed);base.update(b=2*t**10,graph=2/t)
            o={k:v*t**(19 if k=='x' else 9) for k,v in ev(old).items()}
            f={k:v*t**(23 if k=='x' else 13) for k,v in ev(finc).items()}
            sl={k:v*t**(24 if k=='x' else 14) for k,v in ev(sinc).items()}
            raw=physical_difference(n,o,f,base,t**10)
            rawslow=physical_difference(n,{k:o[k]+f[k] for k in prefix.STATE},sl,base,t**10)
            for component in ('Etheta','Ez'):
                checks[f'n{n}_t{t}:independent_fresh_fast_'+component]=raw[component]<=t**21*bound['fast_difference'][component].evaluate(1)
                checks[f'n{n}_t{t}:independent_fresh_slow_'+component]=rawslow[component]<=t**22*bound['slow_difference'][component].evaluate(1)
            checks[f'n{n}_t{t}:pressure_not_falsely_promoted']=raw['pm']<=t**13*bound['fast_difference']['pm'].evaluate(1)
            # Full current source: include covariance, not just an update.
            w=Q(3)*t**10;d=Q(n+4);dr=d*base['graph'];rr=base['r'];b=base['b']
            xv,uv,zv=(o[k] for k in prefix.STATE)
            visc=d*d*(base['graph']**2+1)+rr*dr+rr*rr
            rawgr=d*(t**10+base['c'])*xv+(dr+rr)*(2*b*xv+xv*xv+w)
            rawgr+=t**10*d*(b*zv+base['G']*xv+xv*zv+w)+rr*(2*base['V']*uv+uv*uv+w)+t**10*visc*xv
            normalized=dict(ev(fixed),**ev(old),sigma=Q(1),w=Q(3))
            coarse=prefix.full_mean(n,normalized)['gr']
            checks[f'n{n}_t{t}:independent_full_radial_source_keeps_covariance_and_graph']=rawgr<=t**9*coarse
    a=rates()
    checks['A0_zero_output_needs_current_source_100005']=a['A0_source_order_zero']==100005
    return dict(checks=checks,rates=rows,scope='exact rational physical term budgets; actual high-jet coefficient bounds remain analytical premises')


def fixture(trace=None):
    old=sys.modules.get('budgets');sys.modules['budgets']=haar.fresh
    try:c=load('final_mean_fresh_fixture',HERE.parent/'core_reweighted_fresh_mean_closure/controls.py')
    finally:
        if old is None:sys.modules.pop('budgets',None)
        else:sys.modules['budgets']=old
    return FinalMeanRowClosure(c.fixture(trace))


@lru_cache(None)
def composition_review():
    trace=set();j=fixture(trace);a=j.certificate(0);t,r=a['tangential'],a['radial'];checks={}
    checks['final_stage_does_not_reuse_postwave_row_as_final']=t['stage']==r['stage']=='after-fresh-fast-and-slow' and 'Rfast+Q1 Rslow' in t['temporal_identity']
    checks['sharp_fast_D_is_derived_not_old_B_relabelled']=t['sharp_fast_increment']==j.fast.fast_recipe(3)
    checks['fresh_slow_C_uses_same_scalar_derived_recipe']=t['sharp_slow_increment']['method']=='ScalarGainBudget.slow_increment' and t['sharp_slow_increment']['order']==3 and t['sharp_slow_increment']['embedding']==str(j.embedding)
    checks['fresh_A1_two_integrations_uses_postwave_source_seven']=t['exact_A1']['integration_order']==2 and t['exact_A1']['raw_source_order']==7
    nested=t['exact_A1']['source_recipe']['exact_A1']
    checks['nested_original_A1_high_jet_not_discarded']=nested['raw_source_order']==200008 and 'original_geometry' in nested
    checks['fast_difference_radial_input_retains_old_source_200009']=t['sharp_fast_increment']['radial_rows']['exact_A1']['raw_source_order']==200009
    a0=r['exact_A0']
    checks['A0_source_is_current_final_full_graph_gr']=a0['raw_source_order']==100005 and a0['source_recipe']['method']=='FinalMeanRowClosure.current_radial_source_polynomial' and a0['source_recipe']['component']=='gr'
    checks['A0_source_high_jets_deferred_not_sampled_or_allocated']=a0['source_recipe']['field_input_order']==100008 and not a0['materialized_actual_input_bound']
    checks['A0_uses_e0_raw_profiles_not_A1_operator_norm']=a0['operator']=='A0' and 'k_0' in a0['factors'][1] and 'f_0' in a0['factors'][2]
    checks['source_pole_12p_and_original_floor_width_gap_retained']=a0['source_pole']=='12' and a0['original_geometry']==j.means.reweighting and a0['inverse_output_edge_floor_power']==5
    fields={k:j.embedding*v for k,v in j.means.final(3).items()};wave=j.means.data('wave',3)['U']
    expected=prefix.full_mean(0,dict(j.means.data('fixed',3),**fields,w=j.embedding*wave*wave))['gr']
    checks['current_gr_lower_pole_embedding_precedes_products']=j.current_radial_source_polynomial(0).c==expected.c and a0['source_recipe']['embedding']==str(j.embedding)
    checks['actual_A1_velocities_still_in_fast_increment']=bool(t['sharp_fast_increment']['A1_velocity_factor'])
    checks['radial_compensator_is_actual_postslow_P_recipe']=r['compensator']['method']=='HaarStressJoin.radial_Haar_compensator' and r['compensator']['order']==0
    checks['no_pointwise_pressure_order_promotion']=r['full_pressure_field_not_order_promoted']
    checks['exact_additive_tails_excluded_not_forwarded']='excluded' in t['additive_tail_ledger']
    checks['all_regular_mean_rows_C_not_complete_physical_residual']=a['angular_zero_regular_order']==str(j.orders['C']) and not a['actual_input_validity_or_full_physical_residual']
    try:j.actual_input_certificate()
    except NotImplementedError:checks['deferred_recipe_cannot_certify_actual_inputs']=True
    else:checks['deferred_recipe_cannot_certify_actual_inputs']=False
    # Directly check the finite constants against the independent original
    # moderate integration recipe, avoiding the giant selected integer.
    recipe=oscillatory_recipe(1,2,e=0,source='test',source_rate=Q(1),target_rate=Q(11,10),
        kappa=Q(1,10),geometry=dict(L='4',width='1/2',gap='1'),pole=12)
    raw=ops.flat_budget_data(1,2,4,2,3)
    ours=recipe['raw_output_scalar']*2**recipe['two_power']*3**recipe['three_power']*4*2*3
    checks['moderate_original_Fourier_constant_matches_raw_times_Q1_and_output_conversion']=ours==2*6*raw['constant']
    return dict(checks=checks,A0_recipe=a0,largest_nested_old_source_order=200009,
        largest_original_initial_jet_in_improved_tangential_recipe=200012,
        scope='finite exact composition on manufactured original operator recipes; deferred huge inputs are not evaluated or verified')


@lru_cache(None)
def review():
    parts=dict(compact=compact_review(),rates=rate_review(),composition=composition_review())
    return dict(parts=parts,checks={name+':'+key:v for name,part in parts.items() for key,v in part['checks'].items()})
