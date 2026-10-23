"""Independent full Fourier difference, homogeneous pressure and request checks."""
from fractions import Fraction as Q
from functools import lru_cache
import sympy as s
from bounds import HERE,load,mean,mult,physical,signed,prefix,P,rates,homogeneous_pressure,regular_source_envelope,physical_request,whole_domain,CurrentResidualJoin

finite=load('post_residual_fourier',HERE.parent/'finite_induction_foundation/finite_checks.py')


def fixture():
    class Common:
        reweighting=dict(original_pole='1',L='4',width='1/2',gap='1')
    common=Common()
    class Final:means=common
    class Haar:means=common
    class Current:
        final=Final();haar=Haar()
        def target_inputs(self,n):return dict(stage='same-current-prefix-before-generic-signed')
    m=mean.PostSignedMeanRecipes(Current(),wave_pole=30,original_pole=1,geometry=common.reweighting)
    return CurrentResidualJoin(m,Q(1,200))


@lru_cache(None)
def pressure_review():
    u,R,Z,T=s.symbols('u R Z T',real=True);c,w,lam=s.symbols('c w lam',positive=True)
    n=s.Matrix([1+u,1,0]);t=s.Matrix([1,-1-u,2]);K=s.Matrix([[1,2,3],[4,5,6],[7,8,9]])
    A=(n.dot(K*t)-c*n.diff(u).dot(t));pi=s.I*A/(w*n.dot(n))
    principal=c*t.diff(u)+K*t+lam*t+s.I*w*n*pi
    wrongpi=s.I*n.dot(K*t)/(w*n.dot(n));wrong=c*t.diff(u)+K*t+lam*t+s.I*w*n*wrongpi
    alpha=1+R+Z*T
    checks=dict(
        principal_transversality_identity=n.dot(t)==0,
        full_moving_normal_pressure_cancels_normal_principal_residual=s.simplify(n.dot(principal))==0,
        omitting_moving_normal_pressure_is_detected=s.simplify(n.dot(wrong))!=0,
        homogeneous_pressure_is_not_zero=pi!=0,
        auxiliary_independent_current_target_scalar_preserves_homogeneous_scaling=all(s.simplify(v)==0 for v in (c*(alpha*t).diff(u)+K*(alpha*t)+lam*alpha*t+s.I*w*n*(alpha*pi)-alpha*principal)),
        auxiliary_dependent_scalar_leaves_a_principal_velocity_leak=all(s.simplify(v)==0 for v in (c*((1+u)*t).diff(u)+K*((1+u)*t)+lam*(1+u)*t+s.I*w*n*((1+u)*pi)-(1+u)*principal-c*t)),
        transverse_equation_still_requires_original_projected_solve=any(s.simplify(v)!=0 for v in principal),
    )
    for a in (Q(1),Q(3,2),Q(7)):
        checks[f'bound{a}:both_K_and_moving_normal_pressure_terms_retained']=homogeneous_pressure(a,normal=2,K=3,c=4,normal_time=5,inverse_normal_square=6)==156*a
    return dict(checks=checks,scope='exact normal projection and scalar homogeneity; the full transverse homogeneous solution is a retained original pulse input, not solved by this fixture')


@lru_cache(None)
def fourier_review():
    a=finite.algebra();e=a.exact();waves,m,_=finite.fixture();waves={k:v for k,v in waves.items() if abs(k)<=2}
    zero=finite.zero_vector;plus=finite.plus_vectors
    scale=lambda v,c:[a.scale(x,c) for x in v]
    add=lambda u,v:{k:plus(u.get(k,zero()),v.get(k,zero())) for k in set(u)|set(v)}
    negate=lambda u:{k:scale(v,-1) for k,v in u.items()}
    square=lambda u:finite.convolution(u,u)
    same=lambda u,v:all(u.get(k,zero())==v.get(k,zero()) for k in set(u)|set(v))
    # Nontrivial complex quadratures and both original input harmonics.
    v={k:scale(waves[k],Q(2,7)) for k in (-1,1)}
    f=scale(m,Q(1,5));sl=scale(m,Q(1,11));dm=plus(f,sl)
    old=add(waves,{0:m});after=add(add(waves,v),{0:plus(m,dm)})
    exact=add(square(after),negate(square(old)))
    cross1=finite.convolution(waves,v);cross2=finite.convolution(v,waves);sq=square(v)
    oldm1=finite.convolution({0:m},v);oldm2=finite.convolution(v,{0:m})
    wcur=add(waves,v)
    feedback1=finite.convolution({0:dm},wcur);feedback2=finite.convolution(wcur,{0:dm})
    rhs={}
    for group in (cross1,cross2,sq,oldm1,oldm2,feedback1,feedback2):rhs=add(rhs,group)
    perp=lambda u:{k:v for k,v in u.items() if k!=0}
    # Mean/mean terms are angular zero and are NOT missing from the full
    # mean equations. They disappear only in this Pperp comparison.
    checks=dict(
        exact_current_signed_then_two_mean_nonzero_difference=same(perp(exact),perp(rhs)),
        zero_mode_not_discarded_from_actual_full_system=not same(exact,rhs),
        both_signed_wave_cross_orientations_required=not same(perp(exact),perp(add(rhs,negate(cross1)))) and not same(perp(exact),perp(add(rhs,negate(cross2)))),
        signed_square_is_not_homogeneous_linear_cancellation=not same(perp(exact),perp(add(rhs,negate(sq)))),
        both_mean_feedback_orientations_required=not same(perp(exact),perp(add(rhs,negate(feedback1)))) and not same(perp(exact),perp(add(rhs,negate(feedback2)))),
        signed_addition_creates_third_harmonics=3 in cross1 or 3 in cross2,
        signed_update_does_not_introduce_fifth_or_higher_modes=max(abs(k) for k in rhs)<=3,
        inherited_fourth_harmonic_still_required=4 in square(waves),
        updated_complete_velocities_have_only_first_and_second_harmonics=set(wcur)=={-2,-1,1,2},
        conjugate_output_reality_preserved=all(exact[-k]==[a.conj(x) for x in v] for k,v in exact.items()),
    )
    # Sensitivity to an actual axial repair velocity (distinct from A1 Fz
    # residual tail) in BOTH orientations of the complete current wave.
    repair=scale(m,Q(1,13));full=plus(dm,scale(repair,-1))
    fullfeedback=add(finite.convolution({0:full},wcur),finite.convolution(wcur,{0:full}))
    checks['deleting_actual_repair_velocity_changes_nonzero_feedback']=not same(fullfeedback,add(feedback1,feedback2))
    return dict(checks=checks,nonzero_update_outputs=sorted(perp(rhs)),
        scope='complete cylindrical Fourier curl polynomials and two generic divergence-free mean increments; not actual canonical signed/mean outputs or a C-infinity support enclosure')


def data_fixture(n):
    keys=('b','G','graph','c','r')+prefix.LINEAR+('T','R','oldU','oldM','fastM','slowM','normal_time','inverse_normal_square')
    return {k:P({0:i+2,1:Q(n+1,9)}) for i,k in enumerate(keys)}


@lru_cache(None)
def budget_review():
    checks={};rows=[]
    for n in (0,1,2):
        data=data_fixture(n);out=regular_source_envelope(n,P(7),data)
        for x in (1,4):
            raw=regular_source_envelope(n,Q(7),{k:v.evaluate(x) for k,v in data.items()})
            checks[f'n{n}_x{x}:full_positive_source_polynomial_substitution']=all(v.evaluate(x)==raw[k] for k,v in out.items() if isinstance(v,P))
        # Independent physical epsilon ledger at kappa=1/10, epsilon=t^10.
        # This retains actual graph, frequency, pressure and curl powers;
        # the positive recipe has factored out the justified order D=t^13.
        v={key:poly.evaluate(1) for key,poly in data.items()}
        for t in (Q(1,2),Q(1,8)):
            ep=t**10;d=n+4;H=2;dr=v['graph']*d/t;r=v['r'];normal=v['normal']
            E=(ep*d+v['b']*ep*dr+v['G']*ep*d+2*H*v['Eik']*t**5
               +(v['bR']+v['Bz']+v['b']*r)*ep
               +ep*(d*d*(v['graph']**2/t**2+1)+r*dr+r*r)
               +4*H*normal*t**5*(dr+ep*d)+2*H*normal*t**5*(dr+ep*d+r)
               +4*H*normal*r*t**5)
            pi=homogeneous_pressure(v['T'],normal=normal,K=v['K'],c=v['c'],normal_time=v['normal_time'],inverse_normal_square=v['inverse_normal_square'])
            linear=v['cutoff']*(E*(v['T']*t**9+v['R']*t**13)+(dr+ep*d)*pi*t**14+(v['c']*d+v['K']+4*H*H*normal**2)*v['R']*t**13)
            old=v['oldU']*t**5;new=v['cutoff']*(v['T']*t**9+v['R']*t**13)
            Cww=3*(dr+ep*d+r);Cmw=2*(dr+ep*d+r)+2*H*normal/t**5
            transport=Cww*(2*old*new+new*new)+Cmw*(v['oldM']*t**9*new+(v['fastM']*t**13+v['slowM']*t**14)*(old+new))
            raw=7*t**13+linear+transport
            checks[f'n{n}_t{t}:independent_physical_complete_source_bound_D']=raw<=t**13*out['F'].evaluate(1)
        checks[f'n{n}:old_regular_source_is_added_not_cancelled']=out['F'].c==(P(7)+out['wave_stage']+out['mean_feedback']).c
        for key,group in (('R','principal_on_uncut_curl'),('normal_time','pressure_gradient'),('fastM','mean_feedback'),('slowM','mean_feedback')):
            d=dict(data);d[key]=P()
            checks[f'n{n}:omission_of_{key}_detected']=out[group].evaluate(1)>regular_source_envelope(n,P(7),d)[group].evaluate(1)
        checks[f'n{n}:actual_cut_velocity_uses_cutoff_once']=out['complete_cut_signed'].c==(data['cutoff']*(data['T']+data['R'])).c
        checks[f'n{n}:full_signed_pressure_derived_not_independent_callback']=out['signed_pressure'].c==homogeneous_pressure(data['T'],normal=data['normal'],K=data['K'],c=data['c'],normal_time=data['normal_time'],inverse_normal_square=data['inverse_normal_square']).c
    for k in (Q(1,100000),Q(1,10),Q(249,1000)):
        r=rates(k)
        checks[str(k)+':all_current_nonzero_terms_have_D_or_better']=r['whole_regular']==r['D'] and min(r['groups'].values())==r['D']
        checks[str(k)+':new_mean_C0_does_not_promote_nonzero_source']=r['regular_mean']==r['C0']>r['D'] and r['groups']['signed_complete_linear']==r['D']
        checks[str(k)+':quantitative_seed_relation']=r['D']==r['seed_relation']
        rows.append({key:str(value) for key,value in r.items() if isinstance(value,Q)})
    # Original root methods are really executed at moderate finite order;
    # giant actual input jets are only referenced in the request recipes.
    root=signed.edge.root.fixture();target=signed.PairedTargetBudget(signed.edge.fixture(),signed.edge.fixture())
    env=signed.RootSignedEnvelope(root,target,lambda n:mean.P(n+2))
    checks['signed_wave_uses_original_root_and_complete_curl_recurrence']=all(env.through(n)['L'].evaluate(1)>0 and env.through(n)['R'].evaluate(1)>0 for n in (0,1,2))
    return dict(checks=checks,rates=rows,scope='finite positive-polynomial and rational exponent diagnostics; full original pulse/input inequalities remain conditional')


@lru_cache(None)
def request_review():
    j=fixture();c=j.certificate(0,0);d=c['domains'];q=Q(d['domain']['q']);checks={}
    a0,a1=c['compact']['A0'],c['compact']['A1']
    checks['physical_counts_exceed_normalized_mean_counts']=a0['raw_source_order']==60200007 and a1['raw_source_order']==60100009 and a0['raw_source_order']>200003
    checks['current_physical_target_request_is_A0_source_plus_14']=d['request']['trace']['maxima']['next_signed_target_input']==60200021
    checks['whole_physical_request_needs_larger_common_pole']=q==361200138 and q>Q(j.mean.common_pole(0)['q'])
    checks['physical_A0_and_A1_use_current_4q_and_2q_sources']=Q(a0['raw_kernel']['source_pole'])==4*q and Q(a1['raw_kernel']['source_pole'])==2*q
    checks['entire_nonzero_source_uses_three_q_domain']=Q(d['domain']['full_nonzero_regular_source'])==3*q>=9*j.mean.original_pole
    checks['separate_oriented_primitive_constant_kept']=d['operators']['oriented_target']['pole']==str(2*q) and d['operators']['common_argument'].startswith('S+K_old+K_mean+K_Ie')
    checks['raw_profiles_source_factorials_and_graph_floor_retained']=all(r['raw_kernel']['source_factorial']==r['raw_source_order'] and len(r['raw_kernel']['factors'])>=4 and not r['raw_kernel']['materialized_actual_input_bound'] for r in (a0,a1))
    checks['current_tail_signs_and_pressure_compensation_kept']=a0['exact_tail'].startswith('-A0') and a1['exact_tail'].startswith('+A1') and 'gr_final' in c['full_identity']
    checks['new_signed_tail_has_no_inhomogeneous_source_term']='f=0' in c['regular']['new_tail']
    checks['base_stress_not_counted_twice']='once' in c['base_exterior']
    checks['operators_and_huge_actual_fields_not_materialized']=not d['operators']['pressure_and_fast']['materialized'] and all(not r['evaluated_actual_high_jet_bound'] for r in (a0,a1))
    checks['regular_grammar_only_not_actual_complete_theorem']=c['complete_regular_grammar_mean_tail_join_is_conditional'] and not c['complete_actual_residual_improvement_absolute_decay_or_blowup_certified']
    for h in (Q(1,200),Q(1,1000),Q(99,10000)):
        for m,N in ((0,0),(1,2),(3,7)):
            r=physical_request(m,N,h)
            checks[f'h{h}_m{m}_N{N}:physical_margins_strict']=all(r[key]['strict_margin']>0 for key in ('A0','A1'))
            checks[f'h{h}_m{m}_N{N}:source_jet_trace_not_low_output_pole']=r['trace']['maxima']['next_signed_target_input']==r['A0']['raw_source_order']+14
        out=physical.regular_transfer(0,h)
        checks[f'h{h}:regular_source_still_does_not_certify_absolute_decay']=Q(out['native_physical_exponent'])==-Q(3,2)-h*(Q(1,2)+Q(2,100000))<0 and not out['absolute_decay_certified']
    # Moderate original operators + oriented primitive can be explicitly
    # dominated by the SUM shift. No monotonicity in pole is assumed.
    geometry=j.mean.geometry;ms=mean.mean_shift(3,1,geometry);K=mean.fresh.materialize_shift(ms['recipe'])
    ip=signed.edge.primitive_recipe(12,4,Q(1,2));Ki=sum(a*3**__import__('math').ceil(b) for a,b in ip['branches'].values())
    checks['moderate_sum_shift_dominates_both_distinct_operator_constants']=K+Ki>=K and K+Ki>=signed.edge.rational_upper(ip)
    oldK=K+Ki+7;total=oldK+K+Ki;poly=P({0:2,1:3,2:1})
    checks['archived_lower_shifts_explicitly_retained_without_pole_monotonicity']='K_old' in d['operators']['archived_estimate_shifts'] and all(poly.evaluate(x+total)>=poly.evaluate(x+oldK) for x in (1,4))
    checks['omitting_archived_shift_need_not_dominate_old_input_polynomial']=poly.evaluate(1+oldK)>poly.evaluate(1+K+Ki)
    try:j.certify_actual_inputs()
    except NotImplementedError:checks['actual_input_certification_is_refused']=True
    else:checks['actual_input_certification_is_refused']=False
    return dict(checks=checks,certificate=c,scope='exact finite integer call graph and original deferred raw kernel recipes; illustrative p=1 and h=1/200 do not certify actual inputs or select h')


@lru_cache(None)
def review():
    parts=dict(pressure=pressure_review(),fourier=fourier_review(),budgets=budget_review(),requests=request_review())
    return dict(parts=parts,checks={name+':'+key:bool(value) for name,part in parts.items() for key,value in part['checks'].items()})
