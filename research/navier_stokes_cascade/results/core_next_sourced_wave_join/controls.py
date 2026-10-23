"""Exact fourth-harmonic forced ODE, cut curl, covariance and norm checks."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import comb,factorial
import sympy as s
from budgets import HERE,load,current,generated,tangent,inverse,P,HARMONICS,rates,FiniteHarmonicWave,wave_stage_envelope,next_domains,NextSourcedRecipes,postwave_mean_envelope


def current_fixture():
    class Common:reweighting=dict(original_pole='1',L='4',width='1/2',gap='1')
    common=Common()
    class Final:means=common
    class Haar:means=common
    class Signed:
        final=Final();haar=Haar()
        def target_inputs(self,n):return dict(stage='same-current-before-generic-signed')
    m=current.mean.PostSignedMeanRecipes(Signed(),wave_pole=30,original_pole=1,geometry=common.reweighting)
    return current.CurrentResidualJoin(m,Q(1,200))


def frame_fixture(trace=None):
    p=tangent.normal.frame.inputs(bmin=1,bmax=1,Fmax=1,gmax=2,gmin=1,u=Q(1,100),cmin=1,cmax=1,Dbox=6,Mz=1)
    sp=tangent.normal.selector.fixture();L=tangent.normal.normal_length_certificate(p)['length_rule_floor']
    sp['length']=Q(-(-L.numerator//L.denominator));ell=tangent.normal.selector.select_tail(sp)['ell_tail']
    kw=dict(F=lambda n:P(n+2),G=lambda n:P(n+3),btilde=lambda n:P(n+1),path=lambda n:P({1:1}) if not n else P(int(n==1)),source=lambda n:P(n+1))
    modes={}
    for m in (1,2,3,4):
        _,f,c=tangent.forced_from_frame(p,sp,ell,harmonic=m,**kw);modes[m]=f
    return FiniteHarmonicWave(modes,lambda n:P(n+1),lambda n:P(n+2),lambda n:P(n+1)),c


@lru_cache(None)
def inverse_review():
    z,x=s.symbols('z x',real=True);n=s.Matrix([1,z,0]);nt=n.diff(z);N=n.dot(n)
    K=s.Matrix([[0,1,0],[-1,0,0],[0,0,1]])
    vel=z*(1+z*x+x*x)*s.Matrix([-z,1,x]);checks={}
    # Manufactured source is NOT an actual canonical residual; the same
    # solution is checked for m=+/-1,...,+/-4 on the entire unit box.
    norm=lambda v:sum(abs(Q(c)) for a in v for c in s.Poly(s.expand(a),z,x).coeffs())
    raw=lambda v,j:max(norm(v.diff(z,i).diff(x,j-i)) for i in range(j+1))
    source={};press={}
    for m in HARMONICS:
        pi=s.I*(z+x+z*x*x)/m
        f=-vel.diff(z)-K*vel-m*m*N*vel-s.I*m*n*pi
        source[m]=f;press[m]=pi
        reconstructed=s.I*((n.dot(K*vel)-nt.dot(vel))+n.dot(f))/(m*N)
        checks[f'm{m}:full_forced_equation_and_constraint']=all(s.expand(a)==0 for a in vel.diff(z)+K*vel+m*m*N*vel+s.I*m*n*pi+f) and s.expand(n.dot(vel))==0
        checks[f'm{m}:full_source_normal_moving_pressure']=s.simplify(pi-reconstructed)==0
        checks[f'm{m}:zero_entrance_and_parameter_entrance_jets']=all(s.expand(v.subs(z,0))==0 for v in vel) and all(s.expand(v.subs(z,0))==0 for v in vel.diff(x))
    F=lambda j:P(max(raw(source[m],j) for m in HARMONICS))
    modes={m:inverse.ForcedBudget(normal=lambda j:P(raw(n,j)),moving_normal=lambda j:P(raw(nt,j)),matrix=lambda j:P(3 if j==0 else 0),
        source=F,profile_rate=lambda j:P(),tangent=P(3),length=P(1),inverse_normal_square=1,harmonic=m) for m in (1,2,3,4)}
    for m,f in modes.items():
        for j in (0,1,2):
            checks[f'm{m}_j{j}:full_box_mixed_velocity_majorant']=raw(vel,j)<=f.velocity(j).evaluate(1)
            # Pressure has the explicit epsilon^1/2 factor, epsilon=k=1
            # in this algebraic diagnostic. Sum absolute complex coefficients.
            pv=s.Matrix([s.im(press[m])])
            checks[f'm{m}_j{j}:full_box_mixed_pressure_majorant']=raw(pv,j)<=f.pressure(j).evaluate(1)
    checks['fourth_mode_generator_retains_sixteen_times_damping']=modes[4].damping(0).c==(16*modes[1].damping(0)).c
    checks['third_and_fourth_modes_use_same_Gaussian_ambient_not_short_residence']=all(modes[m].ambient().c==modes[1].ambient().c for m in (3,4))
    # Exact pure-normal forcing: projected force/velocity vanish while pi
    # has only 1/|m| decay. The ambient normal z'=0 has no high-mode decay.
    g=1+z+x;f=s.Matrix([0,0,g]);normal=s.Matrix([0,0,1])
    for m in (3,4):
        pi=s.I*g/m
        checks[f'm{m}:pure_normal_pressure_not_bounded_from_projected_force']=s.simplify(normal.dot(f)-g)==0 and all(s.expand(a)==0 for a in s.I*m*normal*pi+f) and pi!=0
        checks[f'm{m}:pressure_direct_source_is_not_inverse_square_in_m']=Q(1,m)>Q(1,m*m)
    checks['normal_ambient_cannot_use_high_mode_residence']=Q(1,15)<1 and modes[4].ambient().evaluate(1)>=1
    checks['manufactured_source_is_not_flat_entrance_certificate']=any(s.expand(a.subs(z,0))!=0 for a in source[4])
    return dict(checks=checks,scope='exact moving-normal forced equations and coefficient-sum suprema on a unit box, not selected canonical source/flat support enclosures')


@lru_cache(None)
def composed_review():
    w,c=frame_fixture();checks={};rows=[]
    for n in (0,1,2):
        v=w.addition(n);uncut=w.addition(n,cut=False)
        checks[f'n{n}:all_eight_signed_input_modes_are_summed']=all(v[k].c==sum((w.mode(n,m)[k] for m in HARMONICS),P()).c for k in ('T','R','Pi','U'))
        checks[f'n{n}:full_cut_pressure_and_coefficient_curl_are_positive']=v['Pi'].evaluate(1)>0 and v['R'].evaluate(1)>0 and uncut['Pi'].evaluate(1)>0
        f=w.modes[4]
        direct=generated.through(lambda j:sum((comb(j,i)*w.cutoff_raw(i)*f.pressure(j-i) for i in range(j+1)),P()),n)
        checks[f'n{n}:fourth_mode_cut_pressure_retains_all_cutoff_jets']=w.mode(n,4)['Pi'].c==direct.c
        checks[f'n{n}:full_covariance_cross_and_signed_squares_are_retained']=w.covariance(n,P(2)).c==(4*v['U']+v['U']*v['U']).c
        rows.append(dict(order=n,degree=max(v['U'].c)))
    checks['original_normal_floor_and_frame_certificate_recomputed']=all(c['gates'].values()) and c['original_normal_certificate']['inverse_normal_square']==4 and not c['actual_primitive_enclosures_certified']
    checks['old_single_constructor_omits_actual_third_and_fourth_modes']=set(w.modes)=={1,2,3,4}
    # Coefficient recurrence uses only a finite source order; full curl
    # needs one higher potential jet. Record actual function requests.
    trace=set();fs={m:inverse.fixture(trace=trace,harmonic=m) for m in (1,2,3,4)}
    a=FiniteHarmonicWave(fs,lambda n:P(n+1),lambda n:P(n+1),lambda n:P(n+1));a.addition(2)
    checks['order_two_complete_wave_reads_source_through_three']=max(n for name,n in trace if name=='source')==3
    checks['order_two_complete_wave_reads_moving_normal_through_three']=max(n for name,n in trace if name=='moving_normal')==3
    return dict(checks=checks,rows=rows,scope='actual original inverse/frame norm METHODS on manufactured lower polynomials; source is not identified with the current canonical residual')


@lru_cache(None)
def harmonic_review():
    finite=load('next_source_fourier',HERE.parent/'finite_induction_foundation/finite_checks.py')
    a=finite.algebra();e=a.exact();R,Z,T=e.R,e.Z,e.T;zero=e.Poly();waves={};bump=(R-1)**3*(2-R)**3
    for m in (1,2,3,4):
        f=R**2*bump*(1+Z/(m+7)+T/(m+9))
        C=[(f,f/(m+2)),(R*f,-f/(m+3)),(-R*f,f/(m+5))]
        diff=lambda x,j:(x[0].diff(j),x[1].diff(j))
        sh=lambda x:(x[0].shift(-1),x[1].shift(-1))
        cross=[a.add(sh(C[2]),a.neg(C[1])),C[0],a.neg(sh(C[0]))]
        curl=[a.scale(diff(C[1],1),-e.EPS),a.add(a.scale(diff(C[0],1),e.EPS),a.neg(diff(C[2],0))),(e.radial(C[1][0]),e.radial(C[1][1]))]
        waves[m]=[a.add(x,a.scale(a.imag(y),3*m)) for x,y in zip(curl,cross)];waves[-m]=[a.conj(x) for x in waves[m]]
    checks={};components=0
    for m,v in waves.items():
        div=a.add((e.radial(v[0][0]),e.radial(v[0][1])),a.scale((v[2][0].diff(1),v[2][1].diff(1)),e.EPS))
        normal=a.add((v[1][0].shift(-1),v[1][1].shift(-1)),v[2])
        checks[f'm{m}:complete_curl_divergence']=a.add(div,a.scale(a.imag(normal),3*m))==(zero,zero)
        for k,u in waves.items():
            checks[f'pair{m},{k}:input_ratio_and_zero_output_safe']=all(a.scale(a.imag(a.mul(normal,b)),3*k)==a.scale(a.mul(div,b),-Q(k,m)) for b in u);components+=3
    conv=finite.convolution(waves,waves)
    checks['new_complete_wave_has_nonzero_eighth_harmonic_square']=8 in conv and -8 in conv
    cov=[[a.sumc([a.mul(v[i],waves[-m][j]) for m,v in waves.items()],zero) for j in range(3)] for i in range(3)]
    checks['full_covariance_has_all_six_nonzero_entries']=all(any(bool(c) for c in cov[i][j]) for i in range(3) for j in range(i,3))
    low={m:v for m,v in waves.items() if abs(m)<=2}
    lowcov=[[a.sumc([a.mul(v[i],low[-m][j]) for m,v in low.items()],zero) for j in range(3)] for i in range(3)]
    checks['absent_old_high_modes_still_change_all_six_covariances']=all(cov[i][j]!=lowcov[i][j] for i in range(3) for j in range(i,3))
    hc=current.mult.harmonic_constants(HARMONICS)
    checks['fourth_input_requires_H4_ratio4_and_outputs_through8']=hc['max_abs']==4 and hc['max_input_ratio']==4 and hc['nonzero_outputs']==list(range(-8,0))+list(range(1,9))
    return dict(checks=checks,exact_phase_components=components,outputs=sorted(conv),scope='exact finite-regularity complete cylindrical Fourier curl polynomials; no actual inverse/pulse support enclosure')


@lru_cache(None)
def budget_review():
    checks={};rows=[];cp=current.mult.P
    for n in (0,1,2):
        keys=('b','G','graph','c','r')+current.mult.prefix.LINEAR+('T','R','Pi','oldU','newU','oldM')
        data={k:cp({0:i+2,1:Q(n+1,9)}) for i,k in enumerate(keys)}
        r=wave_stage_envelope(n,data)
        for x in (1,3):
            direct=wave_stage_envelope(n,{k:v.evaluate(x) for k,v in data.items()})
            checks[f'n{n}_x{x}:wave_stage_positive_polynomial_substitution']=all(v.evaluate(x)==direct[k] for k,v in r.items() if isinstance(v,cp))
        low=current.mult.increment_recipe(n,dict(data,dM=cp()),(-2,-1,1,2))
        checks[f'n{n}:old_H2_constants_underestimate_H4_linear_and_transport']=r['groups']['mixed'].evaluate(1)==2*low['groups']['mixed'].evaluate(1) and r['principal_on_uncut_curl'].evaluate(1)>low['principal_on_uncut_curl'].evaluate(1) and r['Cww'].evaluate(1)>low['Cww'].evaluate(1)
        checks[f'n{n}:no_fresh_mean_feedback_yet']=r['mean_feedback'].c=={} and r['F'].c==r['wave_stage'].c
        fixed={k:cp(2) for k in ('graph','r','Kmom','Kp','rho')}
        oldrows={k:cp(3) for k in ('Etheta','Ez')};oldscalars={k:cp(4) for k in ('P','Jtheta','Jz')}
        pm=postwave_mean_envelope(n,oldrows,oldscalars,cp(5),cp(6),fixed)
        checks[f'n{n}:old_pressure_and_scalar_memory_not_reset']=all(pm['scalars'][k].evaluate(1)>pm['scalar_change'][k].evaluate(1)>0 for k in oldscalars)
        checks[f'n{n}:full_graph_pressure_source_and_compensation_retained']=pm['complete_pressure_change']['pm'].evaluate(1)>pm['complete_pressure_change']['gr'].evaluate(1)>0
    for k in (Q(1,100000),Q(1,10),Q(249,1000)):
        r=rates(k)
        checks[str(k)+':sourced_D_velocity_not_promoted_to_C0_pressure_order']=r['complete_new_velocity']==r['D']<r['new_pressure_and_potential']==r['C0']
        checks[str(k)+':wave_stage_regular_gain_is_half_minus_kappa']=r['wave_stage_regular']==r['E']==2-3*k and r['E']-r['D']==Q(1,2)-k
        checks[str(k)+':full_tangential_E_but_Haar_target_and_scalars_C0']=r['postwave_full_tangential']==r['E']<r['postwave_Haar_tangential']==r['postwave_scalar']==r['C0']
        rows.append({name:str(value) for name,value in r.items() if isinstance(value,Q)})
    return dict(checks=checks,rates=rows,scope='positive norm and exact exponent bookkeeping; wave-stage gain only, no later fresh means or complete-cycle claim')


@lru_cache(None)
def recipe_review():
    j=NextSourcedRecipes(current_fixture());c=j.certificate(0);d=c['domains'];q=Q(d['common_pole']['q']);checks={}
    checks['same_current_source_program_not_new_fine_source_callback']=c['source']['method']=='current.regular_source_envelope' and c['source']['current_stage_recipe']['method']=='regular_source_envelope' and not c['source']['additive_compact_Gaussian_or_base_residual_forwarded']
    checks['complete_wave_source_jet_and_full_row_jet_ladders']=c['wave']['source_order']==4 and d['current_source_through']==4 and d['new_inverse_normal_through']==5
    checks['padded_old_graph_covers_all_new_requests_without_collapsing_old_DAGs']=q==1200138 and d['previous_complete_trace']['maxima']['next_signed_target_input']==200021 and not d['previous_complete_trace']['old_dependency_DAGs_expanded']
    checks['new_6q_pressure_domain_not_old_4q']=Q(d['ledger']['full_covariance_and_pressure_source'])==6*q>4*q and not d['ledger']['old_2q_4q_proof_automatically_covers_new_6q']
    checks['new_6q_12q_original_kernel_recipes_and_Ie6q_retained']=d['new_operators']['parameter']==str(q) and d['new_oriented_target_primitive']['pole']==str(6*q) and d['ledger']['required_operators']['Kp']==[str(6*q),str(12*q)]
    checks['all_archived_shifts_and_separate_oriented_constants_kept']='K_old_all' in d['common_argument'] and d['old_oriented_target_primitive']['pole']==str(2*q)
    checks['no_new_actual_band_fields_or_huge_constants_materialized']=not d['new_band_physical_field_or_geometry_changed'] and not d['new_operators']['materialized'] and not c['wave']['actual_source_coefficients_or_fields_instantiated']
    checks['new_source_tail_uses_CURRENT_regular_source_only']='CURRENT_Gregular' in c['exact_Gaussian_ledger'] and 'source tails excluded' in c['exact_Gaussian_ledger']
    checks['recomputed_pressure_cannot_reuse_stale_A0_residual']='do not reattach the stale' in c['compact_pressure']
    checks['new_pointwise_radial_A0_not_claimed_from_tangential_pressure_size']=not c['postwave_pointwise_radial_A0_bound_joined']
    checks['whole_path_and_flat_entrance_support_still_required']=c['original_whole_path_support_and_flat_entrance_contracts_required']
    checks['no_fresh_mean_complete_cycle_or_actual_validity_certificate']=not c['next_fresh_means_signed_cycle_or_actual_input_validity_certified'] and not c['absolute_decay_infinite_field_or_blowup_certified']
    # Moderate domain kernel proof is executed, not only a declared ledger.
    new=current.mean.fresh.shift_recipe(6,4,Q(1,2),1)
    checks['moderate_original_shift_covers_new_6q_12q_pressure']=new['covered_operator_poles']['Kp']==['36','72'] and current.mean.fresh.materialize_shift(new)>0
    try:j.materialize_actual_source()
    except NotImplementedError:checks['actual_high_jet_source_materialization_refused']=True
    else:checks['actual_high_jet_source_materialization_refused']=False
    return dict(checks=checks,certificate=c,scope='exact finite/deferred current-source routing, not actual coefficient/profile/smooth-support validity')


@lru_cache(None)
def review():
    parts=dict(inverse=inverse_review(),composed=composed_review(),harmonics=harmonic_review(),budgets=budget_review(),recipes=recipe_review())
    return dict(parts=parts,checks={name+':'+key:bool(v) for name,p in parts.items() for key,v in p['checks'].items()})
