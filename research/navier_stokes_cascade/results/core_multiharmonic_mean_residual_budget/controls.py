"""Exact multiharmonic transport controls and a sharp mean-feedback witness."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from budgets import HERE,load,join,prefix,P,convert,hull,harmonic_constants,increment_recipe,rates,source_pole,NextResidualBudget

finite=load('multi_residual_exact_fourier',HERE.parent/'finite_induction_foundation/finite_checks.py')


def norm(poly):return sum((abs(c)*max(Q(1),Q(2)**r) for (r,z,t,u),c in poly.terms.items()),Q(0))
def cnorm(c):return norm(c[0])+norm(c[1])
def vnorm(v):return sum(map(cnorm,v),Q(0))
def jetnorm(v):
    return vnorm(v)+sum(vnorm([(r.diff(i),s.diff(i)) for r,s in v]) for i in range(4))


def harmonic_review():
    s=finite.algebra();e=s.exact();allwaves,_,_=finite.fixture()
    waves={m:a for m,a in allwaves.items() if abs(m)<=2};checks={};count=0
    for m,a in waves.items():
        diva=s.add((e.radial(a[0][0]),e.radial(a[0][1])),s.scale((a[2][0].diff(1),a[2][1].diff(1)),e.EPS))
        normal=s.add((a[1][0].shift(-1),a[1][1].shift(-1)),a[2])
        checks[f'input{m}:complete_curl_divergence']=s.add(diva,s.scale(s.imag(normal),3*m))==(e.Poly(),e.Poly())
        for n,b in waves.items():
            for j,bi in enumerate(b):
                phase=s.scale(s.imag(s.mul(normal,bi)),3*n)
                reduced=s.scale(s.mul(diva,bi),-Q(n,m))
                checks[f'pair{m},{n}:component{j}:exact_input_ratio']=phase==reduced;count+=1
    # m=1,n=-2 isolates the factor missing from a +/-1-only recipe.
    a,b=waves[1],waves[-2]
    diva=s.add((e.radial(a[0][0]),e.radial(a[0][1])),s.scale((a[2][0].diff(1),a[2][1].diff(1)),e.EPS))
    unit=s.mul(diva,b[0]);twice=s.scale(unit,2)
    checks['single_frequency_ratio_one_underestimates_a_nonzero_pair']=cnorm(twice)==2*cnorm(unit)>cnorm(unit)>0
    actual=finite.convolution(waves,waves)
    U=sum(map(jetnorm,waves.values()),Q(0));Cww=3*(4+4+1)
    checks['independent_full_fourier_square_has_l1_majorant']=sum(map(vnorm,actual.values()),Q(0))<=Cww*U*U
    checks['zero_and_highest_fourth_harmonic_are_retained']=0 in actual and 4 in actual and -4 in actual
    checks['conjugate_complete_transport_reality']=all(actual[-m]==[s.conj(x) for x in a] for m,a in actual.items())
    hc=harmonic_constants(waves)
    checks['all_three_output_sets_use_inputs_not_output_division']=hc['max_abs']==2 and hc['max_input_ratio']==2 and hc['quadratic_outputs']==list(range(-4,5))
    return dict(checks=checks,exact_phase_components=count,outputs=sorted(actual),scope='exact compact finite-regularity cylindrical Fourier polynomials; no actual pulse enclosure')


def witness_review():
    s=finite.algebra();e=s.exact();R=e.R;z=e.Poly();cp=lambda a:(a,z)
    psi=(R-1)**4*(2-R)**4;f=e.radial(psi);a=psi
    mean=[cp(z),cp(z),cp(f)];wave=[cp(z),cp(R*a),cp(-a)]
    forward=finite.transport(mean,wave,1);reverse=finite.transport(wave,mean,0)
    expected=[(z,z),s.imag(cp(3*f*R*a)),s.imag(cp(-3*f*a))]
    # n=(0,1/R,1): both the mean and this wave are exactly divergence free.
    normal=s.add((wave[1][0].shift(-1),wave[1][1].shift(-1)),wave[2])
    checks=dict(mean_feedback_is_exact_nonzero_phase_advection=forward==expected and vnorm(expected)>0,
        reverse_mean_transport_is_zero_in_witness=reverse==[(z,z)]*3,
        wave_normal_constraint_and_full_divergence_exact=normal==(z,z),
        mean_is_divergence_free_and_has_both_zero_velocity_moments=(R*f).integral()==0 and mean[0]==mean[1]==(z,z),
        compact_potential_endpoints_preserve_axis_support=all(psi.at(0,r)==0 for r in (1,2)),
        f_and_wave_overlap_nontrivially=bool(f*a))
    r=rates();delta=r['desired']-r['B'];rows=[]
    for j in (8,16,24,32):
        # epsilon=2^-2^j, k=2^2^(j-1): epsilon^1/2*k=1.
        rows.append(dict(j=j,log2_normalized_ratio=str(delta*2**j),
                         log2_after_degree4_S_cost=str(delta*2**j-8*j)))
    checks['positive_gap_exponential_cannot_be_absorbed_by_fixed_polynomial']=delta>0 and all(Q(x['log2_after_degree4_S_cost'])>0 for x in rows)
    return dict(checks=checks,diagnostic_rows=rows,scope='admissible local mean/wave class witness, NOT the prescribed actual fast/slow outputs; smooth plateau extension is analytical')


def data_fixture(n):
    keys=('b','G','graph','c','r')+prefix.LINEAR+('T','R','Pi','oldU','newU','oldM','dM')
    return {k:P({0:i+1,1:Q(n+1,8)}) for i,k in enumerate(keys)}


def recipe_review():
    checks={};rows=[]
    for n in (0,1,2):
        d=data_fixture(n);out=increment_recipe(n,d)
        for x in (1,3):
            exact=increment_recipe(n,{k:v.evaluate(x) for k,v in d.items()})
            checks[f'order{n}:S{x}:positive_polynomial_rational_substitution']=all(v.evaluate(x)==exact[k] for k,v in out.items() if isinstance(v,P))
        oldkeys=('b','G','graph','c','r')+prefix.LINEAR+('T','R','Pi')
        old={k:d[k] for k in oldkeys};old.update(U=P(),M=P())
        checks[f'order{n}:single_frequency_linear_recipe_matches_historical_operator']=increment_recipe(n,d,(-1,1))['linear'].c==prefix.source_recipe(n,old)['linear'].c
        for key,group in (('Pi','pressure_gradient'),('R','principal_on_uncut_curl'),('dM','mean_feedback')):
            bad=dict(d);bad[key]=P()
            checks[f'order{n}:omission_of_{key}_detected']=out[group].evaluate(1)>increment_recipe(n,bad)[group].evaluate(1)
        one=increment_recipe(n,d,(-1,1))
        checks[f'order{n}:second_mode_increases_mixed_and_damping_costs']=out['groups']['mixed'].evaluate(1)==2*one['groups']['mixed'].evaluate(1) and out['principal_on_uncut_curl'].evaluate(1)>one['principal_on_uncut_curl'].evaluate(1)
        rows.append(dict(order=n,degree=max(out['F'].c)))
    r=rates();fine=rates(fast_order=r['desired'],slow_order=r['C'])
    checks['coarse_fresh_means_do_not_supply_improved_residual_order']=r['guaranteed']==r['B']<r['desired']
    checks['desired_fast_and_slow_input_orders_are_sufficient_only']=fine['guaranteed']==r['desired'] and fine['sufficient_slow_margin']>0
    checks['new_mean_feedback_sets_source_pole_nine_p']=source_pole(Q(3,2))['source']==Q(27,2)
    return dict(checks=checks,rows=rows,rates={k:str(v) for k,v in r.items() if isinstance(v,Q)},scope='exact positive-polynomial and rational exponent bookkeeping')


@lru_cache(None)
def composed_review():
    # Reuse the exact parent fixture constructor while routing its imports
    # to this module's OWN consistently loaded polynomial classes.
    path=HERE.parent/'core_reweighted_fresh_mean_closure/controls.py'
    import sys
    current=sys.modules.get('budgets');sys.modules['budgets']=__import__('budgets').fresh
    try:c=load('multi_residual_parent_fixture',path)
    finally:
        if current is not None:sys.modules['budgets']=current
        else:del sys.modules['budgets']
    trace=set();means=c.fixture(trace);b=NextResidualBudget(means)
    r=b.regular_source(0);u=b.uncut_increment(0);g=b.generated
    raw=b.raw_source(0);converted=b.inverse_source(join.tangent.transfer.inverse)(0)
    high=max(n for name,n in trace if name=='initial')
    checks=dict(full_current_pressure_and_all_scalars_retained=all(r['full_current_mean_rows'][k].evaluate(1)>0 for k in ('gr','P','pm','Jtheta','Jz')),
        raw_regular_source_converts_exactly_to_inverse_interface=raw.c==r['F'].c==converted.c,
        mean_feedback_not_omitted_or_promoted_to_wave_stage=r['F'].evaluate(1)>r['wave_stage'].evaluate(1)>0 and r['mean_feedback'].evaluate(1)>0,
        full_unprojected_source_pressure_and_uncut_curl_positive=u['Pi'].evaluate(1)>0 and u['R'].evaluate(1)>0,
        uncut_principal_bound_does_not_contain_temporal_cutoff=u['T'].c==hull(sum((convert(join.through(g.modes[abs(m)].velocity,0)) for m in g.harmonics),P())).c,
        finite_original_mean_jet_loss_is_31=high==31,
        same_parent_shift_width_and_band_preserved=means.reweighting['width']=='1/2' and 'no physical scale' in means.reweighting['bound_argument'],
        actual_primitive_enclosures_not_promoted=not g.certificate['actual_primitive_enclosures_certified'],
        inverse_interface_explicitly_carries_only_coarse_order_and_nine_pole=r['order_certificate']==dict(guaranteed=rates()['B'],desired=rates()['desired'],improved_order_certified=False,source_pole=Q(9)))
    return dict(checks=checks,highest_initial_jet=high,source_degree=max(raw.c),scope='manufactured retained-input composition; original input validity and improved mean feedback remain uncertified')


@lru_cache(None)
def review():
    parts=dict(harmonics=harmonic_review(),sharp_mean_feedback=witness_review(),polynomial_recipe=recipe_review(),composed_prefix=composed_review())
    return dict(checks={name+':'+k:v for name,row in parts.items() for k,v in row['checks'].items()},parts=parts)
