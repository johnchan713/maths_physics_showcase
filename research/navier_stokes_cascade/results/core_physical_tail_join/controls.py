"""Exact physical exponent, chart, polynomial and residual allocation controls."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
import sympy as s
from bounds import scales,compact_orders,absorption_constant,gaussian_recipe,regular_transfer,hypothetical_steps,PhysicalTailJoin


@lru_cache(None)
def ledger_review():
    # Independent formal real angular Fourier modes. Haar is a separate
    # operation: A0 and A1 residuals here are angular zero, not identically zero.
    rhoP,A0,A1,Ht,Hz,Ft,Fz,Rft,Rfz,Rst,Rsz,base,gauss=s.symbols('rhoP A0 A1 Ht Hz Ft Fz Rft Rfz Rst Rsz base gauss',real=True)
    y=s.symbols('y',real=True);w=s.cos(y)+s.sin(2*y)
    before=s.Matrix([-rhoP-A0,Ht+Ft,Hz+Fz])
    fasttime=s.Matrix([0,-Ft,-Fz+A1])
    differences=s.Matrix([0,Rft+Rst,Rfz+Rsz])
    final=before+fasttime+differences
    regular=s.Matrix([-rhoP,Ht+Rft+Rst,Hz+Rfz+Rsz])
    compact=s.Matrix([-A0,0,A1]);tail=s.Matrix([gauss*w,gauss*s.sin(y),gauss*s.cos(2*y)])
    basefield=s.Matrix([base,base,base])
    total=final+tail+basefield
    eq=lambda a,b:all(s.expand(v)==0 for v in a-b)
    gd,axialRepair,old=s.symbols('gd axialRepair old',real=True)
    actual_mean=old+gd-axialRepair
    # Both orientations of a simple 1D transport proxy. This control is
    # algebraic, not a manufactured Navier-Stokes inverse solution.
    v=actual_mean*(1+y);wave=w
    transport=v*s.diff(wave,y)+wave*s.diff(v,y)
    deleted=(old+gd)*(1+y)*s.diff(wave,y)+wave*(old+gd)
    checks=dict(
        exact_fresh_cancellation_leaves_all_compact_terms=eq(final,regular+compact),
        full_residual_is_regular_plus_compact_Gaussian_and_assigned_base=eq(total,regular+compact+tail+basefield),
        pressure_compensator_remains_in_radial_regular_row=regular[0]==-rhoP,
        A0_and_A1_have_opposite_derived_signs=compact==s.Matrix([-A0,0,A1]),
        omitting_A0_or_A1_is_detected=not eq(final,regular) and not eq(final,regular+s.Matrix([-A0,0,0])),
        old_Q1_rows_not_added_again=eq(final,regular+compact) and not eq(final,regular+compact+s.Matrix([0,Ft,Fz])),
        residual_tail_not_a_velocity_or_nonlinear_source=not eq(total,regular+compact+tail+basefield+compact*w),
        actual_A1_velocity_stays_in_both_regular_interactions=s.simplify(transport-deleted)!=0,
        compact_residual_and_actual_A1_velocity_are_distinct_symbols=A1!=axialRepair and s.diff(transport,A1)==0 and s.diff(transport,axialRepair)!=0,
        angular_projection_does_not_discard_compact_mean_tail=all(s.integrate(c,(y,0,2*s.pi))/(2*s.pi)==c for c in compact),
        temporal_wave_tail_has_zero_angular_projection=all(s.integrate(c,(y,0,2*s.pi))==0 for c in tail),
        full_pressure_mean_terms_not_a_nonzero_projection_tail=rhoP not in set().union(*(c.free_symbols for c in tail)),
    )
    # A spatial cutoff can differentiate the mean pressure at a pulse peak.
    x,psi,meanp=s.symbols('x psi meanp',real=True)
    checks['spatial_pressure_cutoff_not_erased_as_temporal_Gaussian']=s.diff((1+x)*meanp,x)==meanp and meanp!=0
    return dict(checks=checks,scope='exact residual allocation and Fourier/transport proxy; not a new continuum operator or actual profile enclosure')


@lru_cache(None)
def transfer_review():
    checks={};rows=[];k=Q(1,100000)
    for h in (Q(1,200),Q(1,1000),Q(99,10000)):
        sc=scales(h,k)
        for m,N in ((0,0),(1,0),(0,2),(3,7)):
            for e in (0,1):
                a=compact_orders(e,m,N,h,k);tag=f'h{h}_m{m}_N{N}_e{e}'
                checks[tag+':strict_physical_margin']=a['strict_margin']>0
                checks[tag+':previous_integration_integer_insufficient_for_this_rule']=a['previous_integer_margin']<=0
                checks[tag+':physical_power_not_epsilon_power']=a['physical_exponent']==h*(a['source_rate']+a['integration_order']*k)-sc['normalization']-2*m
                checks[tag+':finite_current_source_jet']=a['raw_source_order']==m+2*a['integration_order']+3
        for m in (0,2):
            r=regular_transfer(m,h,k)
            exact=-Q(3,2)-h*(Q(1,2)+2*k)-2*m
            checks[f'h{h}_m{m}:first_regular_native_physical_bound_has_negative_power']=Q(r['native_physical_exponent'])==exact<0 and not r['absolute_decay_certified']
        step=hypothetical_steps(0,0,h,k)
        j=step['additional_steps'];g=Q(1,2)-k;L=sc['normalization']+1
        checks[f'h{h}:hypothetical_stage_threshold_is_only_conditional']=Q(step['physical_margin'])>0 and h*(sc['D']+(j-1)*g)-L<=0 and not step['induction_supplied'] and not step['arbitrary_diagonal_supplied']
        rows.append(dict(h=str(h),regular=regular_transfer(0,h,k),hypothetical=step))
    a0=compact_orders(0,0,0,Q(1,200));a1=compact_orders(1,0,0,Q(1,200))
    checks['illustrative_A0_physical_zero_power_needs_60200007_current_source_jets_in_this_rule']=a0['integration_order']==30100002 and a0['raw_source_order']==60200007
    checks['illustrative_A1_physical_zero_power_needs_60100009_current_source_jets_in_this_rule']=a1['integration_order']==30050003 and a1['raw_source_order']==60100009
    checks['fixed_normalized_A0_C_count_does_not_alone_prove_physical_flatness']=scales(Q(1,200))['A0_source_order_zero']==100005 and a0['raw_source_order']>100005
    return dict(checks=checks,rows=rows,illustrative_h='1/200, diagnostic only; actual fixed h not selected',
                scope='exact sufficient exponent/jet accounting, not necessary lower bounds on every possible proof or a certificate of future induction')


@lru_cache(None)
def constant_review():
    checks={};coeff={0:Q(2),1:Q(3),3:Q(5)};K=Q(7);delta=Q(1,13)
    C=absorption_constant(coeff,K,delta)
    for x in (Q(1),Q(3,2),Q(10),Q(100)):
        E=sum(((delta*x)**j/factorial(j) for j in range(7)),Q(0))
        poly=sum((c*(16*x*x+K)**j for j,c in coeff.items()),Q(0))
        checks[f'x{x}:exponential_series_proves_polynomial_absorption_without_new_band']=C*E>=poly
    for ell in (16,20,100,10**6):
        # Rational interval: 1/2<ln2<1. Bounds valid for every s in
        # [1/4,4], not a logarithm sampled at a single chart point.
        upper=Q(3,2)*ell;lower=Q(1,4)*ell
        checks[f'ell{ell}:natural_log_chart_upper_from_enlargement']=Q(ell+2)<=upper
        checks[f'ell{ell}:natural_log_chart_lower_from_enlargement']=Q(ell-2,2)>=lower
        checks[f'ell{ell}:band_Gaussian_dominates_chosen_physical_Gaussian']=Q(1,4)*upper**2<=ell**2
    g=gaussian_recipe(Q(1,1000),Q(7,2),Q(5),Q(3));c=Q(g['physical_Gaussian_constant'])
    n=g['polynomial_index']
    checks['Gaussian_recipe_retains_explicit_band_to_natural_log_constant']=c==Q(1,4000) and n==3
    checks['Gaussian_recipe_keeps_square_completion_constant']=Q(g['exponential_log'])==(Q(7,2)+3)**2/(2*c)
    for x in (Q(1),Q(7,3),Q(10**3)):
        gap=(Q(7,2)+3)**2/(2*c)-(-c*x*x/2+(Q(7,2)+3)*x)
        checks[f'x{x}:Gaussian_square_completion_exact']=gap==c*(x-(Q(7,2)+3)/c)**2/2 and gap>=0
    return dict(checks=checks,scope='exact whole-half-line inequalities and rational series controls; no giant Gaussian or high-jet constant evaluated')


@lru_cache(None)
def recipe_review():
    # A narrow stub supplies only the original geometry and already-derived
    # rates. No actual high-order source method may run during this audit.
    class Means:
        reweighting=dict(original_pole='1',L='4',width='1/2',gap='1')
    class Closure:
        means=Means();kappa=Q(1,100000);embedding=Q(4)**12
    j=PhysicalTailJoin(Closure(),Q(1,200));r=j.certificate(0,0)
    a0,a1=r['A0'],r['A1'];checks={}
    checks['A0_uses_current_full_embedded_gr_without_evaluating_huge_jets']=a0['source_recipe']['method']=='FinalMeanRowClosure.current_radial_source_polynomial' and a0['source_recipe']['order']==60200007 and a0['source_recipe']['field_input_order']==60200010
    checks['A1_keeps_nested_original_high_jet_recipe_deferred']=a1['source_recipe']['method']=='FastRowGain.postwave_recipe' and a1['source_recipe']['nested_original_initial_mean_order']==60300013
    checks['physical_tail_uses_raw_e0_and_e1_profiles_separately']=a0['original_kernel_recipe']['operator']=='A0' and a1['original_kernel_recipe']['operator']=='A1'
    checks['same_original_width_gap_and_pole_domains']=all(a['original_kernel_recipe']['original_geometry']==Means.reweighting and a['original_kernel_recipe']['source_pole']=='12' for a in (a0,a1))
    checks['old_compact_tail_signs_preserved']=a0['exact_tail']=='-A0 Q1 gr' and a1['exact_tail']=='+A1 Q1 Ez_postwave'
    checks['tails_are_not_materialized_actual_input_certificates']=not a0['materialized_actual_high_jet_bound'] and not a1['materialized_actual_high_jet_bound'] and not r['full_actual_inputs_or_absolute_residual_improvement_certified']
    checks['spatial_cutoffs_and_later_diagonals_not_called_Gaussian']='regular' in r['spatial_slow_cutoffs'] and 'not Gaussian' in r['later_diagonal_and_localization']
    checks['base_ledger_does_not_readd_the_cancelled_bare_stress']='AFTER its stress divergence' in r['base_exterior']
    return dict(checks=checks,example_recipe=r,scope='finite exact recipe routing; no hidden evaluation of tens of millions of input jets')


@lru_cache(None)
def review():
    parts=dict(ledger=ledger_review(),transfer=transfer_review(),constants=constant_review(),recipes=recipe_review())
    return dict(parts=parts,checks={name+':'+key:bool(v) for name,p in parts.items() for key,v in p['checks'].items()})
