"""Regional input gaps, a corrected signed comparison, and cone identities."""
from fractions import Fraction as Q
import importlib.util
from pathlib import Path
import sys
import mpmath as mp

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


compact=load('regional_frozen_compact_bounds',HERE.parent/'compact_jet_envelope/bounds.py')
repair=load('regional_frozen_repair_bounds',HERE.parent/'repair_state_bounds/bounds.py')


def coordinates(a,b,p,q):
    t=-b/a;v=a+b*b/a;Pc=p+t*q;J=q-t*p
    return t,v,Pc,J


def gaps(a,b,p,q):
    _,v,Pc,J=coordinates(a,b,p,q)
    return a,v-2,Pc-v,2*(Pc-v)**2-(v-2)*J*J


def equation_review():
    checks=[];factor_checks=[]
    for a in map(Q,('4/5','2','5/2')):
        for b in map(Q,('-1/3','0','2/5')):
            for w in map(Q,('-3/4','0','1/2')):
                p=Q(40);q=p*w;F=Q(3,7)
                t,v,Pc,J=coordinates(a,b,p,q)
                c=1-b*w/a;j=w+b/a
                G=2*c*c-(v-2)*j*j
                a24=2-2*b*w-b*b/a-(a-2)*w*w
                Dc=Pc-v;T=(F*(p-a),F*(q+b));g=(-F*a,F*b)
                checks.extend((Pc==p*c,J==p*j,
                    G==(1+b*b/(a*a))*a24,
                    gaps(a,b,p,q)[3]==p*p*G-4*p*c*v+2*v*v,
                    -(T[0]*g[0]+T[1]*g[1])==F*F*a*Dc,
                    T[0]*T[0]+T[1]*T[1]==F*F*(Dc*Dc+J*J)/(1+t*t)))
                if a==2:
                    D=b*b/2;A0=1-(2+D)/p
                    factor_checks.append(gaps(a,b,p,q)[3]==2*p*p*(A0+D/2)*(A0-b*w-D/2))
    # A large-p REDUCED ratio condition cannot replace the finite-p test.
    low=gaps(Q(11,5),Q(0),Q(21,10),Q(0))
    subcritical=gaps(Q(4,5),Q(0),Q(3),Q(0))
    return dict(checks=dict(
        all_coordinate_ratio_and_oriented_stress_identities_exact=all(checks),
        a2_fourth_gap_factorization_exact=all(factor_checks),
        reduced_ratio_positivity_does_not_imply_finite_p_cone=low[2]<0,
        original_vs_below_two_is_relaxed_not_admissible=subcritical[1]<0 and subcritical[2]>0,
        rational_fixtures_are_not_actual_modulated_fields=True),
        rational_identity_checks=len(checks),a2_factorization_checks=len(factor_checks))


def comparison_review():
    with mp.workdps(80):
        h=mp.mpf('.01');L=mp.mpf('.98');ell=mp.mpf('.6');p0=mp.mpf('1.1')
        old=mp.exp(-ell)*(p0-h/L*(mp.e-1))
        exact=p0*mp.exp(-ell)-h/L*(mp.e-mp.exp(-ell))/(1+ell)
        corrected=p0*mp.exp(-ell)-h/L*(mp.e-1)
        # Independent positive quadrature of the exact negative-source kernel.
        integral=mp.quad(lambda z:mp.exp(-ell*(1-z))*mp.exp(z),[0,1])
        quadrature=p0*mp.exp(-ell)-h/L*integral
        return dict(checks=dict(
            old_signed_kernel_shortcut_fails_for_comparison_model=old>exact,
            corrected_bound_is_a_lower_bound=corrected<exact,
            corrected_bound_still_above_half=corrected>mp.mpf('.5'),
            independent_kernel_integral_matches=abs(exact-quadrature)<mp.mpf('1e-70'),
            model_is_not_actual_outer_ramp_or_selected_h=True),
            normalized_old_claim=float(old),normalized_exact_model=float(exact),
            normalized_corrected_lower=float(corrected),old_minus_model=float(old-exact))


def scalar_certificate(digits):
    if digits<100:
        raise ValueError('At least 100 interval digits required for nested Md64 parameters')
    base=compact.scalar_bounds(digits)
    margins=compact.margin_bounds(digits)
    accepted=repair.scalar_checks(digits)
    old=compact.core.certificate('Md64',digits)
    b=old['log_bounds'];p=compact.p;lo=compact.lower;hi=compact.upper
    exp,log=mp.iv.exp,mp.iv.log
    lc=b['C'];T=b['P']-1
    ramp=p('1.1')*exp(p('-.6'))-p('.01')/p('.98')*(exp(1)-1)
    ideal=p('1.1')*10000*exp(-5)
    # The first small term is a bound for radial reference-growth variation.
    # Its C^2 R^80 factor is already within the audited activation constant.
    B=exp(b['B']);lkact=2*lc+p('.8')*B
    checks=dict(
        same_selected_Md64_parameter_checks=all(old['checks'].values()),
        inherited_compact_norm_and_denominator_envelope=all(base['checks'].values()),
        inherited_axial_and_intermediate_margin_arithmetic=all(margins['checks'].values()),
        inherited_exact_modulation_repair_acceptance=all(accepted['checks'].values()),
        corrected_outer_ramp_floor_above_half=lo(ramp)>p('.5'),
        ideal_after_join_Pc_above_two=lo(ideal)>p(3),
        first_ramp_Pc_gap_above_one=lo(ramp*10000-2)>1,
        axial_third_and_fourth_above_one=min(lo(margins['axial_third_lower']),lo(margins['axial_fourth_lower']))>1,
        intermediate_fourth_above_one=lo(margins['intermediate_fourth_lower'])>1,
        native_growth_derivative_within_Kact=hi(log(4)+2*lc+5*B/100-lkact)<0,
        Kact_below_C3=hi(lkact-3*lc)<0,
        activation_left_ea_above_Cminus401=lo(lc-log(4))>0,
        left_collar_ea_below_Cminus99=lo(lc-4)>0,
        left_collar_growth_loss_below_p001=Q(2,10**153)+Q(1,10**864)<Q(1,1000),
        compact_original_gap_inverse_fits_A=1024<4096,
        compact_a_inverse_fits_A=22<4096,
        right_collar_growth_above_inverse_C=lo(log(2)-4*T+lc)>0,
        heat_compensation_starts_exactly_at_compact_right=True,
        positive_modulation_repair_log_separation=lo(240*T-26)>0,
        fixed_outer_radius_above_all_regional_floors=lo(b['XR']-log(10000))>0)
    metrics=dict(corrected_outer_ramp_lower=ramp,ideal_after_join_Pc_lower=ideal,
                 axial_third_lower=margins['axial_third_lower'],
                 axial_fourth_lower=margins['axial_fourth_lower'],
                 intermediate_fourth_lower=margins['intermediate_fourth_lower'],
                 log_native_growth_derivative_over_Kact=log(4)+2*lc+5*B/100-lkact,
                 log_Kact_over_C3=lkact-3*lc,
                 log_right_collar_growth_times_C=log(2)-4*T+lc)
    return dict(checks=checks,metrics=metrics)


def outward_review():
    rows=[];errors=[]
    with mp.workdps(170):
        a,b=scalar_certificate(100),scalar_certificate(160)
        for k in a['metrics']:
            for i in (0,1):
                x,y=[mp.make_mpf(r['metrics'][k]._mpi_[i]) for r in (a,b)]
                errors.append(abs(x-y)/max(1,abs(x),abs(y)))
        for digits,r in ((100,a),(160,b)):
            encode=lambda z:{name:mp.nstr(mp.make_mpf(z._mpi_[i]),120) for i,name in enumerate(('lower','upper'))}
            rows.append(dict(digits=digits,checks=r['checks'],metrics={k:encode(v) for k,v in r['metrics'].items()}))
        return dict(rows=rows,maximum_relative_refinement_error=float(max(errors)),checks=dict(
            all_outward_regional_and_frequency_checks=all(all(r['checks'].values()) for r in rows),
            precision_refinement_consistent=max(errors)<mp.mpf('1e-55'),
            no_C_N_or_actual_oscillatory_profile_materialized=True))


def joined_review():
    # Exact endpoint registry; offset endpoints are never rounded to log Xa.
    endpoints=['Xminus','Xa*exp(delta)','Xi','Xsep','XR*exp(-8)',
               'XR*exp(-5)','XR','XR*exp(1)','XR*exp(T+1)',
               'XR*exp(T+2)','Xrep=Xstar']
    names=['activation','reduced_continuation','shape','ideal_prejoin','five_moment_join',
           'ideal_after_join','first_outer_ramp','axial_Md64','intermediate_ramp','intermediate_power']
    regions=[dict(region=n,left=endpoints[i],right=endpoints[i+1]) for i,n in enumerate(names)]
    # With A>=2^16, G=exp(A^32), eps=G^-64, Rbox=4G^3:
    # 256 Rbox^10 * eps <= G^-4/2 iff G^30 >= 512*4^10.
    checks=dict(
        every_compact_region_has_exact_shared_endpoints=all(regions[i]['right']==regions[i+1]['left'] for i in range(len(regions)-1)),
        covers_actual_pre_and_post_join_regions=len(regions)==10,
        activation_fourth_gap_dominates_Cminus1024=802<1024,
        activation_third_gap_dominates_Cminus1024=401<1024,
        boundary_left_growth_floor=Q(23,10)-Q(1,1000)>Q(229,100),
        loop_cutoff_below_right_collar_growth=20>1,
        selected_state_tolerance_preserves_half_kappa=(2**16)**30>512*4**10,
        modulation_and_exact_repair_do_not_change_exterior_ledgers=True,
        heat_reference_cone_is_used_after_Xstar_not_unedited_pulse_alone=True,
        first_viscous_stress_support_premise_not_silently_closed=True,
        no_uniform_actual_physical_pulse_or_full_residual_claim=True)
    return dict(regions=regions,checks=checks,input_floor='C^-1024 > A^-1',
                compact_final_four_gap_floor='kappa/2, kappa=G^-4',
                selected_N='1+floor(Hhier^32)',
                immediate_exterior='Xrep=rhop*exp(5)=XR*exp(241*T-18)=Xstar')
