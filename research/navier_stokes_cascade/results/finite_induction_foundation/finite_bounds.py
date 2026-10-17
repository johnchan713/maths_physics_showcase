"""Current-source half-line estimates, relative to retained analytical maps."""
from dataclasses import dataclass
from fractions import Fraction as F

KAPPA=F(1,100000)
B0=F(3,2)-2*KAPPA
RHO=MU=1-KAPPA
GAIN=F(1,2)-KAPPA


@dataclass(frozen=True)
class Affine:
    slope:F
    constant:F

    def at(self,b):return self.slope*F(b)+self.constant

    def margin(self,target,minimum=B0):
        slope=self.slope-target.slope
        if slope<0:raise ValueError('Negative slope fails on the whole half-line')
        return slope*F(minimum)+self.constant-target.constant


def line(s=0,c=0):return Affine(F(s),F(c))


def encode(v):
    if isinstance(v,F):return str(v)
    if isinstance(v,Affine):return dict(slope=str(v.slope),constant=str(v.constant))
    if isinstance(v,dict):return {k:encode(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [encode(x) for x in v]
    return v


def certificate(gain=GAIN):
    g=F(gain);k=KAPPA
    if g<=0:raise ValueError('Strictly positive gain required')
    wave_target=line(1,g);mean_target=line(1,g+F(1,2))
    wave=dict(sourced_full_linear=line(1,F(1,2)-k),sourced_old_wave=line(1,F(1,2)-k),
        sourced_square=line(2,-k),sourced_actual_mean=line(1,MU-F(1,2)),
        signed_full_linear=line(1,F(1,2)-k),signed_old_wave=line(1,F(1,2)-k),
        signed_square=line(2,-k),signed_actual_mean=line(1,MU-F(1,2)),
        fast_mean_wave=line(1,F(1,2)-k),slow_mean_wave=line(1,F(1,2)))
    averaged=dict(old_complex_cross=line(1,RHO),complete_old_new_curl=line(1,1-k),
        signed_square=line(2,0),current_moment_bump=line(1,F(3,2)),axial_covariance_and_pressure=line(1,F(3,2)))
    fast=dict(slow_time=line(1,F(3,2)-k),base_radial_flux=line(1,F(3,2)-2*k),
        radial_viscosity=line(1,F(3,2)-3*k),old_mean_transport=line(1,MU+F(3,2)-2*k),
        axial_pressure=line(1,F(3,2)-k),new_mean_square=line(2,2-3*k))
    slow=dict(slow_time=line(1,F(3,2)),base_radial_flux=line(1,F(3,2)),radial_viscosity=line(1,F(3,2)),
        old_mean_transport=line(1,MU+F(3,2)-k),axial_pressure=line(1,F(3,2)),new_mean_square=line(2,2))
    defects=dict(accumulated_Haar_cross=line(1,F(3,2)),new_slow_square=line(2,1),
        time_base_viscosity=line(1,F(5,2)),old_radial_mean=line(1,F(7,2)))
    groups=dict(nonzero=wave,averaged_signed=averaged,fast_full_mean=fast,slow_full_mean=slow,postslow_defects=defects)
    margins={group+':'+name:value.margin(wave_target if group=='nonzero' else mean_target)
        for group,rows in groups.items() for name,value in rows.items()}
    invariants=dict(new_wave_above_fixed_primary_error=B0-RHO,fast_mean_above_mu=B0+F(1,2)-k-MU,
        slow_mean_above_old_Haar_one=B0-F(1,2),fast_defect_change_above_C=MU-k,
        old_wave_error_above_half=RHO-F(1,2),strict_order_gain=g)
    if min(margins.values())<0:raise ValueError('Gain blocked by '+', '.join(n for n,v in margins.items() if v<0))
    if min(invariants.values())<=0:raise ValueError('A preserved-state invariant failed')
    return dict(minimum_B=B0,rho=RHO,mu=MU,gain=g,wave_target=wave_target,mean_target=mean_target,
        groups=groups,margins=margins,invariants=invariants,limiting_rows=sorted(n for n,v in margins.items() if v==0),
        proof_domain='every B>=3/2-2*kappa by nonnegative affine slopes and boundary margins')


def stage_order(j):
    if type(j) is not int or j<0:raise ValueError('Nonnegative finite integer stage required')
    return B0+j*GAIN


def harmonic_bound(j):stage_order(j);return 4*2**j


def next_harmonics(values):
    if not values or any(type(m) is not int or not m for m in values):raise ValueError('Nonempty nonzero integer harmonic set required')
    h=set(values)|{-1,1}
    return sorted(h|{m+n for m in h for n in h if m+n})


def physical_loss(h,m):
    h=F(h)
    if not 0<h<F(1,100) or type(m) is not int or m<0:raise ValueError('Construction h and fixed nonnegative jet required')
    return 2*(F(1,2)+h)+F(3,2)+2*m


def physical_power(j,h,m):return F(h)*stage_order(j)-physical_loss(h,m)


def first_stage(target,h=F(1,200),m=0):
    h=F(h);needed=(F(target)+physical_loss(h,m)-h*B0)/(h*GAIN)
    return max(0,-(-needed.numerator//needed.denominator))


def seed_review(parent):
    o=parent['reviews']['current_signed_and_mean_order_join']['orders']
    b=F(o['final_nonzero']);c=b+F(1,2)
    margins=dict(nonzero=b-B0,full_mean=F(o['final_full_tangential'])-c,
        defects=F(o['slow_defects'])-c,old_error=F(o['old_wave_error'])-RHO,
        mean=F(o['mu'])-MU,Haar_tangential=F(o['old_Haar_tangential'])-1,Haar_radial=F(o['old_Haar_radial'])-2)
    return dict(orders=encode(dict(B=b,C=c,rho=RHO,mu=MU,defects=F(o['slow_defects']))),margins=encode(margins),
        checks=dict(current_seed_exactly_supplies_new_B0=b==B0,full_mean_equals_B0_plus_half=F(o['final_full_tangential'])==c,
            defects_strictly_stronger_than_required=F(o['slow_defects'])>c,
            all_current_state_orders_satisfy_hypotheses=all(v>=0 for v in margins.values()),
            actual_flat_velocities_and_additive_Gaussian_ledger_retained=parent['flags']['all_accepted_physical_scales_closed_supports_and_flat_tail_ledger_retained'],
            current_fast_mean_wave_source_ties_D=parent['flags']['fast_mean_wave_source_only_ties_existing_D_not_strictly_higher']))


def half_line_review():
    c=certificate();k=KAPPA
    return dict(certificate=encode(c),gain_obstruction='s<=min(1/2-kappa,rho-1/2,mu-1/2)=1/2-kappa',
        checks=dict(all_complete_cycle_rows_on_whole_half_line=all(v>=0 for v in c['margins'].values()),
            all_preserved_state_margins_strict=all(v>0 for v in c['invariants'].values()),
            seven_nonzero_rows_tie_strongest_available_gain=len([n for n in c['limiting_rows'] if n.startswith('nonzero:')])==7,
            two_generic_covariance_crosses_limit_mean_gain=set(n for n in c['limiting_rows'] if n.startswith('averaged_signed:'))=={'averaged_signed:old_complex_cross','averaged_signed:complete_old_new_curl'},
            fast_viscosity_margin_is_half_minus_two_kappa=c['margins']['fast_full_mean:radial_viscosity']==F(1,2)-2*k>0,
            fresh_defects_beat_required_order=c['margins']['postslow_defects:accumulated_Haar_cross']==F(1,2)+k,
            stale_postfast_target_gap_stays_two_kappa=1-(MU-k)==2*k,
            sharper_gain_preserves_older_conservative_bound=GAIN>F(1,2)-3*k>0),
        exact_rational_budgets_check_estimates_not_continuum_suprema=True)


def physical_review():
    rows=[]
    for m,target in ((0,F(0)),(2,F(1)),(5,F(10))):
        j=first_stage(target,m=m)
        rows.append(dict(jet=m,target=str(target),first_finite_stage=j,power=str(physical_power(j,F(1,200),m)),
            preceding=str(physical_power(j-1,F(1,200),m)) if j else None))
    return dict(rows=rows,physical_bound='for every fixed finite j,m: |R(U_j)|_m<=C_jm*q^(h*B_j-K_m)',
        stage_order='B_j=3/2-2*kappa+j*(1/2-kappa)',K_m='2*A+3/2+2*m',harmonics='M_j<=4*2^j',
        checks=dict(strict_positive_finite_stage_gain=GAIN>0,
            every_fixed_target_has_some_finite_stage=all(F(r['power'])>=F(r['target']) and (r['preceding'] is None or F(r['preceding'])<F(r['target'])) for r in rows),
            full_low_stage_physical_power_still_negative=physical_power(0,F(1,200),0)<0,
            stage_independent_jet_loss=physical_loss(F(1,200),3)==F(851,100),
            finite_harmonics_double_under_exact_product_rule=harmonic_bound(0)==4 and harmonic_bound(1)==8,
            growing_stage_constants_can_prevent_fixed_q_contraction=-F(1,200)*GAIN<0),
        no_uniform_stage_constants_no_j_of_q_substitution_no_single_realization=True,
        actual_thousand_stage_fields_not_numerically_evaluated=True)
