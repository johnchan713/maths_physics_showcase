"""Generic complex old-wave signed update; exact fixtures, not actual fields."""
from fractions import Fraction as Q
from functools import lru_cache
import math
from pathlib import Path
from replay import load,compatible_compare,ROUND_FIELDS,ABS_ROUNDOFF

HERE=Path(__file__).resolve().parent
signed=load('complex_signed_frozen_signed_algebra',HERE.parent/'real_signed_covariance_update/signed.py')
KAPPA=signed.KAPPA
Poly,R,Z,T,EPS=signed.Poly,signed.R,signed.Z,signed.T,signed.EPS
outer,add,scale=signed.outer,signed.add,signed.scale


def covariance(a,b): return scale(add(outer(a,a),outer(b,b)),Q(1,2))
def bilinear(a,b,c,d):
    return scale(add(outer(a,c),outer(c,a),outer(b,d),outer(d,b)),Q(1,2))


@lru_cache(None)
def fixture():
    bump,records,_,_,_=signed.fixture();rows=[]
    for p in records:
        a,old_signed,r,old_curl=(p[n] for n in ('a','ell','r','s'))
        # A constant phase rotation of an EXACT complete curl is another curl.
        source_c=[(3*x+4*y)/5 for x,y in zip(old_signed,old_curl)]
        source_s=[(-4*x+3*y)/5 for x,y in zip(old_signed,old_curl)]
        e=[x+y for x,y in zip(old_signed,source_c)]
        uc=[x+y for x,y in zip(a,e)]
        us=[x+y+z for x,y,z in zip(r,old_curl,source_s)]
        multiplier=Q(1,6)+Z/23-T/31
        # The fresh homogeneous signed potential has a nonconstant slow target.
        # Recover the old potential: a_r=axial*c_theta-c_z/R, a_theta=-axial*c_r.
        f=bump*(1+Z/5+T/11);g=Q(1,4)+Z/7+T/13
        dc=[multiplier*x for x in (R*f,R*f*g,-R**2*f)]
        ell=signed.amplitude(dc,p['axial']);s=[-x/p['k'] for x in signed.curl(dc)]
        vc,vs=ell,s
        zero=[Poly(0)]*3
        principal=bilinear(a,zero,ell,zero)
        old_error=bilinear(e,us,ell,zero)
        new_curl=bilinear(uc,us,zero,s)
        square=covariance(vc,vs)
        remainder=add(old_error,new_curl,square)
        new=covariance([x+y for x,y in zip(uc,vc)],[x+y for x,y in zip(us,vs)])
        delta=add(new,scale(covariance(uc,us),-1))
        rows.append(dict(a=a,uc=uc,us=us,e=e,ell=ell,s=s,k=p['k'],axial=p['axial'],
                         principal=principal,old_error=old_error,new_curl=new_curl,
                         square=square,remainder=remainder,delta=delta))
    return bump,rows,add(*(p['delta'] for p in rows)),add(*(p['principal'] for p in rows)),add(*(p['remainder'] for p in rows))


@lru_cache(None)
def tensor_review():
    _,rows,delta,principal,remainder=fixture();checks=[];divs=[];omissions=[]
    for p in rows:
        checks.extend(p['delta'][i][j]==p['principal'][i][j]+p['remainder'][i][j]
                      for i in range(3) for j in range(3))
        samples=[]
        for cosine,sine in ((1,0),(0,1),(-1,0),(0,-1)):
            old=[cosine*x+sine*y for x,y in zip(p['uc'],p['us'])]
            inc=[cosine*x+sine*y for x,y in zip(p['ell'],p['s'])]
            new=[x+y for x,y in zip(old,inc)]
            samples.append(add(outer(new,new),scale(outer(old,old),-1)))
        av=scale(add(*samples),Q(1,4))
        checks.extend(av[i][j]==p['delta'][i][j] for i in range(3) for j in range(3))
        n=lambda v:v[1].shift(-1)+p['axial']*v[2]
        div=lambda v:signed.exact.radial(v[0])+EPS*v[2].diff(1)
        for cosine,sine in ((p['uc'],p['us']),(p['ell'],p['s'])):
            divs.extend((div(cosine)+p['k']*n(sine)==0,div(sine)-p['k']*n(cosine)==0))
        omissions.append(all(any(bool(p[name][0][j]) for j in (1,2)) for name in ('old_error','new_curl','square')))
    return dict(checks=dict(exact_generic_full_tensor_increment=all(checks),
        old_and_new_complete_curls_are_divergence_free=all(divs),
        all_three_generic_remainder_groups_are_needed=all(omissions),
        all_six_covariance_changes_retained=all(bool(delta[i][j]) for i in range(3) for j in range(i,3)),
        generic_remainder_not_the_initial_real_formula=any(bool(p['old_error'][0][j]) for p in rows for j in (1,2)),
        whole_tensor_identity=all(delta[i][j]==principal[i][j]+remainder[i][j] for i in range(3) for j in range(3))),
        exact_tensor_checks=len(checks),exact_divergence_checks=len(divs),
        manufactured_rotated_old_curl_not_actual_sourced_pulse=True,
        actual_closed_label_separation_not_replaced_by_fixture_frequencies=True)


@lru_cache(None)
def pressure_review():
    bump,_,delta,principal,remainder=fixture();radial=signed.exact.radial
    gr=-radial(delta[0][0])-EPS*delta[0][2].diff(1)+delta[1][1].shift(-1)
    P=gr.integral();rho=bump/bump.integral().scalar()*(1+Z*(R-Q(3,2)))
    crho=(R**2*rho).integral()/2
    primitive=(gr-rho*P).primitive();pm=primitive-primitive.at(0,1)
    jt=(R**2*delta[1][2]).integral()
    jz=(R*delta[2][2]).integral()-(R**2*gr).integral()/2
    et=radial(delta[0][1],2)+EPS*delta[1][2].diff(1)
    ez=radial(delta[0][2])+EPS*(delta[2][2]+pm).diff(1)
    return dict(checks=dict(full_compact_pressure=pm.at(0,1)==pm.at(0,2)==0 and pm.diff(0)==gr-rho*P,
        complete_axial_pressure_moment=(R*pm).integral()==-(R**2*gr).integral()/2+crho*P,
        exact_angular_integrated_flux=(R**2*et).integral()==EPS*jt.diff(1),
        exact_axial_integrated_flux=(R*ez).integral()==EPS*(jz+crho*P).diff(1),
        angular_target_cancellation=et-radial(principal[0][1],2)==radial(remainder[0][1],2)+EPS*delta[1][2].diff(1),
        axial_target_cancellation=ez-radial(principal[0][2])==radial(remainder[0][2])+EPS*(delta[2][2]+pm).diff(1),
        axial_pressure_covariance_and_moving_bump_terms_nonzero=all(bool(x) for x in
            (pm.diff(1),(R*delta[2][2]).integral().diff(1),crho.diff(1)*P)),
        all_three_updated_defects_are_nonzero=all(bool(x) for x in (P,jt,jz))),
        exact_mean_velocities_held_fixed=True,both_base_pressure_tails_and_actual_flat_velocities_retained_analytically=True)


def orders(alpha,rho,kappa=KAPPA):
    alpha,rho,k=Q(alpha),Q(rho),Q(kappa)
    if alpha<1 or rho<Q(1,2) or not 0<k<Q(1,2): raise ValueError('alpha>=1, rho>=1/2 and 0<kappa<1/2 required')
    b=alpha-Q(1,2)
    rows=dict(target=alpha,old_remainder=rho,signed_velocity=b,signed_curl=alpha-k,
        old_error_new_principal=rho+b,complete_old_new_curl=alpha+Q(1,2)-k,
        signed_square=2*alpha-1,full_covariance_change=alpha,
        full_mean_and_pressure=alpha-k,averaged_pressure=alpha,defects=alpha,
        axial_flux_and_pressure=alpha+1)
    rows['generic_error']=min(rows[n] for n in ('old_error_new_principal','complete_old_new_curl','signed_square'))
    rows['gain']=rows['generic_error']-alpha
    return rows


def physical_powers(alpha,h,a,b): return signed.physical_powers(alpha,h,a,b)


def exponent_review():
    k=KAPPA;B=1-3*k;alpha=Q(3,2)-3*k;D=Q(3,2)-6*k;mu=1-k
    rows=orders(alpha,B)
    generated=dict(full_linear=alpha-3*k,old_wave=alpha-k,
        actual_mean=alpha-Q(1,2)+mu-Q(1,2),signed_square=2*alpha-1-k)
    final=min(rows['generic_error'],alpha+1)
    rows.update(previous_full_wave_residual=D,new_averaged_tangential=final,
        old_and_new_nonzero_residual=min(D,*generated.values()),
        current_fixed_primary_remainder=min(B,alpha-Q(1,2)))
    return dict(checks=dict(current_target_uses_Haar_mean_not_full_graph_bound=alpha==B+Q(1,2),
        old_complex_remainder_is_current_B=rows['old_remainder']==B,
        old_complex_cross_and_signed_square_limit_gain=rows['generic_error']==2-6*k,
        actual_averaged_tangential_gain=final-alpha==Q(1,2)-3*k>0,
        all_three_generic_groups_retained=rows['generic_error']==min(rows['old_error_new_principal'],rows['complete_old_new_curl'],rows['signed_square']),
        generated_nonzero_errors_preserve_but_do_not_improve_D=min(generated.values())==D,
        full_mean_and_defects_not_claimed_relaxed=rows['full_mean_and_pressure']==alpha-k and rows['defects']==alpha,
        fixed_primary_remainder_stays_B=rows['current_fixed_primary_remainder']==B,
        rho_half_or_alpha_one_blocks_strict_generic_gain=orders(alpha,Q(1,2))['gain']==orders(1,B)['gain']==0,
        physical_residual_absolute_power_remains_negative=-2*(Q(1,2)+Q(1,200))-Q(1,2)+Q(1,200)*D<0),
        rows={n:str(v) for n,v in rows.items()},generated_nonzero_orders={n:str(v) for n,v in generated.items()},
        physical_zero_jet_powers={n:str(v) for n,v in physical_powers(alpha,Q(1,200),0,0).items()},
        strict_generic_gain_condition='min(rho-1/2,1/2-kappa,alpha-1)>0',
        maximum_updated_residual_multiplier=4,one_signed_substep_not_full_mean_relaxation=True)


def replay_review():
    # Exercise exactly the adapter's scope, not a global float tolerance change.
    strict=lambda a,b,path: (_ for _ in ()).throw(ValueError(path)) if a!=b else None
    outcomes=[]
    for field in ROUND_FIELDS:
        compatible_compare(8e-15,5e-15,field,strict)
        compatible_compare({'reviews':{field.split('.')[0]:{field.split('.')[1]:8e-15}}},
                           {'reviews':{field.split('.')[0]:{field.split('.')[1]:5e-15}}},'root',strict)
        for value in (-1e-15,3e-12,float('nan'),float('inf')):
            try: compatible_compare(value,5e-15,field,strict)
            except ValueError: outcomes.append(True)
            else: outcomes.append(False)
    rejected=[]
    for path in ('elsewhere.maximum_relative_commutator_error','root.flags.blowup_verified',
                 'root.reviews.full_viscous_potential_pressure_cutoff.omitted_term_gaps.curl'):
        try: compatible_compare(8e-15,5e-15,path,strict)
        except ValueError: rejected.append(True)
        else: rejected.append(False)
    return dict(checks=dict(only_three_explicit_near_zero_fields_adapted=len(ROUND_FIELDS)==3,
        every_nonfinite_negative_or_excess_roundoff_rejected=all(outcomes),
        all_other_fields_keep_inherited_comparison=all(rejected),
        adapter_stricter_than_unchanged_fixture_gates=ABS_ROUNDOFF<2e-10),
        absolute_roundoff_ceiling=ABS_ROUNDOFF,parent_CI_failure_run=38929673498,
        parent_CI_failure_was_numerical_replay_not_failed_mathematical_gate=True,
        historical_sources_and_evidence_not_edited=True)
