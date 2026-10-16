"""Fresh mean relaxation of the fixed complex-wave state; exact fixtures."""
from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
KAPPA=Q(1,100000)


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


old=load('complex_mean_frozen_relaxation_algebra',HERE.parent/'post_signed_mean_relaxation/relax.py')
Poly,R,Z,T,EPS=old.Poly,old.R,old.Z,old.T,old.EPS
Harmonic,velocity,plus,avg_product=old.Harmonic,old.velocity,old.plus,old.avg_product
radial,lap=old.radial,old.lap


@lru_cache(None)
def complex_module():
    p=HERE.parent/'complex_signed_covariance_update';sys.path.insert(0,str(p))
    try:
        import covariance
        return covariance
    finally: sys.path.pop(0)


@lru_cache(None)
def fixture():
    base,_,rho,slow,_,_,_,previous_slow,current,_=old.fixture()
    cv=complex_module();_,records,_,_,_=cv.fixture()
    complete=cv.add(*(cv.covariance([x+y for x,y in zip(p['uc'],p['ell'])],
                                   [x+y for x,y in zip(p['us'],p['s'])]) for p in records))
    # The two frozen modules own distinct Poly classes; copy exact Q coefficients.
    w=[[Poly(p.terms) for p in row] for row in complete]
    bump=(R-1)**6*(2-R)**6
    fresh=velocity(Harmonic(0,R*bump*(1+Z/17+T/23)/7,R*bump*(1+Z/19+T/29)/11),
                   Harmonic(0,R*bump*(1+Z/31)/13,R*bump*(1+T/37)/17))
    delta=velocity(Harmonic(R*bump*(1+Z/17+T/29)/7),
                   Harmonic((R**4*bump).diff(0).shift(-2)*(1+Z/13+T/17)/11))
    fast=plus(current,fresh);final=plus(fast,delta)
    haar_before=plus(slow,previous_slow)
    states=[old.state(base,m,w,rho) for m in (current,fast,final)]
    return base,w,rho,haar_before,current,fresh,fast,delta,final,states


@lru_cache(None)
def fast_review():
    base,w,rho,haar,current,u,fast,delta,final,states=fixture()
    beta,v,gamma=current;ub,uv,ug=u
    expected=(-radial(2*avg_product(beta,ub)+avg_product(ub,ub))
              -EPS*(avg_product(beta,ug)+avg_product(gamma,ub)+avg_product(ub,ug)).diff(1)
              +(2*avg_product(v,uv)+avg_product(uv,uv)).shift(-1))
    change=states[1]['gr']-states[0]['gr']
    jt=(R**2*(avg_product(gamma,uv)+avg_product(v,ug)+avg_product(ug,uv))).integral()
    jz=(R*(2*avg_product(gamma,ug)+avg_product(ug,ug))).integral()-(R**2*change).integral()/2
    ci=Q(2,3);sv=-ci*uv.nt();sg=-ci*ug.nt()
    return dict(checks=dict(zero_Haar_actual_fast_increment=all(f.mean==0 for f in u),
        temporal_inverse_cancels_both_current_sources=all((-Q(1)/ci)*s.inverse_nt()==f for s,f in ((sv,uv),(sg,ug))),
        complete_graph_curl_divergence=ub.dr(1)+EPS*ug.diff(1)==Harmonic(),
        preserves_accumulated_old_slow_Haar_mean=all(f.mean==h.mean for f,h in zip(fast,haar)),
        accumulated_Haar_mean_contains_previous_slow_update=all(bool(h.mean) for h in haar) and any(h.mean!=s.mean for h,s in zip(haar,old.fixture()[3])),
        exact_fast_radial_source_change=change==expected,
        exact_fast_pressure_defect_change=states[1]['P']-states[0]['P']==expected.integral(),
        exact_fast_angular_defect_change=states[1]['Jtheta']-states[0]['Jtheta']==jt,
        exact_fast_axial_defect_change=states[1]['Jz']-states[0]['Jz']==jz,
        old_new_fast_correlation_is_nonzero=bool(avg_product(v,uv)),
        stale_prefast_targets_leave_three_nonzero_changes=all(states[1][n]!=states[0][n] for n in ('P','Jtheta','Jz'))),
        current_complex_wave_covariance_held_fixed=True,
        one_circle_manufactured_inverse_not_actual_global_two_torus_solve=True,
        actual_A1_flat_velocity_retained_by_analytical_lemmas=True)


@lru_cache(None)
def slow_review():
    base,w,rho,haar,current,u,fast,d,final,states=fixture()
    b,V,G=base;B,v,g=[f.mean for f in fast];db,dv,dg=[f.mean for f in d]
    linear=(2*V*dv).shift(-1)
    rg=states[2]['gr']-states[1]['gr']-linear
    without_old=(EPS*db.diff(2)-radial(2*b*db+db**2)
        -EPS*(b*dg+G*db+db*dg).diff(1)+(dv**2).shift(-1)+EPS*(lap(db)-db.shift(-2)))
    old_terms=-radial(2*B*db)-EPS*(B*dg+g*db).diff(1)+(2*v*dv).shift(-1)
    jt=(R**2*(g*dv+v*dg+dg*dv)).integral()
    jz=(R*(2*g*dg+dg**2)).integral()-(R**2*rg).integral()/2
    return dict(checks=dict(exact_slow_radial_remainder_with_accumulated_means=rg==without_old+old_terms,
        deleting_old_Haar_mean_changes_remainder=rg!=without_old,
        only_zero_Haar_parts_drop_against_slow_increment=all(avg_product(f,inc)==f.mean*inc.mean for f in fast for inc in d),
        exact_slow_graph_curl_divergence=d[0].dr(1)+EPS*d[2].diff(1)==Harmonic(),
        two_slow_velocity_moments_zero=(R**2*dv).integral()==(R*dg).integral()==0,
        exact_slow_pressure_defect_difference=states[2]['P']-states[1]['P']-linear.integral()==rg.integral(),
        exact_slow_angular_defect_difference=states[2]['Jtheta']-states[1]['Jtheta']-(R**2*(G*dv+V*dg)).integral()==jt,
        exact_slow_axial_defect_difference=states[2]['Jz']-states[1]['Jz']-(2*R*G*dg-R*V*dv).integral()==jz,
        old_swirl_cross_has_nonzero_radial_moment=bool((2*v*dv).shift(-1).integral()),
        old_axial_and_swirl_cross_fluxes_nonzero=bool((R**2*(g*dv+v*dg)).integral()),
        nonlinear_recomputed_defects_are_not_zero=all(bool(x) for x in (rg.integral(),jt,jz))),
        actual_fixed_five_map_cancels_current_postfast_linear_rows_analytically=True,
        manufactured_delta_not_numerical_actual_five_map_solution=True)


@lru_cache(None)
def pressure_review():
    base,w,rho,haar,current,u,fast,d,final,states=fixture()
    crho=(R**2*rho).integral()/2;checks=[];moments=[]
    for m,s in zip((current,fast,final),states):
        moments.extend(((R**2*m[1].mean).integral()==0,(R*m[2].mean).integral()==0))
        checks.extend((s['pm'].at(0,1)==s['pm'].at(0,2)==0,
            s['pm'].diff(0)==s['gr']-rho*s['P'],
            (R*s['pm']).integral()==-(R**2*s['gr']).integral()/2+crho*s['P'],
            (R**2*s['Etheta']).integral()==EPS*s['Jtheta'].diff(1),
            (R*s['Ez']).integral()==EPS*(s['Jz']+crho*s['P']).diff(1)))
    old_w=old.fixture()[1]
    return dict(checks=dict(pressure_and_integrated_fluxes_at_all_three_current_states=all(checks),
        all_six_velocity_moments_zero=all(moments),
        actual_fixture_covariance_updated_beyond_first_real_signed_state=all(bool(w[i][j]) for i in range(3) for j in range(i,3)) and any(w[i][j]!=old_w[i][j] for i in range(3) for j in range(3)),
        moving_bump_derivative_nonzero=bool(crho.diff(1)*states[2]['P']),
        axial_mean_pressure_gradient_nonzero=bool(EPS*(R*states[2]['pm']).integral().diff(1)),
        pressure_changes_after_both_mean_operations=states[0]['pm']!=states[1]['pm']!=states[2]['pm']),
        exact_pressure_flux_checks=len(checks),exact_velocity_moment_checks=len(moments),
        both_base_pressure_tails_and_actual_A0_A1_retained_analytically=True,
        flat_base_residual_not_assigned_shell_support=True)


def physical_powers(h,order,a,b): return old.physical_increment_powers(h,order,a,b)


def exponent_review():
    k=KAPPA;alpha=Q(3,2)-3*k;H=alpha-k;nu=2-6*k;mu=1-k;D=Q(3,2)-6*k
    fast=dict(slow_time=H+1,base_radial_flux=H+1-k,radial_viscosity=H+1-2*k,
        old_mean_transport=mu+H+1-k,changed_axial_pressure=H+1,new_mean_square=2*H+1-k)
    slow=dict(slow_time=alpha+1,base_radial_flux=alpha+1,radial_viscosity=alpha+1,
        old_mean_transport=mu+alpha+1-k,changed_axial_pressure=alpha+1,new_mean_square=2*alpha+1)
    scalar=dict(old_accumulated_swirl=alpha+1,new_slow_square=2*alpha,
        slow_time_and_base_flux=alpha+2,old_radial_mean_products=alpha+3,slow_viscosity=alpha+2)
    rows=dict(alpha=alpha,H=H,nu=nu,mu=mu,D=D,old_Haar_tangential=Q(1),old_Haar_radial=Q(2),
        fast_defect_change=mu+H,postfast_defects=alpha,slow_updated_defects=min(scalar.values()),
        final_full_tangential=min(nu,*fast.values(),*slow.values()),
        final_radial_residual=min(scalar.values()),added_nonzero_wave_source=min(H,alpha),
        final_nonzero_regular=min(D,H,alpha),next_sourced_source_order=D,
        current_fixed_primary_remainder=1-3*k)
    return dict(checks=dict(current_orders_joined_from_generic_signed_parent=alpha==Q(3,2)-3*k and nu==2-6*k,
        averaged_covariance_not_fast_viscosity_limits_full_mean=min(fast.values())>nu,
        fast_correlations_preserve_current_defect_class=mu+H>alpha,
        all_slow_tangential_changes_above_nu=min(slow.values())>nu,
        scalar_gain_limited_by_accumulated_Haar_mean=min(scalar.values())==alpha+1<2*alpha,
        stale_prefast_defect_change_is_below_new_target=mu+H<alpha+1,
        final_full_tangential_order_2_minus_6kappa=rows['final_full_tangential']==nu,
        final_defect_order_5_over_2_minus_3kappa=rows['slow_updated_defects']==Q(5,2)-3*k,
        unchanged_nonzero_residual_remains_limiting=rows['final_nonzero_regular']==D<nu<alpha+1,
        added_mean_wave_sources_keep_smooth_Gaussian_class=H>D and alpha>D,
        total_mean_orders_remain_mu_and_Haar_order_one=min(H,alpha)>1>mu,
        physical_residual_absolute_power_still_negative=-2*(Q(1,2)+Q(1,200))-Q(1,2)+Q(1,200)*D<0),
        rows={n:str(v) for n,v in rows.items()},fast_remainder_orders={n:str(v) for n,v in fast.items()},
        slow_tangential_orders={n:str(v) for n,v in slow.items()},slow_scalar_orders={n:str(v) for n,v in scalar.items()},
        physical_increment_zero_jet_powers={step:{n:str(v) for n,v in physical_powers(Q(1,200),t,0,0).items()}
            for step,t in (('fast',H),('slow',alpha))},
        fixed_derivative_and_support_bounds_retain_analytical_inverse_lemmas=True)
