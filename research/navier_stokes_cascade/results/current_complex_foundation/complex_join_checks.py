"""Independent exact current multiharmonic curl and refreshed-mean controls.

Manufactured finite-regularity data, NOT evaluation of the actual giant fields.
"""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m


@lru_cache(None)
def source():return load('current_complex_first_sourced_checks',HERE.parent/'first_sourced_foundation/sourced_checks.py')


@lru_cache(None)
def means():return load('current_complex_mean_algebra',HERE.parent/'post_signed_mean_relaxation/relax.py')


@lru_cache(None)
def wave_fixture():
    s=source();e=s.exact();R,Z,T,eps=e.R,e.Z,e.T,e.EPS
    bump,primary,_,old,_,Wold,_,_,_=s.covariance_fixture()
    f=bump*(1+Z/5+T/7);factor=(F(2,3)+Z/11+T/13)*(F(1,7)+Z/37-T/41)
    c=[factor*x/3 for x in (R**2*f,R**2*f*(2+Z/7),-R**3*f)]
    L=[-3*(c[2].shift(-1)-c[1]),-3*c[0],3*c[0].shift(-1)]
    rest=[eps*c[1].diff(1),-eps*c[0].diff(1)+c[2].diff(0),-e.radial(c[1])]
    zero=[e.Poly()]*3
    new={1:([x+y for x,y in zip(old[1][0],L)],[x+y for x,y in zip(old[1][1],rest)]),2:old[2]}
    cov=lambda d:s.tensor_scale(s.tensor_add(*(s.tensor_add(s.outer(a,a),s.outer(b,b)) for a,b in d.values())),F(1,2))
    bil=lambda a,b,c,d:s.tensor_scale(s.tensor_add(s.outer(a,c),s.outer(c,a),s.outer(b,d),s.outer(d,b)),F(1,2))
    a=primary[1][0];err=[x-y for x,y in zip(old[1][0],a)]
    principal=bil(a,zero,L,zero);olderror=bil(err,old[1][1],L,zero)
    newcurl=bil(old[1][0],old[1][1],zero,rest);square=cov({1:(L,rest)})
    Wnew=cov(new);delta=s.tensor_add(Wnew,s.tensor_scale(Wold,-1))
    return old,new,Wold,Wnew,L,rest,principal,olderror,newcurl,square,delta


@lru_cache(None)
def wave_review():
    s=source();e=s.exact();old,new,Wold,Wnew,L,r,principal,olderror,newcurl,square,delta=wave_fixture()
    div=True;components=0
    for m,(a,b) in list(old.items())+[(1,(L,r))]:
        div &= e.radial(a[0])+e.EPS*a[2].diff(1)+3*m*(b[1].shift(-1)+b[2])==0
        div &= e.radial(b[0])+e.EPS*b[2].diff(1)-3*m*(a[1].shift(-1)+a[2])==0
        components+=2
    expected=s.tensor_add(principal,olderror,newcurl,square)
    # Independent Fourier convolution, including m=2, instead of a real-phase shortcut.
    count=0;fourier_ok=True
    for data,W in ((old,Wold),(new,Wnew)):
        c={}
        for m,(a,b) in data.items():
            c[m]=[(x/2,-y/2) for x,y in zip(a,b)];c[-m]=[s.conj(x) for x in c[m]]
        for i in range(3):
            for j in range(3):
                p=s.sumc([s.mul(c[m][i],c[-m][j]) for m in c],e.Poly())
                fourier_ok &= p[0]==W[i][j] and p[1]==0;count+=2
    return dict(checks=dict(generic_four_group_full_tensor_identity=all(delta[i][j]==expected[i][j] for i in range(3) for j in range(3)),
        independent_multiharmonic_Fourier_zero_mode=bool(fourier_ok),all_complete_curls_divergence_free=bool(div),
        all_six_tensor_changes_nonzero=all(bool(delta[i][j]) for i in range(3) for j in range(i,3)),
        each_generic_error_group_has_nonzero_radial_tangential_entries=all(any(bool(t[0][j]) for j in (1,2)) for t in (olderror,newcurl,square)),
        old_second_harmonic_not_deleted_from_current_covariance=all(new[2][q]==old[2][q] for q in (0,1)) and any(bool(x) for x in old[2][0]),
        deleting_old_complex_error_changes_increment=any(bool(x) for row in olderror for x in row),
        tensor_increment_not_all_just_target=any(delta[i][j]!=principal[i][j] for i in range(3) for j in range(3))),
        exact_Fourier_components=count,exact_divergence_components=components,
        fixture_potential_zero_extension='C5 polynomial potential; complete curl C4 only, not actual C-infinity',
        fixed_inverse_selected_target_identity_supplied_analytically_not_by_fixture=True)


@lru_cache(None)
def mean_fixture():
    m=means();P,R,Z,T=m.Poly,m.R,m.Z,m.T;b=(R-1)**6*(2-R)**6
    base,_,rho,_,_,_,_,_,current,_=m.fixture()
    # current contains a previous slow update AND previous fast updates.
    fresh=m.velocity(m.Harmonic(0,R*b*(1+Z/43+T/47)/19,R*b*(1+Z/53+T/59)/23),
        m.Harmonic(0,R*b*(1+Z/61)/29,R*b*(1+T/67)/31))
    d=m.velocity(m.Harmonic(R*b*(1+Z/71+T/73)/37),
        m.Harmonic((R**4*b).diff(0).shift(-2)*(1+Z/79+T/83)/41))
    fast=m.plus(current,fresh);final=m.plus(fast,d)
    _,_,Wold,Wnew,*_=wave_fixture()
    tensors=[[[P(p.terms) for p in row] for row in w] for w in (Wold,Wnew)]
    states=[m.state(base,current,tensors[0],rho)]+[m.state(base,v,tensors[1],rho) for v in (current,fast,final)]
    return base,rho,current,fresh,fast,d,final,tensors,states


@lru_cache(None)
def mean_review():
    m=means();R,Z,T,eps=m.R,m.Z,m.T,m.EPS
    base,rho,current,u,fast,d,final,tensors,ss=mean_fixture();b,V,G=base
    beta,v,g=current;ub,uv,ug=u
    fg=(-m.radial(2*m.avg_product(beta,ub)+m.avg_product(ub,ub))
        -eps*(m.avg_product(beta,ug)+m.avg_product(g,ub)+m.avg_product(ub,ug)).diff(1)
        +(2*m.avg_product(v,uv)+m.avg_product(uv,uv)).shift(-1))
    fjt=(R**2*(m.avg_product(g,uv)+m.avg_product(v,ug)+m.avg_product(ug,uv))).integral()
    fjz=(R*(2*m.avg_product(g,ug)+m.avg_product(ug,ug))).integral()-(R**2*fg).integral()/2
    B,v0,g0=[x.mean for x in fast];db,dv,dg=[x.mean for x in d]
    linear=(2*V*dv).shift(-1)
    rg=(eps*db.diff(2)-m.radial(2*b*db+2*B*db+db**2)
        -eps*(b*dg+G*db+B*dg+g0*db+db*dg).diff(1)
        +(2*v0*dv+dv**2).shift(-1)+eps*(m.lap(db)-db.shift(-2)))
    sjt=(R**2*(g0*dv+v0*dg+dg*dv)).integral()
    sjz=(R*(2*g0*dg+dg**2)).integral()-(R**2*rg).integral()/2
    return dict(checks=dict(fast_all_three_actual_velocities_zero_Haar=all(x.mean==0 for x in u),
        fast_and_slow_complete_graph_curl_divergence=all(x[0].dr(1)+eps*x[2].diff(1)==m.Harmonic() for x in (u,d)),
        fast_preserves_accumulated_nonzero_old_Haar_means=all(x.mean==y.mean and bool(y.mean) for x,y in zip(fast,current)),
        current_old_Haar_contains_previous_slow_increment=any(x.mean!=y.mean for x,y in zip(current,m.fixture()[3])),
        exact_fast_full_radial_source_change=ss[2]['gr']-ss[1]['gr']==fg,
        exact_fast_all_three_defect_changes=all(ss[2][n]-ss[1][n]==p for n,p in (('P',fg.integral()),('Jtheta',fjt),('Jz',fjz))),
        stale_prefast_targets_leave_three_nonzero_changes=all(ss[2][n]!=ss[1][n] for n in ('P','Jtheta','Jz')),
        old_new_fast_correlations_nonzero=bool(m.avg_product(v,uv)),
        exact_slow_remainder_with_all_old_means=ss[3]['gr']-ss[2]['gr']-linear==rg,
        old_slow_swirl_cross_cannot_be_deleted=bool((2*v0*dv).shift(-1).integral()),
        exact_slow_pressure_defect_difference=ss[3]['P']-ss[2]['P']-linear.integral()==rg.integral(),
        exact_slow_angular_defect_difference=ss[3]['Jtheta']-ss[2]['Jtheta']-(R**2*(G*dv+V*dg)).integral()==sjt,
        exact_slow_axial_defect_difference=ss[3]['Jz']-ss[2]['Jz']-(2*R*G*dg-R*V*dv).integral()==sjz,
        nonlinear_slow_defect_remainders_not_zero=all(bool(x) for x in (rg.integral(),sjt,sjz)),
        only_Haar_parts_survive_slow_crosses=all(m.avg_product(x,y)==x.mean*y.mean for x in fast for y in d)),
        wave_covariance_held_fixed_after_signed_update=True,
        manufactured_increments_NOT_actual_temporal_inverse_or_five_map_solution=True,
        actual_A1_and_fixed_I4_inverse_supplied_by_retained_analytical_maps=True)


@lru_cache(None)
def pressure_review():
    m=means();R,eps=m.R,m.EPS
    base,rho,current,u,fast,d,final,tensors,ss=mean_fixture();cr=(R**2*rho).integral()/2
    checks=[];moments=[]
    for vel,state in zip((current,current,fast,final),ss):
        moments.extend(((R**2*vel[1].mean).integral()==0,(R*vel[2].mean).integral()==0))
        checks.extend((state['pm'].at(0,1)==state['pm'].at(0,2)==0,
            state['pm'].diff(0)==state['gr']-rho*state['P'],
            (R*state['pm']).integral()==-(R**2*state['gr']).integral()/2+cr*state['P'],
            (R**2*state['Etheta']).integral()==eps*state['Jtheta'].diff(1),
            (R*state['Ez']).integral()==eps*(state['Jz']+cr*state['P']).diff(1)))
    dw=[[q-p for p,q in zip(a,b)] for a,b in zip(*tensors)]
    dg=-m.radial(dw[0][0])-eps*dw[0][2].diff(1)+dw[1][1].shift(-1)
    tail1=eps**2*(m.Z*m.T+m.Z**3);tail2=eps**4*(m.Z+m.T*m.Z**2);pB=R**2+m.Z**2+tail1+tail2
    return dict(checks=dict(all_four_recomputed_pressure_and_flux_states=all(checks),all_eight_velocity_moments_exact=all(moments),
        signed_covariance_change_with_old_means_recomputed=ss[1]['gr']-ss[0]['gr']==dg,
        all_six_current_covariance_entries_retained=all(bool(tensors[1][i][j]) for i in range(3) for j in range(i,3)),
        moving_bump_derivative_not_dropped=all(bool(cr.diff(1)*x['P']) for x in ss),
        every_mean_pressure_has_nonzero_axial_gradient=all(bool(x['pm'].diff(1)) for x in ss),
        pressure_changes_at_each_operation=all(a['pm']!=b['pm'] for a,b in zip(ss,ss[1:])),
        same_full_background_pressure_cancels_only_in_same_subtraction=(pB+ss[3]['pm']).diff(1)-pB.diff(1)==ss[3]['pm'].diff(1),
        leading_only_subtraction_leaves_both_tail_gradients=bool(tail1.diff(1)) and bool(tail2.diff(1))),
        exact_pressure_flux_checks=len(checks),exact_velocity_moment_checks=len(moments),
        toy_pressure_tails_NOT_numerical_evaluation_of_actual_Pi1_Pi2=True,
        complete_pressure_memory_and_flat_velocity_products_retained_analytically=True)


def reviews():return dict(current_multiharmonic_signed_tensor=wave_review(),
    fresh_fast_slow_mean_recomputation=mean_review(),four_state_complete_pressure_and_moments=pressure_review())
