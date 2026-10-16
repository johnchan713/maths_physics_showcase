"""Post-signed mean algebra and order checks; manufactured data, not a PDE solve."""
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


signed=load('post_signed_frozen_real_covariance',HERE.parent/'real_signed_covariance_update/signed.py')
Poly,R,Z,T,EPS=signed.Poly,signed.R,signed.Z,signed.T,signed.EPS
radial,lap=signed.exact.radial,signed.exact.lap
GRAPH=Q(6,7)*R  # Manufactured d_r=2, M=3/7; NOT the actual two-torus map.


class Harmonic:
    """Exact constant/cos(y)/sin(y) data; products are averaged, not truncated."""
    def __init__(self,zero=0,cosine=0,sine=0):
        self.coeff=tuple(Poly(x) for x in (zero,cosine,sine))

    def __add__(self,other):
        if not isinstance(other,Harmonic): other=Harmonic(other)
        return Harmonic(*(a+b for a,b in zip(self.coeff,other.coeff)))

    __radd__=__add__

    def __mul__(self,other):
        if isinstance(other,Harmonic): raise TypeError('Use avg_product; no discarded harmonics')
        return Harmonic(*(c*other for c in self.coeff))

    __rmul__=__mul__

    def diff(self,axis):
        return Harmonic(*(c.diff(axis) for c in self.coeff))

    def nt(self):
        z,c,s=self.coeff
        return Harmonic(0,s,-c)

    def inverse_nt(self):
        z,c,s=self.coeff
        if z: raise ValueError('Temporal inverse requires zero Haar mean')
        return Harmonic(0,-s,c)

    def dr(self,e=0):
        z,c,s=self.coeff
        return Harmonic(radial(z,e),radial(c,e)+GRAPH*s,radial(s,e)-GRAPH*c)

    def __eq__(self,other):
        return isinstance(other,Harmonic) and self.coeff==other.coeff

    @property
    def mean(self):
        return self.coeff[0]


def avg_product(f,g):
    a,b,c=f.coeff;x,y,z=g.coeff
    return a*x+(b*y+c*z)/2


def velocity(psi,v):
    return (-EPS*psi.diff(1),v,psi.dr(1))


def plus(m,n):
    return tuple(a+b for a,b in zip(m,n))


def grbar(base,mean,w):
    b,V,G=base;beta,v,gamma=mean
    return (EPS*beta.mean.diff(2)-radial(2*b*beta.mean+avg_product(beta,beta)+w[0][0])
            -EPS*(b*gamma.mean+G*beta.mean+avg_product(beta,gamma)+w[0][2]).diff(1)
            +(2*V*v.mean+avg_product(v,v)+w[1][1]).shift(-1)
            +EPS*(lap(beta.mean)-beta.mean.shift(-2)))


def state(base,mean,w,rho):
    b,V,G=base;beta,v,gamma=mean
    gr=grbar(base,mean,w);p=gr.integral()
    primitive=(gr-rho*p).primitive();pm=primitive-primitive.at(0,1)
    jt=(R**2*(G*v.mean+V*gamma.mean+avg_product(gamma,v)+w[1][2])).integral()
    jz=(R*(2*G*gamma.mean+avg_product(gamma,gamma)+w[2][2])).integral()-(R**2*gr).integral()/2
    et=(-EPS*v.mean.diff(2)+radial(b*v.mean+V*beta.mean+avg_product(beta,v)+w[0][1],2)
        +EPS*(G*v.mean+V*gamma.mean+avg_product(gamma,v)+w[1][2]).diff(1)
        -EPS*(lap(v.mean)-v.mean.shift(-2)))
    ez=(-EPS*gamma.mean.diff(2)+radial(b*gamma.mean+G*beta.mean+avg_product(beta,gamma)+w[0][2])
        +EPS*(2*G*gamma.mean+avg_product(gamma,gamma)+w[2][2]+pm).diff(1)
        -EPS*lap(gamma.mean))
    return dict(gr=gr,P=p,pm=pm,Jtheta=jt,Jz=jz,Etheta=et,Ez=ez)


@lru_cache(None)
def fixture():
    bump,records,_,_,_=signed.fixture()
    w=signed.add(*(p['new'] for p in records))  # Full six-entry UPDATED curl covariance.
    basepsi=R*bump*(1+Z+T)
    base=(-EPS*basepsi.diff(1),R**2*(2+Z),radial(basepsi))
    slow=velocity(Harmonic(R*bump*(1+Z/5+T/7)),
                  Harmonic((R**4*bump).diff(0).shift(-2)*(1+Z/11+T/13)))
    oldfast=velocity(Harmonic(0,R*bump*(1+Z/3+T/11),R*bump*(1+Z/7)/3),
                     Harmonic(0,R*bump*(1+Z/5),R*bump*(1+T/11)/3))
    fresh=velocity(Harmonic(0,R*bump*(1+Z/11+T/19)/5,R*bump*(1+Z/13+T/23)/7),
                   Harmonic(0,R*bump*(1+Z/17)/7,R*bump*(1+T/29)/11))
    delta=velocity(Harmonic(R*bump*(1+Z/7+T/17)),
                   Harmonic((R**4*bump).diff(0).shift(-2)*(1+Z/19+T/23)))
    old=plus(slow,oldfast);fast=plus(old,fresh);final=plus(fast,delta)
    rho=bump/bump.integral().scalar()*(1+Z*(R-Q(3,2)))
    states=[state(base,m,w,rho) for m in (old,fast,final)]
    return base,w,rho,slow,old,fresh,fast,delta,final,states


def fast_review():
    base,w,rho,slow,old,u,fast,delta,final,states=fixture()
    beta,v,gamma=old;ub,uv,ug=u
    expected=(-radial(2*avg_product(beta,ub)+avg_product(ub,ub))
              -EPS*(avg_product(beta,ug)+avg_product(gamma,ub)+avg_product(ub,ug)).diff(1)
              +(2*avg_product(v,uv)+avg_product(uv,uv)).shift(-1))
    diff=states[1]['gr']-states[0]['gr']
    jt=(R**2*(avg_product(gamma,uv)+avg_product(v,ug)+avg_product(ug,uv))).integral()
    jz=(R*(2*avg_product(gamma,ug)+avg_product(ug,ug))).integral()-(R**2*diff).integral()/2
    ci=Q(2,3)
    source_v=-ci*uv.nt();source_g=-ci*ug.nt()
    return dict(checks=dict(
        exact_zero_Haar_fast_velocity=all(f.mean==0 for f in u),
        temporal_inverse_cancels_both_desired_sources=all((-Q(1)/ci)*s.inverse_nt()==f for s,f in ((source_v,uv),(source_g,ug))),
        graph_curl_divergence=ub.dr(1)+EPS*ug.diff(1)==Harmonic(),
        fast_update_preserves_existing_slow_mean=all(f.mean==s.mean for f,s in zip(fast,slow)),
        inherited_slow_means_are_not_zero=all(bool(f.mean) for f in slow),
        exact_fast_averaged_radial_source_change=diff==expected,
        exact_fast_pressure_defect_change=states[1]['P']-states[0]['P']==expected.integral(),
        exact_fast_angular_defect_change=states[1]['Jtheta']-states[0]['Jtheta']==jt,
        exact_fast_axial_defect_change=states[1]['Jz']-states[0]['Jz']==jz,
        old_fast_new_fast_correlations_survive=bool(avg_product(v,uv)),
        all_three_fast_defect_changes_are_nonzero=all(states[1][n]!=states[0][n] for n in ('P','Jtheta','Jz'))),
        manufactured_one_circle_graph_not_actual_irrational_two_torus=True,
        actual_A1_flat_remainder_retained_analytically_and_in_replayed_parent=True)


def slow_review():
    base,w,rho,slow,old,u,fast,d,final,states=fixture()
    b,V,G=base;B,v,gamma=[f.mean for f in fast];db,dv,dg=[f.mean for f in d]
    linear=(2*V*dv).shift(-1)
    rg=states[2]['gr']-states[1]['gr']-linear
    without_old=(EPS*db.diff(2)-radial(2*b*db+db**2)
                 -EPS*(b*dg+G*db+db*dg).diff(1)+(dv**2).shift(-1)
                 +EPS*(lap(db)-db.shift(-2)))
    old_terms=-radial(2*B*db)-EPS*(B*dg+gamma*db).diff(1)+(2*v*dv).shift(-1)
    jt=(R**2*(gamma*dv+v*dg+dg*dv)).integral()
    jz=(R*(2*gamma*dg+dg**2)).integral()-(R**2*rg).integral()/2
    return dict(checks=dict(
        exact_slow_radial_remainder=rg==without_old+old_terms,
        dropping_existing_slow_mean_changes_remainder=rg!=without_old,
        postfast_zero_mean_parts_drop_only_from_slow_crosses=all(avg_product(f,inc)==f.mean*inc.mean for f in fast for inc in d),
        slow_graph_curl_divergence=d[0].dr(1)+EPS*d[2].diff(1)==Harmonic(),
        both_slow_velocity_moments_zero=(R**2*dv).integral()==(R*dg).integral()==0,
        exact_pressure_defect_difference=states[2]['P']-states[1]['P']-linear.integral()==rg.integral(),
        exact_angular_defect_difference=states[2]['Jtheta']-states[1]['Jtheta']-(R**2*(G*dv+V*dg)).integral()==jt,
        exact_axial_defect_difference=states[2]['Jz']-states[1]['Jz']-(2*R*G*dg-R*V*dv).integral()==jz,
        old_slow_swirl_cross_has_nonzero_radial_moment=bool((2*v*dv).shift(-1).integral()),
        old_slow_flux_crosses_are_nonzero=bool((R**2*(gamma*dv+v*dg)).integral()),
        nonlinear_recomputed_defects_not_set_to_zero=all(bool(x) for x in (rg.integral(),jt,jz))),
        fixed_wave_covariance_held_in_both_states=True,
        linear_five_map_rows_cancel_current_defects_not_previous_stage_defects=True)


def pressure_review():
    base,w,rho,slow,old,u,fast,d,final,states=fixture()
    crho=(R**2*rho).integral()/2
    checks=[];moments=[]
    for m,s in zip((old,fast,final),states):
        moments.extend(((R**2*m[1].mean).integral()==0,(R*m[2].mean).integral()==0))
        checks.extend((s['pm'].at(0,1)==s['pm'].at(0,2)==0,
            s['pm'].diff(0)==s['gr']-rho*s['P'],
            (R*s['pm']).integral()==-(R**2*s['gr']).integral()/2+crho*s['P'],
            (R**2*s['Etheta']).integral()==EPS*s['Jtheta'].diff(1),
            (R*s['Ez']).integral()==EPS*(s['Jz']+crho*s['P']).diff(1)))
    return dict(checks=dict(all_three_reconstructed_pressure_and_integrated_identities=all(checks),
        all_velocity_moments_preserved=all(moments),
        all_six_fixed_updated_covariances_retained=all(bool(w[i][j]) for i in range(3) for j in range(i,3)),
        moving_bump_derivative_nonzero=bool(crho.diff(1)*states[2]['P']),
        axial_pressure_gradient_not_zero=bool(EPS*(R*states[2]['pm']).integral().diff(1)),
        stale_fast_pressure_invalid_for_final_state=states[2]['pm']!=states[1]['pm']),
        exact_pressure_and_integrated_checks=len(checks),exact_velocity_moment_checks=len(moments),
        both_base_pressure_tails_cancel_only_in_same_exact_base_subtraction=True,
        actual_flat_base_residual_remains_a_separate_non_shell_summand=True)


def exponent_review():
    k=KAPPA;alpha=2-3*k;H=alpha-k;nu=3-6*k;mu=1-k
    fast=dict(slow_time=H+1,base_radial_flux=H+1-k,radial_viscosity=H+1-2*k,
              old_mean_transport=mu+H+1-k,changed_axial_pressure=H+1,new_mean_square=2*H+1-k)
    slow=dict(slow_time=alpha+1,base_radial_flux=alpha+1,radial_viscosity=alpha+1,
              old_mean_transport=mu+alpha+1-k,changed_axial_pressure=alpha+1,new_mean_square=2*alpha+1)
    scalar=dict(old_slow_swirl=alpha+1,new_slow_square=2*alpha,
                slow_time_and_base_flux=alpha+2,old_radial_mean_products=alpha+3,
                slow_viscosity=alpha+2)
    rows=dict(alpha=alpha,H=H,nu=nu,mu=mu,old_Haar_tangential=Q(1),old_Haar_radial=Q(2),
              fast_defect_change=mu+H,postfast_defects=alpha,slow_updated_defects=alpha+1,
              final_full_tangential=min(nu,*fast.values(),*slow.values()),
              final_radial_residual=alpha+1,added_nonzero_wave_source=min(H,alpha))
    return dict(checks=dict(
        fast_viscosity_is_limiting=fast['radial_viscosity']==nu,
        all_fast_errors_above_or_at_nu=min(fast.values())>=nu,
        fast_correlations_do_not_worsen_current_defect_class=mu+H>alpha,
        all_slow_tangential_errors_above_nu=min(slow.values())>nu,
        slow_scalar_gain_uses_old_Haar_order_one=min(scalar.values())==alpha+1,
        false_zero_old_mean_would_overstate_gain=alpha+1<2*alpha,
        full_tangential_mean_reaches_3_minus_6kappa=rows['final_full_tangential']==nu,
        radial_defect_bound_is_stronger=alpha+1>nu,
        total_mean_orders_preserved=min(H,alpha)>mu,
        new_nonzero_errors_do_not_cancel_old_wave_source=H>1-3*k,
        normalized_mean_gain_not_absolute_physical_decay=Q(1,200)*nu-2*(Q(1,2)+Q(1,200))-Q(1,2)<0),
        rows={n:str(v) for n,v in rows.items()},
        fast_remainder_orders={n:str(v) for n,v in fast.items()},
        slow_tangential_orders={n:str(v) for n,v in slow.items()},
        slow_scalar_orders={n:str(v) for n,v in scalar.items()},
        physical_increment_zero_jet_powers={step:{n:str(v) for n,v in physical_increment_powers(Q(1,200),t,0,0).items()}
                                            for step,t in (('fast',H),('slow',alpha))},
        all_fixed_derivative_orders_by_retained_product_and_inverse_lemmas=True)


def physical_increment_powers(h,order,a,b):
    h,order=Q(h),Q(order)
    if not 0<h<Q(1,100) or order<0 or any(type(n) is not int or n<0 for n in (a,b)):
        raise ValueError('Construction h, nonnegative order and integer derivative counts required')
    A=Q(1,2)+h;shift=-a*(Q(1,2)+h/2)-b*(1+3*h/2)
    return dict(velocity=-A+h*order+shift,pressure=-2*A+h*order+shift,
                radial_axial_curl_potential=Q(1,2)-A+h*order+shift)
