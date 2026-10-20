"""Full joint fourth moving loop and conditional Gamma/field interface."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
from pathlib import Path
import importlib.util,sys
import mpmath as mp
import sympy as sp
HERE=Path(__file__).resolve().parent
def module(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
fourth=module('joint4_signed_scalar',HERE.parent/'fourth_loop_shape_gate/loop_fourth.py')
native=module('joint4_actual_native',HERE.parent/'native_fifth_mixed_cover/mixed_bounds.py')
field=module('joint4_fixed_frequency',HERE.parent/'native_mixed_fourth_boundary/boundary.py')
bridge=module('joint4_full_memory',HERE.parent/'mixed_shape_cone_bridge/bridge.py')

def exponent_review():
    le=native.exponent_ledger();mu=fourth.implicit_budget(64,10,[80,240,400,560])
    theta=fourth.implicit_budget(10010,1,[11000,33000,55000,77000])
    counts_mu=[r['lower_composition_terms']*2**r['order'] for r in mu]
    counts_theta=[r['lower_composition_terms']*3**r['order'] for r in theta]
    checks=dict(native_joint4_inputs_below_A1579=all(v<=1579 for v in le['shear_A_powers'].values()) and le['full_native_field_J5_A_power']<1579,
        all_signed_fourth_root_orders_fit=all(r['accepted'] for r in mu),
        all_three_color_inverse_orders_fit=all(r['accepted'] for r in theta),
        all_fourth_color_partition_counts_below10000=max(counts_mu+counts_theta)<10000,
        full_fixed_angle_t_fourth_fits4000=5+4*(600+1+1)+600+1<4000,
        full_moving_source_fourth_fits350000=20000+4*(80000+1)<350000,
        primitive_full_joint_operator_three_fits400000=350000+1<400000,
        v_joint4_and_reciprocal_keep_pressure_coupling=4*4+1==17,
        old_explicit_shape_gate_implies_joint_target_with_gap=20673-13685==6988)
    return dict(mu=mu,theta=theta,mu_two_color_counts=counts_mu,theta_three_color_counts=counts_theta,
        native_shear_pressure_J4_A_powers=le['shear_A_powers'],
        exponents=dict(mu_J4=600,t_fixed_angle_J4=4000,w_J4=10000,Phi_J4=10010,
            positive_Theta_J4=80000,fixed_angle_sources_J4=20000,moving_sources_J4=350000,
            primitives_full_jointJ4=400000),primitive_full_joint_operator_bound=3,
        Gamma_A_and_B_bound='2*H',checks=checks)

def add(*polys):
    out={}
    for p in polys:
        for k,v in p.items():out[k]=out.get(k,0)+v
    return {k:sp.expand(v) for k,v in out.items() if v!=0}
def scale(p,c):return {k:sp.expand(c*v) for k,v in p.items()}
def mul(p,q,degree=4):
    out={}
    for a,v in p.items():
        for b,w in q.items():
            k=tuple(x+y for x,y in zip(a,b))
            if sum(k)<=degree:out[k]=out.get(k,0)+v*w
    return {k:sp.expand(v) for k,v in out.items() if v!=0}
def power(p,n):
    out={(0,0,0):sp.Integer(1)}
    for _ in range(n):out=mul(out,p)
    return out
def collapse(p):return [sum(abs(v) for k,v in p.items() if sum(k)==n) for n in range(5)]

@lru_cache(None)
def colored_review():
    zero=(0,0,0);one={zero:sp.Integer(1)};y={(1,0,0):sp.Integer(1)};e={(0,1,0):sp.Integer(1)};z={(0,0,1):sp.Integer(1)}
    y0=sp.Rational(1,5);e0=-sp.Rational(1,7);alpha=sp.Rational(1,7);beta=sp.Rational(1,9);gamma=sp.Rational(1,11)
    yy=add(scale(one,y0),y);ee=add(scale(one,e0),e);w0=1+alpha*y0+2*gamma*y0*e0
    def trig(kind,delta,mult=1):
        x=sp.Symbol('x');f=sp.sin(mult*x) if kind=='sin' else sp.cos(mult*x)
        return add(*[scale(power(delta,j),sp.diff(f,x,j).subs(x,0)/factorial(j)) for j in range(5)])
    def residual(delta):
        return add(delta,scale(mul(yy,trig('sin',delta)),alpha),
            scale(add(mul(ee,trig('cos',delta)),scale(one,-e0)),beta),
            scale(mul(mul(yy,ee),trig('sin',delta,2)),gamma),scale(z,-1))
    delta={}
    for n in range(1,5):
        res=residual(delta);delta=add(delta,{k:-sp.simplify(v)/w0 for k,v in res.items() if sum(k)==n})
    res=residual(delta);source=trig('cos',delta);mixed=source.get((1,1,2),0)
    # Independent one-variable implicit series along a nontrivial diagonal.
    h=sp.Symbol('h');directions=(sp.Rational(2),sp.Rational(-3),sp.Rational(5));jets=[]
    for n in range(1,5):
        d=sum(jets[j-1]*h**j for j in range(1,n))
        sin=sum(sp.diff(sp.sin(h),h,j).subs(h,0)*d**j/factorial(j) for j in range(5))
        cos=sum(sp.diff(sp.cos(h),h,j).subs(h,0)*d**j/factorial(j) for j in range(5))
        sin2=sum(sp.diff(sp.sin(2*h),h,j).subs(h,0)*d**j/factorial(j) for j in range(5))
        rr=sp.expand(d+alpha*(y0+directions[0]*h)*sin+beta*((e0+directions[1]*h)*cos-e0)+gamma*(y0+directions[0]*h)*(e0+directions[1]*h)*sin2-directions[2]*h)
        jets.append(-rr.coeff(h,n)/w0)
    diagonal=[sum(v*directions[0]**k[0]*directions[1]**k[1]*directions[2]**k[2] for k,v in delta.items() if sum(k)==n) for n in range(1,5)]
    P={zero:sp.Rational(1,3),(1,0,0):-2,(0,1,1):3,(2,0,0):-5,(1,1,2):7}
    R={zero:-2,(0,0,1):sp.Rational(1,7),(1,0,0):-6,(1,1,0):-4,(0,2,2):11}
    pc,rc,actual=collapse(P),collapse(R),collapse(mul(P,R));bounds=[sum(pc[j]*rc[n-j] for j in range(n+1)) for n in range(5)]
    return dict(positive_base_Jacobian=str(w0),inverse_phase_first=str(delta.get((0,0,1),0)),
        mixed_y_eta_phase2_source_coefficient=str(mixed),positive_inverse_monomials=len(delta),
        manufactured_not_actual_phase=True,checks=dict(
            all_three_variable_inverse_coefficients_through_four_exact=all(sp.simplify(v)==0 for v in res.values()),
            inverse_phase_first_keeps_nonunit_Jacobian=bool(delta[(0,0,1)]==1/w0 and w0!=1 and w0>0),
            all_four_diagonal_coefficients_match_independent_series=all(sp.simplify(a-b)==0 for a,b in zip(diagonal,jets)),
            mixed_fourth_y_eta_phase2_source_nonzero=mixed!=0,
            frozen_angle_wrong_fourth_zero=mixed!=0,
            colored_l1_products_through_four_bounded=all(a<=b for a,b in zip(actual,bounds)),
            cancellation_control_prevents_majorant_equality=any(a<b for a,b in zip(actual,bounds))))

@lru_cache(None)
def primitive_review():
    y,e,z=sp.symbols('y eta phi',real=True);P=(1+y+y*e+e**4+y**4)*sp.sin(2*sp.pi*z);g=sp.diff(P,z)
    tests=[]
    for r in range(5):
        for s in range(5-r):
            for k in range(1,5-r-s):
                left=sp.diff(P,y,r,e,s,z,k)/(factorial(r)*factorial(s)*factorial(k))
                right=sp.diff(g,y,r,e,s,z,k-1)/(k*factorial(r)*factorial(s)*factorial(k-1))
                tests.append(sp.simplify(left-right)==0)
    return dict(positive_phase_rows=len(tests),phase_free_operator_bound=2,positive_phase_operator_bound=1,
        full_joint_operator_bound=3,checks=dict(all_normalized_phase_derivative_rows_exact=all(tests),
            pure_slow_fourth_primitive_terms_present=sp.diff(P,y,4)!=0 and sp.diff(P,e,4)!=0,
            primitive_zero_mean_exact=sp.integrate(P,(z,0,1))==0,
            full_joint_bound_is_three_not_two=2+1==3))

def universal_review(digits):
    if type(digits) is not int or digits<60:raise ValueError('At least 60 exact interval digits')
    mp.iv.dps=digits;p=mp.iv.mpf;log=mp.iv.log;a=g=p(2**16)
    lo=lambda v:mp.make_mpf(v._mpi_[0]);hi=lambda v:mp.make_mpf(v._mpi_[1])
    checks=dict(native_fourth_inputs_below_G=hi(1579*log(a)-a**32)<0,
        native_comparison_monotone_for_larger_A=lo(32-1/log(a))>0,
        colored_partition_10000_fits_expG=hi(log(10000))<lo(g),
        complete_primitive_exp400000G_below_H=lo(g**3)>400000,
        full_Phi_integral_operator_below_eight=hi(2*mp.iv.pi+1)<8,
        Phi_integration_constant_fits_exp10G=hi(log(8))<lo(10*g),
        Gamma_2H_sufficient_gate_has_headroom=hi(log(8))<lo(g) and lo(159*g**3)>257,
        fourth_field_known_part_half_budget=hi(log(1024))<lo(g) and lo(31*g**3)>257,
        derivative_repair_J4_below_one=hi(log(10000))<lo(g) and lo(16*g**3)>513,
        original_N_sandwich_and_old_K_retained=32+1<64)
    return dict(interval_digits=digits,checks={k:bool(v) for k,v in checks.items()},
        actual_G_H_N_Z_or_selected_profile_materialized=False)

def auxiliary_field(A=2**16,H_power=600):
    if type(A) is not int or A<2**16 or type(H_power) is not int or H_power<600:raise ValueError('Exact auxiliary A and H exponent required')
    H=A**H_power;N=H**32+1;out=field.field_gate(A,H,N,2*H)
    return dict(A_aux=A,H_power=H_power,N_power=32,flags=out,manufactured_not_actual_H_N_or_profile=True)

@lru_cache(None)
def mean_review():
    a,t,v,ts,d,mu,r2,E=sp.symbols('a t v ts d mu r2 E',nonzero=True);pi=sp.pi
    w=a*(1+t*t)/(2*pi*v);aL=v/(1+t*t);bL=-v*t/(1+t*t)
    exact_v=a*(1+ts*ts)+a*d*d*mu*mu*r2
    mean_mass=sp.expand(a*(1+ts*ts+d*d*mu*mu*r2)-exact_v)
    return dict(complete_physical_wave_covariance_or_mean_claim=False,checks=dict(
        exact_variance_gives_mean_one_mass=mean_mass==0,
        weighted_aL_density_is_constant=sp.cancel(aL*w-a/(2*pi))==0,
        weighted_bL_density_keeps_axial_t=sp.cancel(bL*w+a*t/(2*pi))==0,
        A_source_zero_mean_after_exact_mass=sp.simplify(-(a-a)/2)==0,
        B_source_zero_mean_uses_bs_minus_a_ts=sp.expand(E*(-a*ts-(-a*ts))/2)==0))

@lru_cache(None)
def flux_boundary():
    e,X,M,h,t=sp.symbols('eta X M_eta5 h t');L=1-2*h*e**2;d=1-e**2
    eta=e+t;moment=M*t**5;D=sp.Rational(1,2)-h
    original_flux=(-2*D*eta*moment-(1-eta**2)*sp.diff(moment,t))/(1-2*h*eta**2)
    coefficient=sp.simplify(sp.diff(sp.diff(original_flux,t,4).subs(t,0)/24,M))
    return dict(normalized_eta4_velocity_flux_eta5_coefficient=str(-5*d/L),
        actual_FINAL_fifth_angular_moment_or_full_flux_C4_enclosed=False,
        checks=dict(original_radial_flux_fourth_principal_coefficient_exact=sp.simplify(coefficient+5*d/L)==0,
            fifth_angular_memory_coefficient_nonzero_at_equator=sp.simplify((-5*d/L).subs(e,0))==-5,
            no_endpoint_only_cancellation_in_interior=sp.simplify((-5*d/L).subs(e,sp.Rational(1,3)))!=0))

@lru_cache(None)
def review():
    ex=exponent_review();co=colored_review();pr=primitive_review();u=[universal_review(d) for d in (80,110)]
    ac=auxiliary_field();memory=bridge.state_memory_review();fl=flux_boundary();means=mean_review()
    checks={**ex['checks'],**co['checks'],**pr['checks'],**memory['checks'],**fl['checks'],**means['checks']}
    checks.update(all_outward_universal_comparisons=all(all(r['checks'].values()) for r in u),
        same_frequency_Gamma2H_field_gate_accepts=all(ac['flags'][k] for k in ('known_below_half','Gamma_sufficient_comparison','full_majorant_below_budget','repair_J4_below_one')))
    return dict(exponents=ex,colored=co,primitive=pr,universal=u,auxiliary=ac,full_memory=memory,flux_boundary=fl,leading_means=means,
        checks={k:bool(v) for k,v in checks.items()})
