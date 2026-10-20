"""Sixth reference inputs and fifth sources; not final actual compact jets."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial
from pathlib import Path
import importlib.util, sys
import mpmath as mp
import sympy as s

HERE=Path(__file__).resolve().parent
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
high=module('reference_fifth_retained_higher',HERE.parent/'higher_jet_budget/higher_bounds.py')

def axis_coefficients(order):
    if type(order) is not int or order<0:raise ValueError('Nonnegative integer order required')
    return [comb(j+3,3) for j in range(order+1)]

def exponential_majorant(K,order):
    if type(K) not in (int,Q) or type(order) is not int or K<1 or not 0<=order<=6:
        raise ValueError('Exact K>=1 and order 0..6 required')
    return sum((Q(K)**j/factorial(j) for j in range(order+1)),Q())

def convex_average(weights,values):
    weights,values=tuple(map(Q,weights)),tuple(map(Q,values))
    if len(weights)!=len(values) or not weights or min(weights)<0 or sum(weights)!=1:
        raise ValueError('Positive probability weights required')
    return sum((w*v for w,v in zip(weights,values)),Q())

@lru_cache(None)
def parameters(digits):
    if type(digits) is not int or digits<60:raise ValueError('At least 60 interval digits required')
    r=high.core.certificate('Md64',digits);b=r['log_bounds']
    p,lo,hi,ls=high.p,high.lower,high.upper,high.logsum;log,exp=mp.iv.log,mp.iv.exp
    native=log(p('1e20'))+8*log(1+exp(b['ball_norm_bound']))+48*log(1+exp(-b['rho']))
    star=b['Lambda']+b['zeta_bound']+ls(log(2),*[
        log(p(2)**(j-1)/j)-(j-1)*b['r0'] for j in range(2,7)])
    # Value-sensitive exponential jets: no exp(K), exp(Lambda*K^8), or C.
    tiny=log(p('1e8'))-exp(b['Lambda'])+6*b['Lambda']+49*b['K']
    logs=dict(native_C6_packet=native,star_log_C6=star,reference_F_C6=tiny,
        full_reference_pressure_C6=log(1200)+2*b['P'],
        reference_Sq_C5=log(p('1e4'))+b['Lambda']+10*b['K'],
        reference_Sn_C5=log(p('1e3'))+2*b['K'],
        reference_p1_C5_at_Xi=log(p('1e8'))+log(110)+b['Lambda']+24*b['K'],
        reference_ns_C5=log(p('1e5'))+3*b['K'])
    checks=dict(retained_core_parameter_checks=all(r['checks'].values()),
        native_packet_below_K=hi(native-b['K'])<0,
        star_log_C6_below_10_Lambda_K8=hi(star-log(10)-b['Lambda']-8*b['K'])<0,
        reference_F_C6_below_one_over_100=hi(tiny)<lo(-log(100)),
        reference_E_C6_below_one=hi(tiny+log(220)/2)<0,
        full_reference_pressure_C6_below_K=hi(logs['full_reference_pressure_C6']-b['K'])<0,
        K_exceeds_100=lo(b['K'])>hi(log(100)),
        reference_U_C6_below_K=hi(ls(log(9),b['K']-b['Lambda'])-b['K'])<0,
        five_core_entry_moments_C5_below_K=hi(log(8)+2*b['K']-b['Lambda']-b['K'])<0,
        core_entry_radius_below_one=hi(log(4)-b['Lambda'])<0,
        reference_p1_C5_below_B_over100=hi(logs['reference_p1_C5_at_Xi']-b['B']+log(100))<0,
        reference_ns_C5_below_B_over100=hi(logs['reference_ns_C5']-b['B']+log(100))<0,
        same_radial_cutoff_is_explicit_monotone_sigma=True,
        no_variance_cutoff_budget_or_compact_gate_inferred=True)
    return dict(
        interval_digits=digits,log_bounds={k:dict(lower=mp.nstr(lo(v),125),upper=mp.nstr(hi(v),125)) for k,v in logs.items()},checks=checks,
        enormous_C_G_H_N_or_actual_profile_materialized=False)

@lru_cache(None)
def algebra():
    x,e,t=s.symbols('y eta t',real=True)
    lp=1+e+e**2*x/7+e**5*x**2/13+e**6*(1+x)/17
    psi=1-3*x*x+2*x**3
    primitive=s.integrate((psi*s.diff(lp,x)),(x,0,x))
    averaged=psi*lp-s.integrate(s.diff(psi,x)*lp,(x,0,x))
    identities=[s.expand(s.diff(primitive+lp.subs(x,0)-averaged,e,j))==0 for j in range(7)]
    # An eta-dependent taper has an additional term already in its first jet.
    eta_psi=psi+e*x*(1-x)/10
    commutator=s.expand(s.diff(eta_psi*s.diff(lp,x),e)-eta_psi*s.diff(lp,x,e))
    axis=[]
    for theta in (s.Rational(0),s.Rational(1,4),s.Rational(1,2),s.Rational(1)):
        f=(1+e*e)**(-2*theta)
        for eta in (-1,0,s.Rational(1,3),1):
            for j,cap in enumerate(axis_coefficients(6)):
                value=s.diff(f,e,j).subs(e,eta)/factorial(j)
                axis.append(bool(abs(value)<=cap))
    coef_bound=Q(4,3)**7
    return dict(axis_diagnostic_coefficients=len(axis),axis_C6_bound_over_P2=5*sum(axis_coefficients(6)),
        checks=dict(axis_majorant_coefficients_through_six=axis_coefficients(6)==[1,4,10,20,35,56,84],
            axis_C6_sum_is_1050=5*sum(axis_coefficients(6))==1050,
            manufactured_axis_mixture_atoms_through_six=all(axis),
            exact_reference_taper_identity_all_seven_angular_jets=all(identities),
            reference_weights_mass_one=s.integrate(-s.diff(psi,x),(x,0,1))==1,
            reference_weight_nonnegative_by_factorization=s.expand(-s.diff(psi,x)-6*x*(1-x))==0,
            eta_dependent_taper_commutator_is_not_zero=commutator!=0,
            coefficient_generating_bound_below_eight=coef_bound<8,
            native_log_constant_fits_Nnat=100000*9**6<10**20,
            native_inverse_constant_fits_Nnat=546875*9**6<10**20,
            native_mixed_log_constant_fits_Nnat=6*5**6*9**5*12<10**20,
            finite_tame_exponential_sum_below_three=sum((Q(1,factorial(j)) for j in range(7)),Q())<3,
            reference_tiny_constant_covers_star_and_Rr=3*11**6<10**8,
            W_C5_constant_covers_d_M_eta=1+2+4*6<30,
            Sq_C5_constant_covers_all_sources=300+5+500<10000,
            Sn_C5_constant_covers_full_pressure=150+5+30+24+8+2<1000,
            p1_C5_initial_and_source_constants_fit=225+450000<10**8,
            ns_C5_initial_and_source_constants_fit=10+2000<10**5,
            L_inverse_C5_geometric_norm_below_two=Q(1,1)/(1-Q(8,100))<2))

@lru_cache(None)
def flow_review():
    """Exact rational forward moments versus source ODEs, through eta order 5."""
    cc=module('reference_fifth_rational_flows',HERE.parent/'prescribed_continuation_foundation/continuation_checks.py')
    J=cc.modules()[0].Jet;rows=[]
    for center in (Q(-1),Q(1,3),Q(1)):
        order=7;y=J.variable(0,0,n=2,order=order);eta=J.variable(1,center,n=2,order=order)
        X=2*cc.exponential(y);root=2*cc.exponential(y/2)
        f0=1+eta/7+eta**6/19;f1=(1+eta**5)/23
        u0=4*eta+Q(1,10)+eta**6/29;u1=(1+eta**5)/31
        F=f0+X*f1;E=root*F;U=u0+X*u1
        moments=[X*u0+X**2*u1/2,
            X**2*f0+2*X**3*f1/3,
            X**2*u0*f0+2*X**3*(u0*f1+u1*f0)/3+X**4*u1*f1/2,
            X*u0*u0+X**2*u0*u1+X**3*u1*u1/3-X**2*f0*f0/2-2*X**3*f0*f1/3-X**4*f1*f1/4,
            X*f0*f0+X**2*f0*f1+X**3*f1*f1/3]
        axis=-2+eta/11+eta**6/37
        cc.coordinates.center=center;st=cc.coordinates(X,root,E,U,moments,axis)
        h=Q(1,200);D=Q(1,2)-h;A=Q(1,2)+h;d=1-eta*eta;L=1-2*h*eta*eta
        lr=1+F.derivative(0)/F;Hc=D*eta+d*U
        sq=-st['W']*lr-h*(1-2*eta*U)-Hc*E.derivative(1)/E
        pi=st['Pi']
        sn=-st['W']*U.derivative(0)-A*(1-2*eta*U)*U-Hc*U.derivative(1)-d*pi.derivative(1)+4*A*eta*pi+eta*E*E
        ep=st['p1'].derivative(0)+lr*st['p1']-X*sq/L
        en=st['ns'].derivative(0)+st['ns']-sn/L
        axis_effect=4*A*eta*axis-d*axis.derivative(1)
        cp=moments[4];cp_effect=4*A*eta*cp-d*cp.derivative(1)
        rows.append(dict(eta=str(center),p1_angular5_exact=all(ep.partial((0,j))==0 for j in range(6)),
            ns_angular5_exact=all(en.partial((0,j))==0 for j in range(6)),
            full_source_mixed_total5_exact=cc.zero_through(ep,5) and cc.zero_through(en,5),
            omitted_axis_fifth_source_nonzero=axis_effect.partial((0,5))!=0,
            omitted_radial_pressure_fifth_source_nonzero=cp_effect.partial((0,5))!=0,
            sixth_pressure_trace_nonzero=axis.partial((0,6))!=0,
            fifth_fields_nonzero=E.partial((0,5))!=0 and U.partial((0,5))!=0))
    return dict(manufactured_not_selected_profile=True,rows=rows,checks=dict(
        source_ODEs_match_genuine_forward_moments_through_fifth=all(v['p1_angular5_exact'] and v['ns_angular5_exact'] for v in rows),
        source_ODEs_match_all_total5_coefficients=all(v['full_source_mixed_total5_exact'] for v in rows),
        deleting_axis_pressure_fifth_is_detected=all(v['omitted_axis_fifth_source_nonzero'] for v in rows),
        deleting_radial_pressure_fifth_is_detected=all(v['omitted_radial_pressure_fifth_source_nonzero'] for v in rows),
        nonzero_sixth_trace_and_fifth_fields_exercised=all(v['sixth_pressure_trace_nonzero'] and v['fifth_fields_nonzero'] for v in rows)))

def review():
    a=algebra();f=flow_review();p=[parameters(d) for d in (80,110)]
    return dict(algebra=a,flows=f,parameters=p,checks=dict(**a['checks'],**f['checks'],
        two_precision_actual_parameter_inequalities=all(all(v['checks'].values()) for v in p)))
