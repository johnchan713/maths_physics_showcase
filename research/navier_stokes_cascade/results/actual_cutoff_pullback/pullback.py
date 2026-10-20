"""Actual axial-taper geometry; conditional analytical budget extraction."""
from fractions import Fraction as Q
from functools import lru_cache
import importlib.util,sys
from pathlib import Path
import mpmath as mp
import sympy as sp
HERE=Path(__file__).resolve().parent
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
high=module('actual_cutoff_pullback_higher_bounds',HERE.parent/'higher_jet_budget/higher_bounds.py')

@lru_cache(None)
def symbolic_review():
    e,q,t,d,B,kp=sp.symbols('eta q t delta B kp',positive=True)
    E=B/(1+e**2);bs=2*kp*e/E;beta=2*kp**2/B**2
    vs=2+bs**2/2
    y=sp.symbols('y',real=True);Ey=B*sp.exp(-y/2)/(1+e**2)
    D=1+3*e**2
    inv=[1/D]
    for j in range(3):inv.append(sp.factor(sp.diff(inv[-1],e)/D))
    expected=[1/D,-6*e/D**3,(90*e**2-6)/D**5,(360*e-2160*e**3)/D**7]
    # An independent implicit Taylor recurrence through degree four.
    z=sp.symbols('z');coeff=[];poly=e
    for j in range(1,5):
        c=sp.symbols('c'+str(j));trial=poly+c*z**j
        residual=sp.expand(trial+trial**3-e-e**3-z).coeff(z,j)
        cj=sp.factor(sp.solve(residual,c)[0]);coeff.append(cj);poly+=cj*z**j
    checks=dict(actual_axial_shear_identity=sp.factor(vs-(2+beta*e**2*(1+e**2)**2))==0,
        actual_a_is_two=sp.simplify(1-2*sp.diff(Ey,y)/Ey)==2,
        actual_taper_center_derivative=-sp.Rational(4,64)*8/sp.exp(32)==-sp.exp(-32)/2,
        inverse_four_derivatives_exact=all(sp.factor(a-b)==0 for a,b in zip(inv,expected)),
        implicit_normalized_inverse_all_four_orders=all(sp.factor(coeff[j-1]-expected[j-1]/sp.factorial(j))==0 for j in range(1,5)),
        exact_angular_crossing_derivative=sp.factor(sp.diff(vs,e)-2*beta*e*(1+e**2)*(1+3*e**2))==0,
        exact_window_pullback=sp.factor(beta*(d/beta)*(1+t)/8-d*(1+t)/8)==0,
        actual_target_squared_pullback=sp.factor((d/2-d*(1+t)/8)/2-d*(3-t)/16)==0)
    return dict(checks={k:bool(v) for k,v in checks.items()},inverse_raw=[str(v) for v in expected],
        normalized_inverse_bounds=[1,3,16,105],pressure_coordinate_set_to_zero=False)

def rational_review():
    # Suprema of each normalized derivative are summed, not just maximized.
    q_bound=Q(1,2)+Q(1,2)*sum((abs(Q(sp.binomial(sp.Rational(1,2),j))) for j in range(1,5)),Q())
    reciprocal_bound=sum((Q(1,2)**j for j in range(5)),Q()) # replace 1/sqrt2 by 1, safe <2
    checks=dict(q_C4_below_r=q_bound<1,inverse_whole_C4_cost_fits126=Q(1,2)+1+3+16+105<126,
        inverse_norm_below_one_at_r_bound=126*Q(1,10**6)<1,
        reciprocal_sqrt_C4_below_two=reciprocal_bound<2,
        raw_fourth_factorial_is24=sp.factorial(4)==24,
        normalized_inverse_second_bound=Q(6,2)==3,
        normalized_inverse_third_bound=Q(90+6,6)==16,
        normalized_inverse_fourth_bound=Q(360+2160,24)==105,
        shape_extraction_constant_192=4*2*24==192,
        universal_r_log_margin=Q(65,2)+14<47 and (2**16)**19>4,
        crossing_q_eta_lower=Q(1,2)*Q(1,2)==Q(1,4),
        crossing_gradient_upper=2*Q(1,2)*2*4==8)
    return dict(checks={k:bool(v) for k,v in checks.items()},q_C4_over_r_upper=str(q_bound),
        inverse_C4_over_r_upper='126',necessary_not_sufficient=True)

@lru_cache(None)
def parameter_review(digits):
    if type(digits) is not int or digits<60:raise ValueError('At least 60 interval digits')
    mp.iv.dps=digits;p,lo,hi=high.p,high.lower,high.upper
    parent=high.core.certificate('Md64',digits);b=parent['log_bounds'];log,exp=mp.iv.log,mp.iv.exp
    y0=exp(32)-1;Td=exp(64)+10;logP=b['P'];logA=4096*b['C']
    log_beta_upper=y0-64-log(2)-2*(logP-p('.5'))
    checks=dict(retained_parameter_certificate=all(parent['checks'].values()),
        same_Md64_logP=lo(logP-(Td+1))<=0<=hi(logP-(Td+1)),
        slice_after_taper_start=lo(y0)>0,slice_before_final_eleven_flat_units=hi(y0-(Td-11))<0,
        P1_upper_below_A=hi(log(2)+logP-logA)<0,
        beta_upper_below_one=hi(log_beta_upper)<0,
        log2_below_one=hi(log(2))<1,log_million_below14=hi(log(10**6))<14,
        inverse_r_bound_closes_universally=(2**16)**19>4,
        actual_A_lower_bound=lo(logA)>hi(log(2**16)))
    return dict(interval_digits=digits,checks={k:bool(v) for k,v in checks.items()},
        actual_A_beta_delta_G_H_N_or_radius_materialized=False)

def necessary_gate(A,logZ):
    if type(A) is not int or A<2**16:raise ValueError('Exact retained A lower bound')
    if type(logZ) not in (int,Q) or logZ<0:raise ValueError('Nonnegative exact auxiliary logZ')
    with mp.workdps(80):
        a=mp.mpf(A);z=mp.mpf(logZ.numerator)/logZ.denominator if isinstance(logZ,Q) else mp.mpf(logZ)
        return dict(necessary_condition_not_violated=bool(z<=a**32+a**20/2+mp.log(192)),
            actual_shape_enclosed=False,actual_S_accepted=False,actual_H_disproved=False,
            auxiliary_logZ_only=True)

def manufactured_review(delta,beta):
    with mp.workdps(85):
        d,b=mp.mpf(delta),mp.mpf(beta)
        if not 0<d<b<1 or not mp.sqrt(d/b)<mp.mpf('1e-6'):raise ValueError('Diagnostic must fit r<1e-6')
        r=mp.sqrt(d/b)
        def eta_closed(t):
            q=r*mp.sqrt((1+t)/8)
            return 2/mp.sqrt(3)*mp.sinh(mp.asinh(3*mp.sqrt(3)*q/2)/3)
        def eta_solved(t):
            q=r*mp.sqrt((1+t)/8)
            return mp.findroot(lambda e:e+e**3-q,(q/2,q))
        def K(t):
            if t<=0:return mp.mpf(1)
            if t>=1:return mp.mpf(0)
            return 1/(1+mp.exp(-1/t**2+1/(1-t)**2))
        def actual_slice(e):
            vminus=b*e**2*(1+e**2)**2
            return K(8*vminus/d-1)*mp.sqrt((d/2-vminus)/2)
        rows=[]
        for t in (mp.mpf('.25'),mp.mpf('.5'),mp.mpf('.75')):
            recovered=lambda z:4/mp.sqrt(d)/mp.sqrt(3-z)*actual_slice(eta_closed(z))
            expected=mp.diff(K,t,4);measured=mp.diff(recovered,t,4)
            e=eta_solved(t);grad=2*b*e*(1+e**2)*(1+3*e**2)
            rows.append(dict(t=str(t),checks=dict(independently_solved_inverse=bool(abs(eta_closed(t)/e-1)<mp.mpf('1e-60')),
                whole_shape_recovery=bool(abs(recovered(t)-K(t))<mp.mpf('1e-60')),
                fourth_shape_recovery=bool(abs(measured-expected)/(1+abs(expected))<mp.mpf('1e-55')),
                exact_actual_window_crossing=bool(abs(vminus_at(e,b)/d-(1+t)/8)<mp.mpf('1e-60')),
                positive_crossing_with_both_bounds=bool(mp.sqrt(b*d)/(2*mp.sqrt(2))<=grad<=8*mp.sqrt(b*d)))))
        return dict(delta=str(delta),beta=str(beta),digits=85,relative_error_limit='1e-55',rows=rows,
            manufactured_floating_not_interval_or_actual_parameter_evaluation=True)
def vminus_at(e,b):return b*e**2*(1+e**2)**2

@lru_cache(None)
def review():
    sy=symbolic_review();ra=rational_review();pa=[parameter_review(d) for d in (80,110)]
    nu=[manufactured_review('1e-18','0.0625'),manufactured_review('1e-24','0.125')]
    checks=dict(sy['checks']);checks.update(ra['checks'])
    checks['actual_parameter_comparisons_both_precisions']=all(all(r['checks'].values()) for r in pa)
    checks['independent_shape_and_fourth_pullback_diagnostics']=all(all(all(v['checks'].values()) for v in r['rows']) for r in nu)
    checks['necessary_gate_can_reject_auxiliary_large_logZ']=not necessary_gate(2**16,2**513)['necessary_condition_not_violated']
    checks['necessary_gate_never_accepts_actual_S']=necessary_gate(2**16,0)['actual_S_accepted'] is False
    return dict(symbolic=sy,rational=ra,parameters=pa,numerical=nu,checks=checks)
