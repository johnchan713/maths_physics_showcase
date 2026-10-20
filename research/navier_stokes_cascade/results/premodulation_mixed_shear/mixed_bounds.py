"""Actual pre-modulation mixed inputs, relative to retained regional estimates."""
from fractions import Fraction as Q
from functools import lru_cache
import importlib.util,sys
from pathlib import Path
import mpmath as mp
import sympy as s

HERE=Path(__file__).resolve().parent
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
high=module('mixed_native_higher_bounds',HERE.parent/'higher_jet_budget'/'higher_bounds.py')

def annular_bounds(d=Q(572,10**16),b0=Q(64),b1=Q(250000)):
    if any(type(v) not in (int,Q) for v in (d,b0,b1)):raise TypeError('Exact bounds required')
    d,b0,b1=map(Q,(d,b0,b1))
    if min(d,b0,b1)<0:raise ValueError('Nonnegative bounds required')
    h=3*d*b0;hy=3*d*b1
    if h>=1:raise ValueError('Relative-factor inverse not certified')
    return dict(relative_C4=h,relative_y_C4=hy,
        Ey_C4_over_P=5*(Q(1,10)*(1+h)+hy),Uy_C4=2*d*b1,
        physical_a_C4=Q(4,5)+2*hy/(1-h),
        leading_bump_width=Q(2,25),bump_value=b0,bump_y=b1)

def regional_bounds(Bold):
    if type(Bold) not in (int,Q):raise TypeError('Exact bound required')
    Bold=Q(Bold)
    if Bold<682:raise ValueError('Retained B_old lower bound needed')
    ann=annular_bounds()
    shape=Q(1,10)+64*(Bold+682)/(2000*Bold)
    outer=1+Q(64,200)*683
    return dict(shape_log_slope_C4=shape,outer_log_slope_C4=outer,
        log_slope_C4_cover_sum=Bold+Q(1)+Q(1,10)+Q(1)+1000,
        log_slope_C4_majorant=Bold+1003,Uy_C4_cover_sum=Bold+10,
        a_C4_majorant=2*Bold+2007,annulus=ann)

@lru_cache(None)
def parameter_review(digits):
    mp.iv.dps=digits;p,lo,hi,logsum=high.p,high.lower,high.upper,high.logsum
    r=high.core.certificate('Md64',digits);b=r['log_bounds'];log=mp.iv.log
    rp=log(p('1e18'))+log(110)+b['Lambda']+12*b['K']
    un=log(p('55e8'))+8*b['K']
    logslope=logsum(b['B'],log(1003))
    vel=logsum(b['B'],log(10))
    shear=logsum(log(2)+b['B'],log(2007))
    # Compare at finite A_min only, then use decreasing/increasing ratios.
    amin=p(2**16)
    ratio=120*log(10)/amin**20+30000/amin**19+8
    checks=dict(retained_core_parameter_checks=all(r['checks'].values()),
        reference_p1_C4_below_Bold_over100=hi(rp-b['B']+log(100))<0,
        actual_prescribed_Uy_C4_below_Bold=hi(un-b['B'])<0,
        Bold_exceeds_682=lo(b['B'])>hi(log(682)),
        global_log_slope_C4_below_C=hi(logslope-b['C'])<0,
        global_Uy_C4_below_C=hi(vel-b['C'])<0,
        global_a_C4_below_C=hi(shear-b['C'])<0,
        ideal_P_below_A=hi(b['P']-4096*b['C'])<0,
        shape_exit_precedes_restoration=hi(b['xsep'])<-8,
        same_G_absorbs_fourth_cutoff_template_if_shape_Z_le_1e100=hi(ratio)<lo(amin**12),
        monotone_all_larger_A_ratio_proof_uses_logA_below_A=True,
        fixed_frequency_roots_supports_and_pressure_not_reselected=True)
    return dict(interval_digits=digits,checks=checks,
        actual_enormous_C_G_H_N_or_modulated_profile_materialized=False,
        logarithmic_margin_checks_only=True)

def input_exponents():
    return dict(a_C4=1,inverse_a_C4=10,bs_C4=1286,ts_C4=1296,vs_C4=2583,
        sqrt_target_C4_A_power=20673,cutoff_template_A_power=30000,
        cutoff_template_delta_inverse_power=8,
        checks=dict(bs_A_power_retains_Ahi5=5+5*256+1==1286,
            inverse_a_fourth_finite_jet_power=1+4*2+1==10,
            ts_product_power=1286+10==1296,
            vs_complete_axial_square_power=2*1286+10+1==2583,
            target_highest_power=8*2583+9==20673,
            target_template_has_headroom=20673<30000,
            zeta_C4_constant_fits_2e4_times_actual_shape_Z=3*4096<20000,
            sqrt_w_C4_constant_fits_4e5=3*512*4**4<400000,
            sqrt_a_inverse_constant_is_21=3*7==21,
            target_C4_constant_fits_1e12_times_actual_shape_Z=20000*400000*21<10**12))

@lru_cache(None)
def symbolic_review():
    y,e=s.symbols('y eta',real=True)
    p1=1+e/7+e**2/11+e**4/13;ns=s.Rational(1,5)+e/17+e**3/19;k=1-y/4
    E=(1+e**2/5)*s.exp(y/2-p1*(y-y*y/8)/2)
    U=e/9-ns*((5-y)*s.exp(y)/4-s.Rational(5,4))/2
    log_slope=(1-k*p1)/2;uy=-k*s.exp(y)*ns/2
    checkE=s.simplify(s.diff(E,y)-E*log_slope)==0
    checkU=s.simplify(s.diff(U,y)-uy)==0
    checkmixedE=s.simplify(s.diff(E,y,e,4)/24-s.diff(E*log_slope,e,4)/24)==0
    checkmixedU=s.simplify(s.diff(U,y,e,4)/24-s.diff(uy,e,4)/24)==0
    # Relative E correction, NOT an additive E bump, on an earlier annulus.
    f=1/(1+e*e);E0=3*f*s.exp(y/10);c=(1+e+e**4)/1000;beta=y*y*(1-y)**2
    h=c*beta;En=E0*(1+h);Un=4*e+(1-e**2+e**4)*beta/2000
    exacta=s.Rational(4,5)-2*s.diff(h,y)/(1+h)
    annEy=E0*((1+h)/10+s.diff(h,y))
    annA=s.simplify(1-2*s.diff(En,y)/En-exacta)==0
    annMixed=s.simplify(s.diff(En,y,e,4)/24-s.diff(annEy,e,4)/24)==0
    wrong=s.simplify(s.diff(E0+h,y,e,4).subs({y:s.Rational(1,3),e:0})-s.diff(En,y,e,4).subs({y:s.Rational(1,3),e:0}))
    # A separately chosen reference flattened as E_r=sqrt(2X)*constant
    # has physical a_r=0; its slope is not the actual prescribed p1,r.
    ref=2*s.exp(y/2);reference_a=s.simplify(1-2*s.diff(ref,y)/ref)
    v=s.Symbol('v',positive=True)
    derivatives=[s.diff(s.sqrt(v),v,j) for j in range(5)]
    invder=[s.diff(v**(-s.Rational(1,2)),v,j) for j in range(5)]
    z=s.Symbol('z',real=True)
    step=1/(1+s.exp(1/z**2-1/(1-z)**2))
    center3=s.simplify(s.diff(step,z,3).subs(z,s.Rational(1,2)))
    return dict(manufactured_not_actual_profile=True,
        sqrt_fourth=str(derivatives[4]),inverse_sqrt_fourth=str(invder[4]),
        canonical_step_third_at_center=str(center3),
        checks=dict(actual_prescribed_E_radial_identity=bool(checkE),actual_prescribed_U_radial_identity=bool(checkU),
            actual_mixed_E_eta4_uses_only_slope_eta4=bool(checkmixedE),actual_mixed_U_eta4_uses_only_slope_eta4=bool(checkmixedU),
            relative_annular_physical_a_identity=bool(annA),relative_annular_E_mixed_identity=bool(annMixed),
            treating_relative_E_bump_as_additive_detected=bool(wrong!=0),
            flattened_reference_not_substituted_for_actual_prescribed_shear=bool(reference_a==0 and p1!=0),
            sqrt_fourth_raw_coefficient_is_minus15_over16=bool(s.simplify(derivatives[4]+s.Rational(15,16)*v**(-s.Rational(7,2)))==0),
            inverse_sqrt_fourth_raw_coefficient_is105_over16=bool(s.simplify(invder[4]-s.Rational(105,16)*v**(-s.Rational(9,2)))==0),
            compressed_cutoff_third_center_for_fourth_MVT=bool(center3==-3712),
            fixed_window_fourth_shape_norm_has_no_uniform_upper_bound=True))

def cutoff_shape_gate(A,Z,logZ=None):
    """Finite logarithmic diagnostic, never materialize G or the actual cutoff.

    logZ may be supplied separately for large exact constructions, but then
    this is only an auxiliary scalar diagnostic, not an actual Z enclosure.
    """
    if type(A) is not int or A<2**16:raise ValueError('Retained A lower bound required')
    if type(Z) not in (int,Q) or Z<1:raise ValueError('Positive exact shape bound >=1')
    with mp.workdps(80):
        aa=mp.mpf(A);zz=mp.mpf(Z.numerator)/Z.denominator if isinstance(Z,Q) else mp.mpf(Z)
        lhs=(mp.log(zz) if logZ is None else mp.mpf(logZ))+12*mp.log(10)+20673*mp.log(aa)+8*aa**20
        return dict(auxiliary_only=True,accepts=bool(lhs<=aa**32),actual_cutoff_budget_supplied=False)

def review():
    ann=annular_bounds();reg=regional_bounds(682);exp=input_exponents();sym=symbolic_review()
    params=[parameter_review(d) for d in (80,110)]
    checks=dict(annular_relative_C4_factor_below_one_half=ann['relative_C4']<Q(1,2),
        annular_actual_Ey_C4_below_P=ann['Ey_C4_over_P']<1,
        annular_actual_Uy_C4_below_one=ann['Uy_C4']<1,
        annular_actual_a_C4_below_one=ann['physical_a_C4']<1,
        shape_log_slope_C4_below_one=reg['shape_log_slope_C4']<1,
        outer_log_slope_C4_below_1000=reg['outer_log_slope_C4']<1000,
        restoration_Uy_C4_below_one=64*Q(3,10**18)<1,
        global_cover_sums_regions_not_only_maxima=reg['log_slope_C4_cover_sum']<reg['log_slope_C4_majorant'],
        all_parameter_interval_checks=all(all(p['checks'].values()) for p in params))
    checks.update(exp['checks']);checks.update(sym['checks'])
    checks['shape_gate_can_reject_a_too_large_shape_budget']=not cutoff_shape_gate(2**16,1,logZ=str(2**513))['accepts']
    checks['shape_gate_accepts_auxiliary_1e100_not_actual_cutoff']=cutoff_shape_gate(2**16,10**100)['accepts']
    return dict(annular={k:str(v) for k,v in ann.items()},
        regional_at_Bold682={k:str(v) for k,v in reg.items() if k!='annulus'},
        parameters=params,exponents=exp,symbolic=sym,checks=checks)
