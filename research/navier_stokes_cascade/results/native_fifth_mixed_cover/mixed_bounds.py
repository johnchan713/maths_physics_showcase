"""Actual PREMOD mixed fifth cover; no FINAL or actual-cutoff acceptance."""
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
high=module('fifth_mixed_higher',HERE.parent/'higher_jet_budget/higher_bounds.py')
pressure=module('fifth_mixed_full_pressure',HERE.parent/'native_mixed_fourth_boundary/boundary.py')

def add_power(*v):
    if not v or any(type(x) is not int or x<0 for x in v):raise ValueError('Nonnegative exact powers required')
    if len(v)>=10**9:raise ValueError('Too many terms for retained minimum C')
    return max(v)+1

@lru_cache(None)
def step_algebra():
    t,z=sp.symbols('t z');P=t;rows=[]
    for k in range(7):
        norm=sum(abs(c) for c in sp.Poly(P,t).all_coeffs())
        rows.append(dict(order=k,coefficient_l1=int(norm),bound=2**k*factorial(k),accepted=bool(norm<=2**k*factorial(k))))
        P=sp.expand(t*(1-t)*sp.diff(P,t))
    a=sp.symbols('a1:7');poly=sum(a[j-1]*z**j/factorial(j) for j in range(1,7))
    # Coefficient of z^6 in exp(poly) times 6! is the complete Bell polynomial.
    bell=sp.Poly(sp.expand(sum(poly**k/factorial(k) for k in range(7))),z).coeff_monomial(z**6)*factorial(6)
    bell=sp.Poly(sp.expand(bell),*a)
    return dict(logistic_rows=rows,Bell6_coefficient_sum=int(sum(bell.coeffs())),
        maximum_edge_derivative_constant=2*factorial(7),
        checks=dict(all_logistic_derivative_l1_caps=all(r['accepted'] for r in rows),
            complete_sixth_Bell_sum203=sum(bell.coeffs())==203,
            Bell6_coefficients_positive=all(c>0 for c in bell.coeffs()),
            Bell6_weighted_order_six=all(sum((i+1)*p for i,p in enumerate(mon))==6 for mon,_ in bell.terms()),
            Bell6_factor_count_at_most_six=all(sum(mon)<=6 for mon,_ in bell.terms())))

def step_comparisons(digits):
    if type(digits) is not int or digits<60:raise ValueError('At least 60 exact interval digits')
    mp.iv.dps=digits;p=mp.iv.mpf;exp=mp.iv.exp;log=mp.iv.log
    lo=lambda v:mp.make_mpf(v._mpi_[0]);hi=lambda v:mp.make_mpf(v._mpi_[1])
    middle=p(2)**6*factorial(6)*203*p('1e12')**6
    edge=p(406)*10080**6*4**18*exp(p(16)/9-16)
    return dict(interval_digits=digits,checks=dict(
        central_all_six_phase_derivatives_below_1e12=2*factorial(7)*4**8<10**12,
        central_Bell_bound_below_1e150=hi(middle)<lo(p('1e150')),
        edge_n6_geometric_ratio_below_half=hi(64*exp(-14))<lo(p('.5')),
        edge_maximum_t9_exp_minus_t_decreases=16>9,
        edge_Bell_bound_below_1e150=hi(edge)<lo(p('1e150')),
        fixed_sigma6_constant_fits_C17=hi(150*log(10)-17*log(p('1e9')))<0,
        relative_width_point08_J5_fits_C20=hi(156*log(10)-20*log(p('1e9')))<0),
        actual_huge_C_or_profile_materialized=False)

def exponent_ledger():
    rows=[]
    def row(name,raw,selected):
        rows.append(dict(name=name,raw_C_power=raw,selected_C_power=selected,accepted=raw<=selected));return selected
    lp=row('reference_positive_log_J4',add_power(2,75+4),128)
    u=row('reference_U_J4',add_power(1,75+4),128)
    e=row('reference_E_J4_value_sensitive',4*lp+1,1024)
    m=row('reference_M_J4_full_incoming',add_power(2,2+u),132)
    pi=row('reference_full_Pi_J4',add_power(1,2*e),2050)
    w=row('reference_W_J3',add_power(0,2+1+(m+1)),140)
    hc=row('reference_Hc_J3',add_power(1,1+u),130)
    sq=row('reference_Sq_J3',add_power(w+lp,1+u,hc+lp+1),272)
    sn=row('reference_Sn_J3',add_power(w+u+1,1+2*u,hc+u+1,1+pi+1,1+pi,1+2*e),2060)
    p1=row('reference_p1_J4_original_ODE',add_power(16,2+1+sq,lp+128),512)
    ns=row('reference_ns_J4_original_ODE',add_power(9,1+sn,192),4096)
    a=row('actual_prescribed_a_J4',add_power(100+p1,100),1024)
    b=row('actual_prescribed_b_J4_two_cutoffs',100+100+2+ns,8192)
    f=row('actual_positive_log_E_J5',add_power(1,a),2048)
    ee=row('actual_continuation_E_J5_value_sensitive',5*f+1,16384)
    uu=row('actual_continuation_U_J5',add_power(1,b),16384)
    cover=[16384,16384,128,128,256]
    total=row('SUM_five_region_field_J5',add_power(*cover),4096*5)
    return dict(rows=rows,region_C_powers=cover,full_native_field_J5_A_power=5,
        checks=dict(all_positive_exponent_rows_fit=all(r['accepted'] for r in rows),
            five_regional_norms_SUMMED_not_maximized=total==4096*5 and add_power(*cover)==16385,
            new_native_fifth_still_below_unchanged_Bn=5<256,
            radial_derivative_order_reduction_has_factor_at_most_five=5<2**16,
            field_to_shear_J4_constant_fits_one_A=11<2**16,
            native_shear_J4_A_power262=5+256+1==262,
            inverse_a_J4_A_power1054=1+4*(1+262)+1==1054,
            ts_J4_A_power1316=262+1054==1316,
            vs_J4_A_power1579=max(262,2*262+1054)+1==1579,
            full_pressure_J4_A1547_is_retained_not_rederived_or_reset=1547>256),
        shear_A_powers=dict(a=262,bs=262,inverse_a=1054,ts=1316,vs=1579,p1=1547,p2=1547))

@lru_cache(None)
def parameters(digits):
    if type(digits) is not int or digits<60:raise ValueError('At least 60 exact interval digits')
    r=high.core.certificate('Md64',digits);b=r['log_bounds'];log,exp=mp.iv.log,mp.iv.exp
    p,lo,hi=high.p,high.lower,high.upper;lc=b['C']
    checks=dict(retained_core_parameter_checks=all(r['checks'].values()),
        C_at_least_1e9=hi(log(p('1e9'))-lc)<0,
        K_old_below_C=hi(b['K']-lc)<0,
        native_star_packet_10LambdaK8_below_C2=hi(log(10)+b['Lambda']+8*b['K']-2*lc)<0,
        core_native_U_and_log_packet_below_C4=hi(log(p('1e6'))+b['Lambda']+8*b['K']-4*lc)<0,
        actual_pure_angular_log_CE_five_below_C=hi(b['B']-lc)<0,
        P_and_twice_P_below_C2=hi(log(2)+b['P']-2*lc)<0,
        P_below_C=hi(b['P']-lc)<0,
        reference_X_and_inverse_J3_below_C2=hi(log(330)+b['Lambda']-2*lc)<0,
        one_hundred_K_native_slope_constant_below_C4=hi(log(100)+b['K']-4*lc)<0,
        earlier_sigma5_cutoff_J3_C75=hi(log(p('1e102'))+60*lc-75*lc)<0,
        actual_sigma5_cutoff_J4_C100=hi(log(p('1e103'))+80*lc-100*lc)<0)
    g=p(2**16);checks.update(
        largest_new_native_power_below_G=hi(1579*log(g)-g**32)<0,
        ratio_A32_over_logA_increases=lo(32-1/log(g))>0,
        old_explicit_Z_gate_implies_new_joint_target_gate=20673-13685==6988 and lo(log(g))>0)
    return dict(interval_digits=digits,checks={k:bool(v) for k,v in checks.items()},
        actual_C_G_H_N_Z_or_Gamma_materialized=False)

def target_majorant(A,Z,delta_inverse):
    if any(type(v) not in (int,Q) for v in (A,Z,delta_inverse)):raise ValueError('Exact auxiliary inputs required')
    A,Z,d=map(Q,(A,Z,delta_inverse))
    if A<2**16 or Z<1 or d<1:raise ValueError('Outside retained envelope')
    V=A**1579;aa=A**262
    raw=20000*Z*d**4*V**4*400000*d**4*V**4*21*A**5*aa**4
    majorant=10**12*Z*d**8*A**13685
    return dict(coefficient=20000*400000*21,power=8*1579+5+4*262,
        raw_below_majorant=raw<majorant,Z_linear=True,
        actual_Z_enclosed=False,actual_Gamma_or_FINAL_field_accepted=False)

@lru_cache(None)
def annular_principal():
    y,e,z=sp.symbols('y eta sigma6');d=sp.Rational(2,25)
    c=(1+e+e**2+e**3+e**4+e**5)/1000
    # Relative width-.08 bump has NO extra unit-mass 1/d prefactor.
    beta=1+3*y+5*y**2+7*y**3+11*y**4+z*y**5/(factorial(5)*d**5)
    E0=(1+y/10+y**2/200)/(1+e**2)
    E=E0*(1+c*beta);U=4*e+c*beta
    ec=sp.diff(E,y,5).subs({y:0,e:0})/120
    uc=sp.diff(U,y,5).subs({y:0,e:0})/120
    expected=sp.Rational(1,1000)/(120*d**5)
    return dict(checks=dict(relative_E_sigma6_principal_exact=sp.diff(ec,z)==expected,
        additive_U_sigma6_principal_exact=sp.diff(uc,z)==expected,
        dropping_sigma6_changes_fifth_radial_jet=sp.diff(ec,z)!=0,
        incorrect_unit_mass_prefactor_rejected=sp.diff(ec,z)!=expected/d,
        width_point08_not_later_point1=d!=sp.Rational(1,10)),
        sigma6_coefficient=str(expected),sixth_angular_root_derivative_required=False)

@lru_cache(None)
def reference_ns_conversion():
    y=sp.Symbol('y');X=sp.exp(y);E=sp.Function('E')(y);p2=sp.Function('p2')(y)
    Sn=sp.Function('Sn')(y);L=sp.Symbol('L',nonzero=True);ns=p2*E/X
    flow=X*Sn/(L*E)-sp.diff(E,y)*p2/E
    residual=sp.simplify((sp.diff(ns,y)+ns-Sn/L).subs(sp.diff(p2,y),flow))
    return dict(checks=dict(original_p2_to_ns_full_source_ODE_exact=residual==0,
        dropping_X_in_reference_ns_conversion_changes_flow=sp.simplify(ns-p2*E)!=0,
        dropping_E_in_reference_ns_conversion_changes_flow=sp.simplify(ns-p2/X)!=0))

@lru_cache(None)
def prescription_fixtures():
    cc=module('fifth_mixed_prescription_jets',HERE.parent/'prescribed_continuation_foundation/continuation_checks.py')
    J=cc.modules()[0].Jet;n=6;y=J.variable(0,0,n=2,order=n);e=J.variable(1,0,n=2,order=n)
    sy,se=sp.symbols('y eta');rows=[]
    pj=2+e/7+e**2/11+y*(1+e)/13+y**2*(1+e**2)/17+y**3/19
    nj=1+e/23+y*(1+e**2)/29+y**2/31+y**3*e/37
    ps=2+se/7+se**2/11+sy*(1+se)/13+sy**2*(1+se**2)/17+sy**3/19
    ns=1+se/23+sy*(1+se**2)/29+sy**2/31+sy**3*se/37
    for name,kj,bj,tj,ks,bs,ts,stop in (
        ('activation',1-y/4,1,0,1-sy/4,1,0,False),
        ('axial_taper',Q(1,10),1-y/3,0,sp.Rational(1,10),1-sy/3,0,False),
        ('shear_blend',Q(1,10),1,y/4,sp.Rational(1,10),1,sy/4,True)):
        aj=(1-tj)*kj*pj+Q(4,5)*tj;ass=(1-ts)*ks*ps+sp.Rational(4,5)*ts
        X=2*cc.exponential(y);bb=cc.jet(0,n) if stop else -kj*bj*X*nj/2
        Ej=2*cc.exponential(e/10+e**2/19+y/2-cc.primitive(aj)/2)
        Uj=4*e+e**5/41+cc.primitive(bb)
        fs=se/10+se**2/19+sy/2-sp.integrate(ass,(sy,0,sy))/2
        Es=sp.Poly(sp.expand(2*sum(fs**k/factorial(k) for k in range(6))),sy,se)
        ex=sum(sy**k/sp.factorial(k) for k in range(6))
        bss=0 if stop else -ks*bs*ex*ns
        Us=sp.Poly(sp.expand(4*se+se**5/41+sp.integrate(bss,(sy,0,sy))),sy,se)
        indices=[(2,3),(3,2),(4,1),(5,0)]
        exact=lambda jet,poly:all(jet.partial(ix)==factorial(ix[0])*factorial(ix[1])*Q(poly.coeff_monomial(sy**ix[0]*se**ix[1])) for ix in indices)
        rows.append(dict(stage=name,E_missing_four_indices_exact=exact(Ej,Es),U_missing_four_indices_exact=exact(Uj,Us),
            E_pure_radial_fifth_nonzero=Ej.partial((5,0))!=0,
            U_pure_radial_fifth_nonzero=Uj.partial((5,0))!=0,
            shear_E_total4_exact=cc.zero_through(Ej.derivative(0)-(1-aj)*Ej/2,4),
            shear_U_total4_exact=cc.zero_through(Uj.derivative(0)-bb,4)))
    return dict(manufactured_not_actual_parameter_profile=True,rows=rows,checks=dict(
        eight_missing_coefficients_match_independent_polynomial_expansion=all(r['E_missing_four_indices_exact'] and r['U_missing_four_indices_exact'] for r in rows),
        all_prescribed_shear_total4_jets_exact=all(r['shear_E_total4_exact'] and r['shear_U_total4_exact'] for r in rows),
        fifth_radial_E_nonzero_in_every_stage=all(r['E_pure_radial_fifth_nonzero'] for r in rows),
        fifth_radial_U_nonzero_before_axial_stop=all(r['U_pure_radial_fifth_nonzero'] for r in rows[:2]),
        stopped_U_has_zero_positive_radial_jets=not rows[-1]['U_pure_radial_fifth_nonzero']))

@lru_cache(None)
def review():
    sa=step_algebra();sc=[step_comparisons(d) for d in (80,110)];le=exponent_ledger();pa=[parameters(d) for d in (80,110)]
    an=annular_principal();fx=prescription_fixtures();pr=pressure.pressure_identity();ns=reference_ns_conversion();ta=target_majorant(2**16,7,Q(3,2))
    checks={**sa['checks'],**le['checks'],**an['checks'],**fx['checks'],**pr['checks'],**ns['checks']}
    checks.update(sixth_step_outward_comparisons=all(all(r['checks'].values()) for r in sc),
        actual_parameter_comparisons=all(all(r['checks'].values()) for r in pa),
        joint_target_power_and_constant_exact=ta['power']==13685 and ta['coefficient']==168000000000,
        joint_target_auxiliary_positive_majorant=ta['raw_below_majorant'])
    return dict(step=sa,step_intervals=sc,exponents=le,parameters=pa,annular=an,fixtures=fx,
        complete_pressure_identity=pr,reference_ns=ns,target=ta,checks=checks)
