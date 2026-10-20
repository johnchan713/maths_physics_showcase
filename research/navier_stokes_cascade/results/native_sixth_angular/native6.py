"""Native-only sixth angular fields, full memory, and fifth R5 inputs."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb,factorial
from pathlib import Path
import importlib.util,sys
import mpmath as mp
import sympy as s
HERE=Path(__file__).resolve().parent
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m
ref=module('native6_retained_reference',HERE.parent/'reference_fifth_transfer/reference_bounds.py')
high=ref.high

@lru_cache(None)
def parameters(digits):
    if type(digits) is not int or digits<60:raise ValueError('At least 60 exact interval digits')
    r=high.core.certificate('Md64',digits);b=r['log_bounds']
    p,lo,hi,ls=high.p,high.lower,high.upper,high.logsum;log,exp=mp.iv.log,mp.iv.exp
    N7=log(p('1e30'))+10*log(1+exp(b['ball_norm_bound']))+64*log(1+exp(-b['rho']))
    star=b['Lambda']+b['zeta_bound']+ls(log(2),*[
        log(p(2)**(j-1)/j)-(j-1)*b['r0'] for j in range(2,8)])
    tiny=log(p('1e10'))-exp(b['Lambda'])+7*b['Lambda']+64*b['K']
    tau=log(p('27.5'))+b['Lambda'];B=exp(b['B']);lc=b['C']
    J=log(p('1e11'))+b['Lambda']+27*b['K']
    actual=ls(log(20)+b['Lambda']+8*b['K'],log(p('.9')*tau),log(p('55e8'))+b['Lambda']+26*b['K'])
    shape=201*B+log(3)+6*(log(204)+b['B'])
    core=exp(b['Lambda']+b['zeta_bound'])+log(3)+7*(log(11)+b['Lambda']+8*b['K'])+b['K']+log(2)/2
    target=p('.255')*p('1e-16')+p('3e-18')+p('3e-18')**2
    logs=dict(native_packet_N7=N7,star_log_C7=star,reference_F_C7=tiny,
        full_reference_pressure_C7=log(2000)+2*b['P'],
        reference_p1_C6_Xi=log(p('1e8'))+log(110)+b['Lambda']+26*b['K'],
        reference_ns_C6=log(p('1e5'))+3*b['K'],actual_log_CE_C6=actual,J6=J,
        core_Cscaled_regular_C6=core,shape_CE_C6=shape,
        normalized_actual_base_target_C6=log(p('2e6'))+412*B-2*lc)
    checks=dict(retained_core_parameter_checks=all(r['checks'].values()),
        retained_reference_C5_inequalities=all(ref.parameters(digits)['checks'].values()),
        native_seventh_packet_below_same_Kold=hi(N7-b['K'])<0,
        star_log_C7_below_10LambdaK9=hi(star-log(10)-b['Lambda']-9*b['K'])<0,
        reference_F_C7_below_point01=hi(tiny)<lo(-log(100)),
        reference_E_C7_below_one=hi(tiny+log(220)/2)<0,
        full_axis_and_reference_pressure_C7_below_Kold=hi(logs['full_reference_pressure_C7']-b['K'])<0,
        reference_U_C7_below_Kold=hi(ls(log(9),b['K']-b['Lambda'])-b['K'])<0,
        reference_p1_C6_fits_same_B_over100=hi(logs['reference_p1_C6_Xi']-b['B']+log(100))<0,
        reference_ns_C6_fits_same_B_over100=hi(logs['reference_ns_C6']-b['B']+log(100))<0,
        actual_prescription_log_CE_C6_below_J6=hi(actual-J)<0,
        actual_J6_below_same_B_over100=hi(J-b['B']+log(100))<0,
        actual_native_g_tail_below_j_over2=hi(b['K']-b['Lambda'])<lo(log(p('5e-19'))),
        actual_width_g_tail_below_j_over2=hi(log(p('1e8'))+3*b['K']-20*lc)<lo(log(p('5e-19'))),
        sixth_shape_tame_sum_fits_old_Eb=hi(shape-204*B)<0,
        old_Eb_still_below_one=hi(204*B-lc)<0,
        core_regular_scaled_C6_below_expB=hi(core-B)<0,
        actual_Cp_C6_log_length_fits_412B=hi(log((tau+2000*B)/2+2)-4*B)<0,
        same_C6_base_target_below_epsilon_over4=hi(logs['normalized_actual_base_target_C6'])<lo(log(p('1e-16')/4)),
        same_C6_root_target_below_d6=hi(target)<p('2.86e-17'),
        same_C6_root_contraction=hi(8*1000**2*100*target)<1,
        same_C6_root_coefficients_below_5p72e14=hi(2000*target)<p('5.72e-14'),
        eta_dependent_Q_C6_norm_below_one=hi(25*exp(p('1.2'))/256)<1,
        actual_first_slope_U_C5_below_B=hi(log(p('55e5'))+3*b['K']-b['B'])<0,
        B_and_P_below_actual_C=hi(b['B']-lc)<0 and hi(b['P']-lc)<0,
        fixed_width_first_radial_cost_fits_C20=hi(154*log(10)-20*log(p('1e9')))<0,
        original_compact_Xa_floor_retained=hi(-log(4)+b['Lambda']-4096*lc)<0)
    return dict(interval_digits=digits,log_bounds={k:dict(lower=mp.nstr(lo(v),125),upper=mp.nstr(hi(v),125)) for k,v in logs.items()},
        checks={k:bool(v) for k,v in checks.items()},actual_C_G_H_N_Z5_or_FINAL_profile_materialized=False)

def envelope(A=2**16,P=16):
    if any(type(v) not in (int,Q) for v in (A,P)):raise ValueError('Exact envelope parameters required')
    A,P=Q(A),Q(P)
    if A<2**16 or not 16<=P<=A:raise ValueError('Outside native envelope')
    Ecover=1+28*P+29*P+1000*A;Ucover=36
    inv=6*A**16;W=1+28*A**9
    Ns=A**11+45*A**9+32*A**8
    Qs=A**10+43*A**27
    p1=2*A*Qs;p2=2*A**18*Ns
    R=1+2*A**23+A**2+A**31
    return dict(A=str(A),P=str(P),R5_bound_power=32,pressure_C5_bound_power=31,
        native_field_C6_bound_powers=[2,1],native_first_radial_C5_bound_power=5,
        full_moment_C6_bound_power=8,inverse_E_C5_bound_power=17,shear_C5_bound_power=23,
        actual_variance_cutoff_or_FINAL_gate_verified=False,
        checks=dict(sum_of_all_native_E_C6_regions_fits_A2=Ecover<A**2,
            sum_of_all_native_U_C6_regions_fits_A=Ucover<A,
            native_first_radial_region_SUM_fits_A5=A**3+3*A**2+A<A**5,
            inverse_E_C5_fits_A17=inv<A**17,
            native_shear_C5_fits_A23=1+2*A**22<A**23,
            six_normalized_moment_derivative_cost_W_fits_A10=W<A**10,
            complete_Ns_with_axial_pressure_fits_A12=Ns<A**12,
            complete_Qs_fits_A29=Qs<A**29,
            complete_p1_C5_fits_A31=p1<A**31,
            complete_p2_C5_fits_A31=p2<A**31,
            actual_native_R5_common_sum_fits_A32=R<A**32,
            all_five_moments_plus_regular_Cp_and_axis_fit_A8=2*A**5+A**6+A**3<A**8))

@lru_cache(None)
def algebra():
    eta,y=s.symbols('eta y',real=True)
    ell=1+eta/7+eta**6*(1+y)/13+eta**7*y*y/17
    psi=1-3*y*y+2*y**3
    actual=ell.subs(y,0)+s.integrate(psi*s.diff(ell,y),(y,0,y))
    average=psi*ell-s.integrate(s.diff(psi,y)*ell,(y,0,y))
    f_inv=1+eta*eta;f_inv2=f_inv**2
    norm=lambda v:sum(abs(s.diff(v,eta,j).subs(eta,1))/factorial(j) for j in range(7))
    noncommuting=s.expand(s.diff((psi+eta*y/11)*s.diff(ell,y),eta)-(psi+eta*y/11)*s.diff(ell,y,eta))
    return dict(axis_C7_coefficient_sum=5*comb(11,4),normalizer_C6_costs=[str(norm(f_inv)),str(norm(f_inv2))],
        checks=dict(seventh_coefficient_geometric_factor_below_ten=Q(4,3)**8<10,
            sixth_differentiated_n_factor_below_eighteen=Q(7,4)*Q(4,3)**8<18,
            native_log_inverse_C7_constants_fit_N7=10**7*11**7<10**30,
            native_mixed_C6_constant_fits_N7=7*5**7*11**6*18<10**30,
            axis_seventh_positive_mixture_cost1650=5*comb(11,4)==1650,
            positive_taper_all_eight_angular_jets_exact=all(s.expand(s.diff(actual-average,eta,j))==0 for j in range(8)),
            positive_taper_mass_one=s.integrate(-s.diff(psi,y),(y,0,1))==1,
            positive_taper_weight_nonnegative=s.expand(-s.diff(psi,y)-6*y*(1-y))==0,
            eta_dependent_taper_would_not_commute=noncommuting!=0,
            seventh_finite_exponential_sum_below_three=sum((Q(1,factorial(j)) for j in range(8)),Q())<3,
            T7_constant_covers_all_seven_exp_jets=3*11**7<10**10,
            sixth_angular_derivative_cost_is_seven=factorial(7)//factorial(6)==7,
            reference_W6_constant40_covers_all_terms=1+4+4*7<40,
            reference_Sq6_constant_covers_all_terms=400+5+500<10**4,
            reference_Sn6_constant_covers_both_pressure_terms=200+5+35+28+8+2<1000,
            reference_p1_C6_variation_constant_fits=225+900000<10**8,
            reference_ns_C6_variation_constant_fits=10+2000<10**5,
            sixth_shape_log_f_still_below_six=Q(7,10)+sum((Q(2,j) for j in range(1,7)),Q())<6,
            sixth_shape_exp_sum_still_below1000=sum((Q(6)**j/factorial(j) for j in range(7)),Q())<1000,
            sixth_f_majorant_sum28=sum(range(1,8))==28,
            normalizer_C6_costs_still5_and25=norm(f_inv)==5 and norm(f_inv2)==25,
            annular_energy_constant_still8125=25*(101+160+64)==8125,
            native_pressure_fifth_uses_moment_sixth=factorial(6)//factorial(5)==6))

@lru_cache(None)
def root_review():
    e=s.Symbol('eta',real=True);q=(1+e*e)**2/31
    c=s.Rational(1,1000)*(1+e/7+e**2/11+e**3/13+e**4/17+e**5/19+e**6/23)
    z=7*c+q*c*c;rows=[]
    for center in (-1,s.Rational(1,3),1):
        cj=[s.diff(c,e,j).subs(e,center)/factorial(j) for j in range(7)]
        qj=[s.diff(q,e,j).subs(e,center)/factorial(j) for j in range(7)]
        lower=sum(qj[i]*cj[j]*cj[k] for i in range(7) for j in range(7) for k in range(7)
            if i+j+k==6 and (i,j,k) not in ((0,6,0),(0,0,6)))
        frozen=sum(qj[0]*cj[j]*cj[6-j] for j in range(1,6))
        rows.append(dict(eta=str(center),exact=bool(s.diff(z,e,6).subs(e,center)/720==(7+2*qj[0]*cj[0])*cj[6]+lower),
            frozen_wrong=bool(s.expand(lower-frozen)!=0),Q4_c1_c1_nonzero=bool(qj[4]*cj[1]**2!=0),higher_Q_zero=bool(qj[5]==qj[6]==0)))
    return dict(manufactured_not_actual_root=True,rows=rows,checks=dict(
        sixth_full_eta_dependent_root_recurrence_exact=all(r['exact'] for r in rows),
        freezing_earlier_Q_eta_detected_at_sixth=all(r['frozen_wrong'] for r in rows),
        Q4_c1_c1_term_present_at_sixth=all(r['Q4_c1_c1_nonzero'] for r in rows),
        no_unneeded_fifth_sixth_Q_derivatives=all(r['higher_Q_zero'] for r in rows)))

@lru_cache(None)
def flow_review():
    cc=module('native6_independent_rational_flows',HERE.parent/'prescribed_continuation_foundation/continuation_checks.py')
    J=cc.modules()[0].Jet;refs=[];stages=[]
    for center in (Q(-1),Q(1,3),Q(1)):
        order=8;y=J.variable(0,0,n=2,order=order);eta=J.variable(1,center,n=2,order=order)
        X=2*cc.exponential(y);root=2*cc.exponential(y/2)
        f0=1+eta/7+eta**7/19;f1=(1+eta**6)/23
        u0=4*eta+Q(1,10)+eta**7/29;u1=(1+eta**6)/31
        F=f0+X*f1;Er=root*F;Ur=u0+X*u1
        mom=[X*u0+X**2*u1/2,X**2*f0+2*X**3*f1/3,
            X**2*u0*f0+2*X**3*(u0*f1+u1*f0)/3+X**4*u1*f1/2,
            X*u0*u0+X**2*u0*u1+X**3*u1*u1/3-X**2*f0*f0/2-2*X**3*f0*f1/3-X**4*f1*f1/4,
            X*f0*f0+X**2*f0*f1+X**3*f1*f1/3]
        axis=-2+eta/11+eta**7/37;cc.coordinates.center=center
        st=cc.coordinates(X,root,Er,Ur,mom,axis)
        h=Q(1,200);D=Q(1,2)-h;Ap=Q(1,2)+h;d=1-eta*eta;L=1-2*h*eta*eta
        lr=1+F.derivative(0)/F;hc=D*eta+d*Ur
        sq=-st['W']*lr-h*(1-2*eta*Ur)-hc*Er.derivative(1)/Er
        pi=st['Pi'];sn=-st['W']*Ur.derivative(0)-Ap*(1-2*eta*Ur)*Ur-hc*Ur.derivative(1)-d*pi.derivative(1)+4*Ap*eta*pi+eta*Er*Er
        ep=st['p1'].derivative(0)+lr*st['p1']-X*sq/L
        en=st['ns'].derivative(0)+st['ns']-sn/L
        refs.append(dict(eta=str(center),reference_sources_angular6_exact=all(ep.partial((0,j))==en.partial((0,j))==0 for j in range(7)),
            nonzero_axis_seventh=axis.partial((0,7))!=0,
            dropped_axis_sixth_detected=(4*Ap*eta*axis-d*axis.derivative(1)).partial((0,6))!=0,
            dropped_Cp_sixth_detected=(4*Ap*eta*mom[4]-d*mom[4].derivative(1)).partial((0,6))!=0))
        for name,kappa,beta,theta,stop in [('activation',1-y/4,1,0,False),('axial_taper',Q(1,10),1-y/3,0,False),('shear_blend',Q(1,10),1,y/4,True)]:
            a=(1-theta)*kappa*st['p1']+Q(4,5)*theta
            b=cc.jet(0,order) if stop else -kappa*beta*X*st['ns']/2
            ebase=cc.jet({(0,j):v for (i,j),v in Er.terms.items() if i==0},order)
            ubase=cc.jet({(0,j):v for (i,j),v in Ur.terms.items() if i==0},order)
            E=ebase*cc.exponential(y/2-cc.primitive(a)/2);U=ubase+cc.primitive(b)
            incoming=[cc.jet({(0,j):v for (i,j),v in m.terms.items() if i==0},order) for m in mom]
            actual=[m+cc.primitive(v) for m,v in zip(incoming,cc.densities(X,root,E,U))]
            states=cc.coordinates(X,root,E,U,actual,axis)
            stages.append(dict(eta=str(center),stage=name,
                both_shear_equations_total6_exact=cc.zero_through(E.derivative(0)-(1-a)*E/2,6) and cc.zero_through(U.derivative(0)-b,6),
                five_actual_density_flows_total6_exact=all(cc.zero_through(m.derivative(0)-v,6) for m,v in zip(actual,cc.densities(X,root,E,U))),
                full_actual_pressure_flow_total6_exact=cc.zero_through(states['Pi'].derivative(0)-E*E/2,6),
                all_sixth_incoming_moment_traces_nonzero=all(m.partial((0,6))!=0 for m in incoming),
                actual_reference_moments_not_reused=any(not cc.zero_through(m-v,6) for m,v in zip(actual,mom))))
    return dict(manufactured_not_selected_native_profile=True,reference_rows=refs,actual_stage_rows=stages,
        checks=dict(seventh_input_sixth_reference_source_identity=all(r['reference_sources_angular6_exact'] for r in refs),
            seventh_axis_trace_retained=all(r['nonzero_axis_seventh'] for r in refs),
            deleting_axis_sixth_source_detected=all(r['dropped_axis_sixth_detected'] for r in refs),
            deleting_Cp_sixth_source_detected=all(r['dropped_Cp_sixth_detected'] for r in refs),
            actual_all_three_prescription_stage_shears_total6_exact=all(r['both_shear_equations_total6_exact'] for r in stages),
            actual_all_five_density_flows_total6_exact=all(r['five_actual_density_flows_total6_exact'] for r in stages),
            actual_full_pressure_flow_total6_exact=all(r['full_actual_pressure_flow_total6_exact'] for r in stages),
            six_order_incoming_moment_traces_not_zeroed=all(r['all_sixth_incoming_moment_traces_nonzero'] for r in stages),
            actual_moments_not_copied_from_reference=all(r['actual_reference_moments_not_reused'] for r in stages)))

@lru_cache(None)
def review():
    parts=dict(algebra=algebra(),root=root_review(),flows=flow_review(),envelope=envelope())
    for r in parts.values():r['checks']={k:bool(v) for k,v in r['checks'].items()}
    parts['parameters']=[parameters(d) for d in (80,110)]
    checks={k:v for r in parts.values() if isinstance(r,dict) for k,v in r['checks'].items()}
    checks['all_outward_Md64_comparisons']=all(all(r['checks'].values()) for r in parts['parameters'])
    parts['checks']=checks;return parts
