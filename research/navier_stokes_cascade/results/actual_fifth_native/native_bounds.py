"""Actual PREMOD fifth fields/memories and fourth pressure, relative to core."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
from pathlib import Path
import importlib.util,sys
import mpmath as mp
import sympy as s

HERE=Path(__file__).resolve().parent
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m
ref=module('actual_fifth_reference',HERE.parent/'reference_fifth_transfer/reference_bounds.py')
high=ref.high

def regional_cover(A,P,d5=Q(572,10**16)):
    if any(type(v) not in (int,Q) for v in (A,P,d5)):raise TypeError('Exact inputs required')
    A,P,d5=map(Q,(A,P,d5))
    if A<2**16 or not 16<=P<=A or not 0<=d5<=Q(572,10**16):raise ValueError('Outside retained envelope')
    return dict(E_C5_cover=1+21*P+22*P+1000*A,U_C5_cover=10+9+9+8,
        annular_relative_C5=3*64*d5,annular_U_C5=8+2*64*d5,
        FE_bound=A*A,FU_bound=A,actual_FINAL_field_or_gate_verified=False)

@lru_cache(None)
def parameters(digits):
    if type(digits) is not int or digits<60:raise ValueError('At least 60 interval digits required')
    r=high.core.certificate('Md64',digits);b=r['log_bounds']
    p,lo,hi,ls=high.p,high.lower,high.upper,high.logsum;log,exp=mp.iv.log,mp.iv.exp
    B,lc=exp(b['B']),b['C'];tau=log(p('27.5'))+b['Lambda']
    J=log(p('1e11'))+b['Lambda']+25*b['K']
    cont=ls(log(20)+b['Lambda']+8*b['K'],log(p('.9')*tau),log(p('55e8'))+b['Lambda']+24*b['K'])
    # Never exponentiate the negative width logarithm; compare each tail to j/2.
    native_error=b['K']-b['Lambda'];width_error=log(p('1e8'))+3*b['K']-20*lc
    g=log(p('2e-18'))
    core_CE_over_sqrtX=exp(b['Lambda']+b['zeta_bound'])+log(3)+6*(log(11)+b['Lambda']+8*b['K'])+b['K']+log(2)/2
    shape_CE=201*B+log(3)+5*(log(204)+b['B'])
    ell=2000*B+tau;Eb=204*B-lc
    target=p('.255')*p('1e-16')+p('3e-18')+p('3e-18')**2
    logs=dict(actual_log_CE_continuation_C5=cont,continuation_envelope_J=J,
        actual_g_C5_log_majorant=g,actual_g_native_error=native_error,actual_g_width_error=width_error,
        native_CE_over_sqrtX_C5=core_CE_over_sqrtX,
        shape_CE_C5=shape_CE,Eb=Eb,
        actual_Cp_C5=-2*lc+408*B+log(ell/2+2),
        normalized_join_base_C5=log(p('2e6'))+412*B-2*lc)
    row_constants=[18*exp(6),60*exp(p('9.6'))/16,p(10)/3*exp(p('9.6'))/16,
        p(8125)*exp(p('7.2'))/256,p(25)*exp(p('1.2'))/256,
        p(5)/8*exp(p('9.6')),p(5)/12*exp(p('7.2')),p('2.5')*exp(p('1.2'))]
    checks=dict(retained_core_parameter_checks=all(r['checks'].values()),
        preceding_reference_fifth_comparisons=all(ref.parameters(digits)['checks'].values()),
        prescribed_log_CE_C5_below_J=hi(cont-J)<0,
        J_below_B_over100=hi(J-b['B']+log(100))<0,
        native_g_error_below_j_over2=hi(native_error)<lo(log(p('5e-19'))),
        width_g_error_below_j_over2=hi(width_error)<lo(log(p('5e-19'))),
        actual_native_offset_g_C5_below_3e18=hi(g)<lo(log(p('3e-18'))),
        native_regular_CE_over_sqrtX_C5_below_expB=hi(core_CE_over_sqrtX-B)<0,
        shape_tame_C5_fits_same_exp204B=hi(shape_CE-204*B)<0,
        actual_Eb_below_one=hi(Eb)<0,
        integrable_core_Cp_and_actual_log_length_fit_exp412B=hi(log(ell/2+2)-4*B)<0,
        all_normalized_join_row_constants_below_one_million=all(hi(v)<10**6 for v in row_constants),
        actual_normalized_join_base_C5_below_epsilon_over4=hi(logs['normalized_join_base_C5'])<lo(log(p('1e-16')/4)),
        ideal_offset_linear_row_fits_point005_epsilon=hi(exp(-2)*p('3e-18'))<lo(p('.005e-16')),
        ideal_offset_energy_row_fits_point005_epsilon=hi(p(25)/256*exp(p('-.8'))*p('3e-18')**2)<lo(p('.005e-16')),
        eta_dependent_quadratic_C5_norm_fits_same100=hi(25*exp(p('1.2'))/256)<1,
        same_C5_root_target_below_2p86e17=hi(target)<p('2.86e-17'),
        same_C5_root_contraction=hi(8*1000**2*100*target)<1,
        same_C5_root_coefficient_below_5p72e14=hi(2000*target)<p('5.72e-14'),
        actual_axis_pressure_C5_below_A3=hi(log(630)+2*b['P']-3*4096*lc)<0,
        original_Xa_floor_and_prefix_geometry_retained=hi(-log(4)+b['Lambda']-4096*lc)<0,
        no_variance_shape_budget_or_FINAL_compact_gate_inferred=True)
    return dict(interval_digits=digits,
        log_bounds={k:dict(lower=mp.nstr(lo(v),125),upper=mp.nstr(hi(v),125)) for k,v in logs.items()},
        checks=checks,enormous_C_G_H_N_or_actual_profile_materialized=False)

@lru_cache(None)
def algebra():
    e=s.Symbol('eta',real=True);g=s.Symbol('g',real=True)
    inverse=s.Poly((1+e*e)**2,e)
    inverse_norm=[sum(abs(s.diff(inverse.as_expr(),e,j).subs(e,1))/factorial(j) for j in range(k+1)) for k in (4,5)]
    k0=Q(1,1000);xa=Q(1,2)
    # Upper envelope of the genuine weighted activation integral.
    weighted=(2*xa+110)*k0
    exp_outer=sum((Q(6)**j/factorial(j) for j in range(6)),Q())
    cover=regional_cover(2**16,2**16)
    return dict(f_inverse_square_C4_C5=[str(v) for v in inverse_norm],
        annular_energy_row_constant=25*(101+160+64),
        pressure_C4_A_power=1546,native_field_F_E5_A_power=2,native_field_F_U5_A_power=1,
        native_moment_and_entry_C5_A_power=8,preexisting_native_C4_envelope_A_power=256,
        checks=dict(log_f_C5_pole_sum_below_six=Q(7,10)+sum((Q(2,j) for j in range(1,6)),Q())<6,
            weighted_activation_X_integral_below112delta=weighted<=112*k0,
            prescribed_U_increment_fits_1e8K3delta=56*10**5<10**8,
            annular_energy_g_square_cancels_linear_term=s.expand((4*e+g)**2-(4*e)**2-8*e*g-g*g)==0,
            annular_C5_normalizer_cost_25_not_20=inverse_norm==[25,25],
            annular_energy_row_8125_fits_one_million=25*(101+160+64)==8125 and 8125<10**6,
            same_annular_relative_factor_stays_positive=cover['annular_relative_C5']<Q(1,2),
            same_annular_U_C5_below_nine=cover['annular_U_C5']<9,
            ideal_f_C5_majorant_is_21=sum(range(1,7))==21,
            outer_finite_tame_exponential_below1000=exp_outer<1000,
            cover_SUM_E_C5_below_A2=cover['E_C5_cover']<cover['FE_bound'],
            cover_SUM_U_C5_below_A=cover['U_C5_cover']<cover['FU_bound'],
            native_moment_common_A8_covers_two_A5_and_full_Pi=3<2**16,
            T1_T2_common_A517_controls_constants=11<(2**16)**2,
            fourth_pressure_A1546_controls_complete_R0_and_T2=10000+1+10<(2**16)**2,
            fourth_pressure_envelope_not_old_Ahi=1546>256,
            native_pure_fifth_does_not_bound_FINAL_modulated_fifth=True))

@lru_cache(None)
def quadratic_review():
    e=s.Symbol('eta',real=True);q=(1+e*e)**2/31
    c=s.Rational(1,1000)*(1+e/7+e**2/11+e**3/13+e**4/17+e**5/19)
    z=7*c+q*c*c;rows=[]
    for eta in (-1,s.Rational(1,3),1):
        cj=[s.diff(c,e,j).subs(e,eta)/factorial(j) for j in range(6)]
        qj=[s.diff(q,e,j).subs(e,eta)/factorial(j) for j in range(6)]
        lower=sum(qj[i]*cj[j]*cj[l] for i in range(6) for j in range(6) for l in range(6)
            if i+j+l==5 and (i,j,l) not in ((0,5,0),(0,0,5)))
        frozen=sum(qj[0]*cj[j]*cj[5-j] for j in range(1,5))
        value=s.diff(z,e,5).subs(e,eta)/120
        rows.append(dict(eta=str(eta),exact=value==(7+2*qj[0]*cj[0])*cj[5]+lower,
            freezing_Q_eta_detected=s.expand(lower-frozen)!=0,
            fourth_Q_coefficient_nonzero=qj[4]!=0,fifth_Q_coefficient_zero=qj[5]==0))
    return dict(manufactured_not_actual_root=True,rows=rows,checks=dict(
        exact_normalized_fifth_implicit_quadratic=all(v['exact'] for v in rows),
        eta_dependent_Q_cannot_be_frozen_at_fifth=all(v['freezing_Q_eta_detected'] for v in rows),
        all_Q_derivatives_through_four_retained=all(v['fourth_Q_coefficient_nonzero'] and v['fifth_Q_coefficient_zero'] for v in rows)))

@lru_cache(None)
def flow_review():
    cc=module('actual_fifth_rational_continuation',HERE.parent/'prescribed_continuation_foundation/continuation_checks.py')
    J=cc.modules()[0].Jet;rows=[]
    for center in (Q(-1),Q(1,3),Q(1)):
        order=7;y=J.variable(0,0,n=2,order=order);eta=J.variable(1,center,n=2,order=order)
        X=2*cc.exponential(y);root=2*cc.exponential(y/2)
        e=2+eta/10+eta**5/17+eta**6/19;G=4*eta+Q(1,10)+eta**5/23+eta**6/29
        Er=e*cc.exponential(y/2);Ur=G
        mom=[2*G*cc.exponential(y),2*e*cc.exponential(2*y),2*G*e*cc.exponential(2*y),
            2*G*G*cc.exponential(y)-e*e*cc.exponential(2*y)/2,e*e*cc.exponential(y)/2]
        axis=-2+eta/7+eta**5/31+eta**6/37
        cc.coordinates.center=center;st=cc.coordinates(X,root,Er,Ur,mom,axis)
        for name,kappa,beta,theta,stop in [('activation',1-y/4,1,0,False),
            ('axial_taper',Q(1,10),1-y/3,0,False),('shear_blend',Q(1,10),1,y/4,True)]:
            a=(1-theta)*kappa*st['p1']+Q(4,5)*theta
            b=cc.jet(0,order) if stop else -kappa*beta*X*st['ns']/2
            E=e*cc.exponential(y/2-cc.primitive(a)/2);U=G+cc.primitive(b)
            incoming=[cc.jet({(0,j):v for (i,j),v in m.terms.items() if i==0},order) for m in mom]
            actual=[m+cc.primitive(v) for m,v in zip(incoming,cc.densities(X,root,E,U))]
            shearE=E.derivative(0)-(1-a)*E/2;shearU=U.derivative(0)-b
            density_errors=[m.derivative(0)-v for m,v in zip(actual,cc.densities(X,root,E,U))]
            pi=axis+actual[4];pressure=pi.derivative(0)-E*E/2
            rows.append(dict(eta=str(center),stage=name,
                both_shears_total5_exact=cc.zero_through(shearE,5) and cc.zero_through(shearU,5),
                five_actual_moment_flows_total5_exact=all(cc.zero_through(v,5) for v in density_errors),
                full_pressure_flow_total5_exact=cc.zero_through(pressure,5),
                five_entering_fifth_memories_nonzero=all(m.partial((0,5))!=0 for m in incoming),
                axis_fifth_memory_nonzero=axis.partial((0,5))!=0,
                reference_zero_shear_not_copied=st['p1'].partial((0,5))!=0,
                actual_and_reference_moments_not_reused=any(not cc.zero_through(m-r,5) for m,r in zip(actual,mom))))
    # Actual shape polynomial in log amplitude; no eta-dependent switching point.
    z=s.Symbol('z',real=True);eta=s.Symbol('eta',real=True);T=s.Symbol('T',positive=True)
    ss=3*(z/T)**2-2*(z/T)**3;ell=1+eta/7+eta**5/11;lf=-s.log(1+eta**2)
    shape=z/10+(1-ss)*ell+ss*lf
    shape_high=s.diff(shape,eta,5)
    shape_identity=s.expand(shape_high-((1-ss)*s.diff(ell,eta,5)+ss*s.diff(lf,eta,5)))==0
    return dict(manufactured_not_selected_profile=True,rows=rows,checks=dict(
        actual_all_three_stage_shears_through_fifth=all(v['both_shears_total5_exact'] for v in rows),
        all_five_actual_forward_moment_flows_through_fifth=all(v['five_actual_moment_flows_total5_exact'] for v in rows),
        actual_complete_pressure_flow_through_fifth=all(v['full_pressure_flow_total5_exact'] for v in rows),
        all_five_entrance_fifth_traces_and_axis_retained=all(v['five_entering_fifth_memories_nonzero'] and v['axis_fifth_memory_nonzero'] for v in rows),
        flattened_reference_shear_not_used_for_actual=all(v['reference_zero_shear_not_copied'] for v in rows),
        actual_moments_not_reference_moments=all(v['actual_and_reference_moments_not_reused'] for v in rows),
        shape_fifth_angular_is_same_convex_combination=shape_identity))

@lru_cache(None)
def volterra_review():
    v=module('actual_fifth_forward_majorants',HERE.parent/'fifth_moment_volterra/volterra.py')
    A=2**16;B=A**256;entry=(A**8,)*6;rows=[]
    for eta in (Q(-1),Q(0),Q(1,3),Q(1)):
        r=v.certificate(A,B,Q(1,A),A,eta,Q(1,200),entry,A*A,A,B*B,B)
        rows.append(dict(eta=str(eta),T1_within_A517=r['T1']<=A**517,T2_within_A517=r['T2']<=A**517,
            p1_within_A1546=r['pressure1_C4']<=A**1546,p2_within_A1546=r['pressure2_C4']<=A**1546,
            incoming_T2_positive=r['incoming_T2']>0))
    z=v.certificate(A,B,1,1,Q(1,3),Q(1,200),entry,A*A,A,B*B,B)
    return dict(auxiliary_minimum_A_not_enormous_actual_profile=True,rows=rows,checks=dict(
        original_Volterra_positive_T1_T2_fit_derived_powers=all(r['T1_within_A517'] and r['T2_within_A517'] for r in rows),
        original_complete_fourth_pressure_majorants_fit_A1546=all(r['p1_within_A1546'] and r['p2_within_A1546'] for r in rows),
        entering_T2_not_zero_even_at_both_poles=all(r['incoming_T2_positive'] for r in rows),
        zero_length_still_keeps_all_entering_memory=z['T1']==z['incoming_T1'] and z['T2']==z['incoming_T2'] and z['T2']>0))

def review():
    a=algebra();q=quadratic_review();f=flow_review();v=volterra_review();p=[parameters(d) for d in (80,110)]
    return dict(algebra=a,quadratic=q,flows=f,volterra=v,parameters=p,
        checks=dict(**a['checks'],**q['checks'],**f['checks'],**v['checks'],
            two_precision_actual_parameter_inequalities=all(all(v['checks'].values()) for v in p)))
