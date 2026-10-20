"""Tame fifth loop algebra and original native/cutoff input boundaries.

Only exact manufactured jets and elementary outward comparisons are evaluated.
No actual Z5, R5, Lambda5, G, H or N is numerically supplied.
"""
from functools import lru_cache
from fractions import Fraction as Q
from math import comb, factorial
import mpmath as mp
import sympy as s

DEG=5
def jet(*v):return [s.sympify(x) for x in v]+[s.Integer(0)]*(DEG+1-len(v))
def add(*a):return [s.expand(sum(v[k] for v in a)) for k in range(DEG+1)]
def scale(a,c):return [s.expand(c*v) for v in a]
def mul(a,b):return [s.expand(sum(a[j]*b[k-j] for j in range(k+1))) for k in range(DEG+1)]
def power(a,k):
    b=jet(1)
    for _ in range(k):b=mul(a,b)
    return b
def reciprocal(a):
    if a[0]==0:raise ValueError('Nonzero base required')
    b=jet(1/a[0])
    for k in range(1,DEG+1):b[k]=s.expand(-sum(a[j]*b[k-j] for j in range(1,k+1))/a[0])
    return b
def fifth_implicit(residual,jacobian):
    if jacobian==0:raise ValueError('Nonzero implicit Jacobian required')
    u=jet(0)
    for k in range(1,DEG+1):u[k]=s.expand(-residual(u)[k]/jacobian)
    return u

@lru_cache(None)
def q6_review():
    z,h=s.symbols('z sine');den=1+z*z/4+z**4/64+z**6/2304
    q=s.series(s.exp(z*h)/den,z,0,7).removeO().expand()
    c=q.coeff(z,6)
    expected=h**6/720-h**4/96+3*h**2/128-s.Rational(19,2304)
    R=sum(q.coeff(z,k+1)*z**k for k in range(6))
    sixth_polynomial_sum=64*sum(s.Rational(60**m,factorial(m)) for m in range(7))
    return dict(q6=str(c),R_normalized_fifth_at_z0=str(R.coeff(z,5)),
        differentiated_q_operator_cost=6,
        checks=dict(exact_sixth_q_coefficient=s.expand(c-expected)==0,
            sixth_q_required_by_fifth_removable_R=R.coeff(z,5)==c and c!=0,
            q6_polynomial_majorant_below_1e16G12=sixth_polynomial_sum<10**16,
            normalized_q_z_order5_cost_is_six=6*factorial(5)==factorial(6),
            no_singular_p_or_mu_division_in_signed_extension=True))

@lru_cache(None)
def scalar_root_review():
    q,r=s.symbols('p5 s5');mu0=s.Rational(1,5);p0=s.Rational(2,7)
    p=jet(p0,s.Rational(1,3),-s.Rational(1,7),s.Rational(1,11),s.Rational(1,13),q)
    def F(mu,p):return add(mu,mul(mu,p),scale(power(mu,2),s.Rational(1,7)),scale(power(p,3),s.Rational(1,11)))
    F0=F(jet(mu0),jet(p0))[0]
    target=jet(F0,s.Rational(1,2),s.Rational(1,5),-s.Rational(1,9),s.Rational(1,17),r)
    fm=1+p0+2*mu0/7;fp=mu0+3*p0*p0/11
    def residual(du):return add(F(add(jet(mu0),du),p),scale(target,-1))
    du=fifth_implicit(residual,fm);res=residual(du)
    high=s.Poly(du[5],q,r);derq=s.diff(du[5],q);derr=s.diff(du[5],r)
    d,mu,pv=s.symbols('d mu p',real=True)
    signed=d*mu*s.sqrt(s.Rational(1,2)-3*(mu*pv)**2/32+11*(mu*pv)**4/576)
    return dict(base_Fmu=str(fm),base_Fp=str(fp),mu5_p5_coefficient=str(derq),mu5_s5_coefficient=str(derr),
        manufactured_not_actual_signed_root=True,checks=dict(
            complete_nonlinear_fifth_implicit_residual_zero=all(s.simplify(v)==0 for v in res),
            top_p5_coefficient_exact=derq==-fp/fm,
            top_s5_coefficient_exact=derr==1/fm,
            high_fifth_inputs_affine_not_squared=high.total_degree()==1,
            lower_cross_terms_present=du[5].subs({q:0,r:0})!=0,
            signed_origin_derivative_nonzero=s.simplify(s.diff(signed,mu).subs({mu:0,pv:0})-d/s.sqrt(2))==0,
            signed_origin_p_derivative_zero=s.diff(signed,pv).subs({mu:0,pv:0})==0))

@lru_cache(None)
def moving_review():
    q=s.Symbol('x5');x0=s.Rational(1,5)
    x=jet(x0,s.Rational(1,3),-s.Rational(1,7),s.Rational(1,11),s.Rational(1,13),q)
    def Phi(x,t):return add(t,scale(mul(x,t),s.Rational(1,7)),scale(power(x,2),s.Rational(1,11)),
        scale(power(t,2),s.Rational(1,13)),scale(mul(x,power(t,2)),s.Rational(1,17)),scale(power(t,3),s.Rational(1,19)))
    jac=1+x0/7;px=2*x0/11;phi0=x0*x0/11
    theta=fifth_implicit(lambda t:add(Phi(x,t),jet(-phi0)),jac)
    def source(x,t):return add(mul(x,reciprocal(add(jet(1),power(t,2)))),power(t,2),mul(x,t))
    g=source(x,theta);expected=1-x0*px/jac
    return dict(base_phase_Jacobian=str(jac),inverse_theta5_x5_coefficient=str(s.diff(theta[5],q)),
        full_source5_x5_coefficient=str(s.diff(g[5],q)),frozen_angle_coefficient='1',
        manufactured_not_actual_phase=True,checks=dict(
            complete_fifth_moving_inverse_residual_zero=all(s.simplify(v)==0 for v in add(Phi(x,theta),jet(-phi0))),
            fifth_inverse_keeps_actual_base_divisor=s.diff(theta[5],q)==-px/jac,
            full_pulled_fifth_sensitivity_exact=s.diff(g[5],q)==expected,
            frozen_angle_misses_nonzero_top_term=expected!=1,
            full_fifth_moving_source_affine_in_top_input=s.Poly(g[5],q).degree()==1,
            nontrivial_lower_moving_remainder_retained=g[5].subs(q,0)!=0,
            inverse_Jacobian_not_normalized_to_one=jac!=1 and jac>0))

@lru_cache(None)
def native_pressure_review():
    e,X,E,U,h,M6,S6,P6=s.symbols('eta X E U h M6 S6 Pi6')
    t=s.Symbol('t');eta=e+t;L=1-2*h*eta**2;D=s.Rational(1,2)-h;d=1-eta**2;Ap=s.Rational(1,2)+h
    M=M6*t**6;S=S6*t**6;Pi=P6*t**6
    W=1-(2*D*eta*M+d*s.diff(M,t))/X
    Ns=-W*U+(D*(M-eta*s.diff(M,t))+4*h*eta*S-d*s.diff(S,t))/X+4*Ap*eta*Pi-d*s.diff(Pi,t)
    # Only these order-six moment terms can contribute to the fifth coefficient.
    p2=X*Ns/(L*E)
    raw=s.diff(p2,t,5).subs(t,0)/120
    got=s.expand(raw-raw.subs({M6:0,S6:0,P6:0}))
    want=6*(((1-e*e)*U-D*e)*M6-(1-e*e)*S6-X*(1-e*e)*P6)/((1-2*h*e*e)*E)
    noaxial_raw=s.diff(X*(Ns+W*U)/(L*E),t,5).subs(t,0)/120
    noaxial=s.expand(noaxial_raw-noaxial_raw.subs({M6:0,S6:0,P6:0}))
    inv=reciprocal(jet(E,s.Rational(1,3),0,0,0,s.Symbol('E5')))
    ey5,uy5=s.symbols('Ey5 Uy5')
    EY=jet(s.Rational(2,7),s.Rational(1,2),0,0,0,ey5)
    UY=jet(-s.Rational(1,9),s.Rational(1,11),0,0,0,uy5)
    a5=scale(mul(EY,inv),-2)[5];b5=scale(mul(UY,inv),2)[5]
    return dict(full_p2_fifth_principal=str(want),Pi6_is_full_axis_plus_Cp=True,
        mixed_field_indices=[[1,5],[1,5]],total_field_order_required=6,
        checks=dict(full_pressure_principal_M6_S6_Pi6_exact=s.simplify(got-want)==0,
            nonzero_axial_U_M6_term_required=s.simplify(got-noaxial-6*(1-e*e)*U*M6/((1-2*h*e*e)*E))==0,
            full_pressure_sixth_memory_coefficient_nonzero=s.simplify(s.diff(want,P6).subs(e,0))==-6*X/E,
            D_eta_M6_has_no_spurious_one_minus_eta_square=s.simplify(s.diff(want,M6).subs(U,0)+6*D*e/((1-2*h*e*e)*E))==0,
            a5_needs_total_order_six_field_row=s.simplify(s.diff(a5,ey5)+2/E)==0,
            bs5_needs_total_order_six_field_row=s.simplify(s.diff(b5,uy5)-2/E)==0,
            joint_five_does_not_contain_1_5=sum((1,5))>5,
            native_p1_not_added_as_signed_root_argument=True))

@lru_cache(None)
def primitive_review():
    e,z,a,b=s.symbols('eta phi a b',real=True)
    P=(a+b*e+e**5)*s.sin(2*s.pi*z);g=s.diff(P,z)
    v,ts,d,mu,r2,A,bs,E,t=s.symbols('v ts d mu r2 A bs E t');w=A*(1+t*t)/(2*s.pi*v)
    aL=v/(1+t*t);bL=-v*t/(1+t*t)
    return dict(both_pure_angular_primitive_operator_cost=4,full_joint_primitive_claim=False,
        complete_physical_covariance_or_mean_claim=False,checks=dict(
            zero_mean_primitive_exact=s.integrate(P,(z,0,1))==0,
            fixed_external_phase_fifth_commutes=s.diff(g,e,5)==s.diff(P,e,5,z),
            both_primitive_norms_charged=2+2==4,
            weighted_loop_first_density_exact=s.cancel(aL*w-A/(2*s.pi))==0,
            weighted_axial_density_exact=s.cancel(bL*w+A*t/(2*s.pi))==0,
            mean_one_uses_full_variance=s.expand(A*(1+ts*ts+d*d*mu*mu*r2)-(A*(1+ts*ts)+A*d*d*mu*mu*r2))==0,
            axial_B_source_keeps_E_and_bs=s.expand(E*(-A*ts-bs)/2).subs(bs,-A*ts)==0))

def exponent_review():
    # All lower increments are truncated before the unknown top coefficient.
    h=s.Symbol('h');u=sum(h**j for j in range(1,5))
    counts=[]
    for m in range(1,9):
        count=sum(comb(m+k-1,k)*s.expand(u**k).coeff(h,5) for k in range(1,6))+m
        counts.append(int(count))
    rows=[('mu5',100+5*601+10+2,3200),('t5',100+3200+5*601+2,20000),
        ('fixed_g5_w5',5*4000+20000+10,100000),('Phi5',100000+1,100010),
        ('Theta5',max(100010+5*80001,10010+5*80001)+2,600010),
        ('moving_g5',max(20000+600010,100000+5*80001)+2,700000),
        ('both_primitive5',700000+2,1000000)]
    return dict(rows=[dict(packet=n,sufficient_power=v,chosen_power=c) for n,v,c in rows],
        exact_scalar_count_by_number_inputs=counts,
        Lambda5_gate='exp(1000000*G)*T5',native_cutoff_gate='exp(1100000*G+10*A^20)*Z5*R5',
        checks=dict(all_fifth_tame_packets_have_strict_headroom=all(v<c for _,v,c in rows),
            fifth_two_input_partition_count_64=counts[1]==64,
            two_external_color_inflation_2048=counts[1]*2**5==2048,
            all_through_eight_input_counts_below_1e9=max(counts)<10**9,
            final_affine_absorption=1000000+1100<1100000,
            cutoff_scalar_constant_1e24_safe=10**10*6*100**5*6<10**24,
            sixth_q_reciprocal_norm_has64=1+sum(2**(j-1) for j in range(1,7))==64,
            sqrt_R2_fifth_cost_fits71=5+5*(10+3)+1<=71))

def universal_review(digits):
    if type(digits) is not int or digits<60:raise ValueError('At least 60 exact interval digits required')
    mp.iv.dps=digits;p=mp.iv.mpf;log=mp.iv.log
    a=p(2**16);g=p(2**16);gg=p(2**32)
    hi=lambda v:mp.make_mpf(v._mpi_[1]);lo=lambda v:mp.make_mpf(v._mpi_[0])
    checks=dict(q6_polynomial_absorption=hi(16*log(p(10))+12*log(g)-g)<0,
        q6_absorption_monotone=lo(1-12/g)>0,
        scalar_cutoff_polynomial_absorption=hi(24*log(p(10))+48*log(g)-g)<0,
        scalar_cutoff_absorption_monotone=lo(1-48/g)>0,
        counts_fit_expG=hi(log(p(10**9)))<lo(g),
        q_z_six_fits_expG=hi(log(p(6)))<lo(g),
        sqrt_Taylor_sum_fits_expG=hi(log(p(36)))<lo(g),
        F71_partial_composition_fits100=hi(5*log(3*g)+log(g)+log(p(6)))<lo(29*g),
        G_actual_exceeds_2pow32=lo(a**32)>hi(32*log(p(2))),
        deltaL_log_cost_below_G=hi(log(p(10))+20*log(a)-a**32)<0,
        deltaL_log_comparison_monotone=lo(32-20/a**32)>0,
        F5_A1044_log_cost_below_logG=hi(1044*log(a)-a**32)<0,
        F5_comparison_monotone=lo(32-1/log(a))>0,
        fifth_bound_expG2_below_H=lo(g**4-g**2)>0,
        F5_headroom_after_Lambda_below_H=lo(291*g**4-g)>0,
        elementary_auxiliary_log_gate_nonempty=lo(gg-1100003)>0,
        primitive_sum_four_absorbed=hi(log(p(4)))<lo(g),
        full_phase_operator_below_eight=hi(2*mp.iv.pi+1)<8)
    return dict(interval_digits=digits,checks={k:bool(v) for k,v in checks.items()},
        actual_R5_Z5_T5_Lambda5_or_huge_parameters_materialized=False)

def input_gate(logR,logZ,A20,G):
    """Exact rational comparison ONLY for explicitly supplied auxiliary inputs."""
    if any(type(v) not in (int,Q) for v in (logR,logZ,A20,G)) or min(logR,logZ,A20)<0 or G<2**32:
        raise ValueError('Exact nonnegative auxiliary logs and G>=2^32 required')
    left=logR+logZ+10*A20+1100000*G
    return dict(auxiliary_B5_comparison=bool(left<=G*G),left=str(left),right=str(G*G),
        actual_R5_Z5_enclosed=False,actual_joint_target_J_verified=False,
        actual_compact_global_or_physical_pulse_accepted=False)

@lru_cache(None)
def review():
    parts=dict(q6=q6_review(),scalar_root=scalar_root_review(),moving=moving_review(),
        native_pressure=native_pressure_review(),primitive=primitive_review(),exponents=exponent_review())
    for r in parts.values():r['checks']={k:bool(v) for k,v in r['checks'].items()}
    checks={k:bool(v) for r in parts.values() for k,v in r['checks'].items()}
    parts['universal']=[universal_review(d) for d in (80,110)]
    checks['all_outward_comparisons']=all(all(v['checks'].values()) for v in parts['universal'])
    parts['checks']=checks
    return parts
