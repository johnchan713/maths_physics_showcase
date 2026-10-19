"""Global axial/source envelopes relative to retained regional lemmas."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb,factorial
import importlib.util
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
def rat(x):
    if isinstance(x,(bool,float)):raise TypeError('Exact rational inputs required')
    return F(x)
def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [serial(v) for v in x]
    return x
@lru_cache(None)
def moments():
    p=HERE.parent/'signed_pressure_source_moments/moments.py'
    s=importlib.util.spec_from_file_location('global_axial_pointwise_moments',p)
    m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m

def source_bounds(h,nu,Xv,U3,R1,R2):
    """U3>=sup|U|3; R1>=sup|U_X|2; R2>=sup|U_XX|1, on [0,Xv].

    No source moment cancellation is used. The construction keeps nu=1.
    """
    h,nu,Xv,U3,R1,R2=map(rat,(h,nu,Xv,U3,R1,R2))
    if not 0<h<=F(1,100) or nu<0 or Xv<=0 or min(U3,R1,R2)<0:
        raise ValueError('Need 0<h<=.01, nonnegative viscosity/norms, positive support radius')
    A=F(1,2)+h;D=F(1,2)-h;l=1-2*h;RL=1/l+4*h/l**2
    v2=36*U3
    VX=RL*((4*A+6)*U3+4*Xv*R1)
    VXX=RL*((4*A+10)*R1+4*Xv*R2)
    w=RL*(VX+4*D*v2)+v2*(VX+v2/2)+U3*RL*(4*VX+6*v2)+2*nu*VXX
    CU=64*Xv**2*R2+704*Xv*R1+1664*U3
    return dict(U_C3=U3,U_X_C2=R1,U_XX_C1=R2,v_C2=v2,V_C1=Xv*v2,
        VX_C1=VX,VXX_C1=VXX,w0_C1=w,shifted_U_C1=CU,
        absolute_source_mass=Xv*w,absolute_source_first_mass=Xv**2*w/2)

def global_bounds(h,nu,Xv,Xa,Kaxis,Ahi,log_first,log_second_difference):
    """Log budgets bound |U_y|2 and |U_yy-U_y|1 after the axis-side region.

    The proof supplies H^2 and H^3*N, at the unchanged N. Exact rational
    calls here are auxiliary envelope controls, not actual profile values.
    """
    Xa,Kaxis,Ahi,log_first,log_second_difference=map(rat,(Xa,Kaxis,Ahi,log_first,log_second_difference))
    if Xa<=0 or min(Kaxis,Ahi,log_first,log_second_difference)<0 or rat(Xv)<Xa:
        raise ValueError('Positive ordered radii and nonnegative retained region bounds required')
    U3=max(Kaxis,2*Ahi,F(352))
    R1=max(Kaxis,log_first/Xa);R2=max(Kaxis,log_second_difference/Xa**2)
    return source_bounds(h,nu,Xv,U3,R1,R2)

def scalar_axial_bound(X,U3,R1,R2):
    X,U3,R1,R2=map(rat,(X,U3,R1,R2))
    if min(X,U3,R1,R2)<0:raise ValueError('Nonnegative radius and full mixed norm bounds required')
    return 64*X**2*R2+704*X*R1+1664*U3

def exact_budget():
    # Moderate rational hierarchy control: it checks arithmetic, not actual values.
    H=F(2)**64;Ahi=F(2)**16;N=H**32+1
    mod=4*Ahi*H*N**2;repair=H**20/N;M=H**4*N**2
    first=Ahi+H+1;second=2*Ahi+2*H*N+2
    # Reference R0 value and first three derivatives from its three exact pieces.
    start=(F(1,50),F(1),F(32)*50,F(1216)*50**2)
    tail=(F(11),F(1)+11*32,F(2)*32+11*1216,F(3)*1216+11*220000)
    raw=[max(a,b) for a,b in zip(start,tail)]
    R0=sum(v/factorial(j) for j,v in enumerate(raw))
    beta=sum(F(10)**100/F(3,10)**j/factorial(j) for j in range(4))
    later=3*16*(2*R0+2*beta)
    bump_mass=sum(F(10)**100*10**(j+1)/factorial(j) for j in range(4))
    # Every normalized pullback coefficient is N^(-a-1), at most 1/N.
    rows=[dict(slow=a,angular=k,phase=c,weight=f'N^(-{a+1})')
        for a in range(4) for k in range(4-a) for c in range(4-a-k)]
    return dict(moderate_hierarchy_not_selected_profile=True,
        reference_R0_derivative_caps=serial(raw),reference_R0_J3_upper=str(R0),
        reference_beta_J3_upper=str(beta),reference_U_J3_upper='1e110',
        compact_positive_order='J3 native E,U; no small radial errors asserted',
        pullback_weight_rows=rows,checks=dict(
            weighted_pullback_all_twenty_terms_bounded_by_H_over_N=len(rows)==20 and all(N**(-r['slow']-1)<=1/N for r in rows),
            H_over_N_less_than_one=H/N<1,
            original_and_repair_fit_same_mixed_bound=mod+repair<M,
            final_repair_C3_mixed_is_below_one=repair<1,
            selective_first_log_jet_does_not_need_N=first<H**2,
            selective_second_log_jet_retains_N=second<H**2*N,
            log_second_difference_includes_lower_derivative=first+second<H**3*N,
            fixed_frequency_upper_gives_H70=N<2*H**32 and M<H**70,
            actual_unit_mass_width_point_one_third_jet_cap=bump_mass<4*F(10)**104,
            bump_derivative_constant_absorbed_by_C12=4*F(10)**104<(F(10)**9)**12,
            R0_three_joined_pieces_fit_one_e7=R0<F(10)**7,
            actual_width_point_three_reference_bump_fits_one_e103=beta<F(10)**103,
            full_reference_axial_pulse_fits_one_e110=later<F(10)**110,
            reference_taper_three_derivatives_fit_one_e7=8*(1+32+(1216+32)/2+(220000+3*1216+64)/6)<F(10)**7,
            actual_H_above_two_to_1024_dominates_reference=2**1024>10**110))

def symbolic_review():
    import sympy as s
    X,e,h,nu=s.symbols('X eta h nu',real=True);U=s.Function('U')(X,e);M=s.Function('M')(X,e)
    A=s.Rational(1,2)+h;D=s.Rational(1,2)-h;d=1-e**2;L=1-2*h*e**2
    V=(2*e*X*U-2*D*e*M-d*s.diff(M,e))/L
    def reduce(z):
        return s.simplify(z.subs(s.diff(M,X,e),s.diff(U,e)).subs(s.diff(M,X),U))
    VX=(2*A*e*U+2*e*X*s.diff(U,X)-d*s.diff(U,e))/L
    VXX=(2*(A+1)*e*s.diff(U,X)+2*e*X*s.diff(U,X,2)-d*s.diff(U,X,e))/L
    a,b=-2*e*X/L,d/L;k=s.symbols('k',real=True);c=2*k*e/L;cl=2*(k-D)*e/L
    Z=lambda f,cc:a*s.diff(f,X)+b*s.diff(f,e)+cc*f
    full=(a*a*s.diff(U,X,2)+2*a*b*s.diff(U,X,e)+b*b*s.diff(U,e,2)
        +(a*s.diff(a,X)+b*s.diff(a,e)+a*(c+cl))*s.diff(U,X)
        +(b*s.diff(b,e)+b*(c+cl))*s.diff(U,e)+(b*s.diff(c,e)+cl*c)*U)
    Y,y=s.symbols('Y y');p=Y**3+Y**2*e+Y*e**2+e**3
    raw=s.diff(p.subs(Y,s.log(X)),X,2)
    mu=s.symbols('mu',nonzero=True);vp=-mu/X
    plateau=vp*(-vp/2)
    return dict(checks=dict(
        exact_VX_identity=reduce(s.diff(V,X)-VX)==0,
        exact_VXX_needs_only_U_third_total_jets=s.simplify(s.diff(VX,X)-VXX)==0,
        full_shifted_axial_operator=s.simplify(Z(Z(U,c),cl)-full)==0,
        omitted_axial_coefficient_derivatives_detected=s.simplify(Z(Z(U,c),cl)-(full-b*s.diff(c,e)*U))!=0,
        log_second_retains_lower_derivative=s.simplify(raw-(s.diff(p,Y,2)-s.diff(p,Y)).subs(Y,s.log(X))/X**2)==0,
        reserved_plateau_keeps_nonzero_pressure_source=s.simplify(plateau+mu**2/(2*X**2))==0 and plateau!=0))

def pullback_review():
    import sympy as s
    y,e,p=s.symbols('y eta phase');f=y*y*e*p+p**3+e**3+y*e*e
    rows=[]
    for N in (2,7,31):
      beta=f.subs(p,N*y)/N
      for r in range(4):
       for k in range(4-r):
        direct=s.diff(beta,y,r,e,k)
        expanded=sum(comb(r,c)*s.Rational(N)**(c-1)*s.diff(f,y,r-c,e,k,p,c).subs(p,N*y) for c in range(r+1))
        rows.append(dict(N=N,radial=r,angular=k,exact=s.expand(direct-expanded)==0))
    return dict(manufactured_finite_polynomial_not_actual_loop=True,rows=rows,
        checks=dict(all_thirty_independent_chain_controls=all(r['exact'] for r in rows),
            third_pure_phase_loss_is_N_squared=s.diff((p**3).subs(p,7*y)/7,y,3)==6*7**2))

def manufactured_review():
    m=moments();rows=[]
    for index in (0,1):
      for eta in (F(-1),F(-1,3),F(0),F(2,5),F(1)):
       for viscosity in (F(0),F(1),F(7,3)):
        U,V,v,e,L,d=m.fields(index,eta,F(1,200));A=F(101,200);D=F(99,200)
        vx=(2*A*e*U+2*e*m.X*U.dx()-d*U.de())/L
        vxx=(2*(A+1)*e*U.dx()+2*e*m.X*U.dx().dx()-d*U.dx().de())/L
        direct=(v+m.X*v.dx()+D*e*v.de())/L+v*(v/2+m.X*v.dx())+U*(-2*e*(v+m.X*v.dx())+d*v.de())/L-2*viscosity*(2*v.dx()+m.X*v.dx().dx())
        lowered=(vx+D*e*v.de())/L+v*(vx-v/2)+U*(-2*e*vx+d*v.de())/L-2*viscosity*vxx
        for R in (F(0),F(1,4),F(1)):
          pair=lambda f:f.value(R).pair()
          B=source_bounds(F(1,200),viscosity,1,1000,1000,1000)
          p=pair(direct);rows.append(dict(index=index,eta=str(eta),viscosity=str(viscosity),X=str(R),checks=dict(
              velocity_first_identity=pair(V.dx())==pair(vx),
              velocity_second_identity=pair(V.dx().dx())==pair(vxx),
              exact_pointwise_source_and_first_angular_derivative=pair(direct)==pair(lowered),
              auxiliary_triangle_enclosure=sum(abs(q) for q in p)<=B['w0_C1'])))
    U,V,v,e,L,d=m.fields(0,F(1,3),F(1,200))
    vxx=V.dx().dx().value(F(1,4)).pair()
    return dict(auxiliary_polynomial_fields_not_actual_core=True,case_count=len(rows),rows=rows,
        checks=dict(all_ninety_exact_pointwise_controls=all(all(r['checks'].values()) for r in rows),
            both_axis_and_angular_endpoints_retained=any(r['X']=='0' and r['eta']=='-1' for r in rows) and any(r['X']=='0' and r['eta']=='1' for r in rows),
            viscosity_is_not_a_pointwise_zero_term=any(q!=0 for q in vxx)))

def coefficient_review():
    h=F(1,100);l=1-2*h
    inv2=1/l+6*h/l**2+16*h*h/l**3
    # Raw coefficient C1 bounds are a<=8X,b<=8,c<=16, aX<=8,
    # a_eta<=16X,b_eta<=16,c_eta<=32, with |k|,|k-D|<=2.
    return dict(reciprocal_L_C2_upper=str(inv2),checks=dict(
        reciprocal_L_C2_strictly_below_two=inv2<2,
        UXX_coefficient=8*8==64,
        UX_mixed_and_transport_coefficients=2*8*8*2+(8*8+8*16+8*32)==704,
        U_pure_angular_and_lower_coefficients=8*8*6+(8*16+8*32)*2+(8*32+16*16)==1664))

def review():
    s,e,c,p,f=symbolic_review(),exact_budget(),coefficient_review(),manufactured_review(),pullback_review()
    fixture=global_bounds(F(1,200),1,11,1,7,13,101,303)
    return dict(symbolic=s,mixed_and_reference=e,axial_coefficients=c,manufactured=p,pullback=f,
        auxiliary_global_envelope=serial(fixture),checks=dict(
            symbolic_full_source_order_and_axial_identities=all(s['checks'].values()),
            mixed_pullback_repair_and_reference_envelopes=all(e['checks'].values()),
            independent_fixed_frequency_chain_controls=all(f['checks'].values()),
            differentiated_axial_coefficients=all(c['checks'].values()),
            independent_exact_manufactured_source=all(p['checks'].values())))
