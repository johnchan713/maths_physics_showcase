"""Nested analytic-core second jets; actual global repair remains an input."""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
from pathlib import Path
import importlib.util
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
def module(folder,file):
    p=HERE.parent/folder/file;prior=sys.path[:];sys.path.insert(0,str(p.parent))
    try:
        s=importlib.util.spec_from_file_location('axis_second_'+folder,p)
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=prior

def nested_core(S,rho,Lambda,h):
    """Bounds remain symbolic in the first Picard exponential M1core."""
    S,rho,Lambda,h=map(rat,(S,rho,Lambda,h))
    if S<2 or not 0<rho<=1 or Lambda<=0 or not 0<h<=F(1,100):
        raise ValueError('Retained core norm, strip, radial scale and h required')
    first=module('axis_first_coefficient','axis_bounds.py').core_envelope(S,rho,Lambda)
    delta=first['delta'];d2=delta/4;x=F(9,2)/Lambda;b=1+x
    B3=2*(256*S+9+S/Lambda)*(1+Lambda)**3
    H0=B3*module('second_coefficient_source','source_bounds.py').leading_angular_upgrade(delta)
    # X disks of radius 1/(4Lambda) stay in 4.75/Lambda<5/Lambda.
    # eta disks of radius delta/2 about delta/2 stay in delta<3delta/2.
    H1_factor=6*(1+4*Lambda)**3*(1+2/delta+8/delta**2)
    op=module('second_coefficient_source','source_bounds.py').operator_and_source(x,b,1,1,h)
    zfactor=2*b*op['Q'];betafactor=6*zfactor*zfactor/d2
    return dict(delta=delta,delta2=d2,Xa=4/Lambda,first_complex_radius=5/Lambda,
        second_complex_radius=x,first_jet_X_disk=1/(4*Lambda),second_jet_X_disk=1/(4*Lambda),
        xi_path_upper=b,leading_radial_three=B3,H0=H0,H1_M1core_factor=H1_factor,
        first_core=first,source_q=op['Q'],source_qxi=op['Qxi'],
        source_outer_eta_strip=delta/2,solution_eta_strip=3*d2/2,
        z2_factor=zfactor,beta2_factor=betafactor,
        Tcore='max(1,H0,H1_M1core_factor*M1core); SAME first analytic seed solution',
        M2core='z=2*b*q*Tcore^2;B=6*z^2/delta2;M2core=(z+B)*exp(B)',
        radial_factors=tuple(F(factorial(r))*(4*Lambda)**r for r in range(4)),
        Ceta1_factor=1+2/d2,
        scope='analytic seed only; identified with actual second inner solution on 0<=X<=Xa; no analytic continuation of flat activation or cutoff asserted')

def joined_jets(S,rho,Lambda,h,x,b):
    """Separate true core and collar exponentials; collar depends on actual T."""
    core=nested_core(S,rho,Lambda,h);x,b=map(rat,(x,b))
    if x<F(41,10)/rat(Lambda):raise ValueError('Use x>=4.1/Lambda enclosing the entire retained tc<=.01 collar')
    a=min(F(1),core['Xa'])
    collar=module('second_coefficient_source','source_bounds.py').second_collar(x,b,h,a,core['delta'])
    core_factors=tuple(core['Ceta1_factor']*v for v in core['radial_factors'])
    return dict(core=core,collar=collar,xi_positive_lower=a,overlap=a*a<=core['Xa'],
        core_Ceta1_factors=core_factors,collar_affine=collar['radial_Ceta1_affine'],
        collar_T_powers=collar['T_powers'],
        expression='J2_r<=core_Ceta1_factors[r]*M2core+Tcollar^p[r]*(A_r*M2collar+B_r)',
        domain='pre-cut second inner solution, 0<=X<=Xa*exp(tc/4), real eta; x and b must enclose that actual interval',
        common_smaller_eta_strip=core['delta2']/2)

def cut_rows(core_factors,collar_affine,T,kappa):
    """Coefficients of M2core,M2collar,1 for the SAME second cutoff."""
    m=module('axis_first_coefficient','axis_bounds.py');T=rat(T)
    core_factors=tuple(map(rat,core_factors));collar_affine=tuple(tuple(map(rat,r)) for r in collar_affine)
    if T<1 or len(core_factors)!=4 or min(core_factors)<0 or len(collar_affine)!=4 or any(len(r)!=2 or min(r)<0 for r in collar_affine):
        raise ValueError('Four nonnegative core/affine rows and T>=1 required')
    powers=(0,0,2,4)
    rows=tuple((core_factors[r],T**powers[r]*collar_affine[r][0],T**powers[r]*collar_affine[r][1]) for r in range(4))
    return m.cut_product(rows,kappa)

def symbolic_review():
    """Literal current-order matrix and full regular source parity; all entries."""
    import sympy as s
    xi,e,h=s.symbols('xi e h');X=s.Symbol('X');A=s.Rational(1,2)+h;D=1-A
    L=1-2*h*e*e;d=1-e*e;lam=4*h;b=-1-h;c=-A;kp=-2*A+lam
    FF,UU,v,FX,UX,Fe,Ue=s.symbols('FF UU v FX UX Fe Ue')
    JF=FF+xi*xi*FX;JU=xi*xi*UX
    ZF=(2*b*e*FF-2*e*xi*xi*FX+d*Fe)/L;ZU=(2*c*e*UU-2*e*xi*xi*UX+d*Ue)/L
    Hc=D*e+d*UU;transport=xi*(s.Rational(1,2)-e*UU)/L+xi*v/2
    a0=s.zeros(6);a1=s.zeros(6)
    a0[0,4]=a0[1,5]=1;a0[2,5]=-1;a0[3,0]=4*xi*FF
    a0[4,0]=2*((b+lam)*(2*e*UU-1)/L+v)
    a0[4,1]=2*(2*(A-lam)*e*JF/L+ZF);a0[4,2]=-4*(D+lam)*e*JF/L
    a0[5,0]=-8*e*xi*xi*FF/L
    a0[5,1]=2*((c+lam)*(2*e*UU-1)/L+ZU+2*(A-lam)*e*JU/L)
    a0[5,2]=-4*(D+lam)*e*JU/L;a0[5,3]=4*kp*e/L
    a0[4,4]=a0[5,5]=2*transport
    a1[4,0]=2*Hc/L;a1[4,1]=a1[4,2]=-2*d*JF/L
    a1[5,1]=2*(Hc-d*JU)/L;a1[5,2]=-2*d*JU/L;a1[5,3]=2*d/L
    p0=s.diag(1,1,1,1,-1,-1);kernel=s.diag(*s.symbols('k0:6'));diag=s.diag(0,0,2,0,3,1)
    # P2,theta2,axial2 are holomorphic functions of X from the FULL regular
    # formulas in the parent, not constants whose pressure/viscosity is removed.
    pp,tt,zz=(s.Function(n)(xi*xi,e) for n in ('P2','theta2','axial2'))
    f=s.Matrix([0,0,0,2*xi*pp,2*tt,2*(zz-2*e*xi*xi*pp/L)])
    wrong=f+s.Matrix([0,0,0,0,xi,0])
    return dict(current_lambda=str(lam),current_pressure_exponent=str(kp),
        nonzero_A0=[[i,j] for i in range(6) for j in range(6) if a0[i,j]!=0],
        nonzero_A1=[[i,j] for i in range(6) for j in range(6) if a1[i,j]!=0],
        checks=dict(current_full_matrix_reflects= a0.subs(xi,-xi)==-p0*a0*p0,
            current_full_derivative_matrix_reflects=a1.subs(xi,-xi)==-p0*a1*p0,
            full_second_source_reflects=f.subs(xi,-xi)==-p0*f,
            pressure_current_coupling_is_minus1_plus2h=s.simplify(kp-(-1+2*h))==0,
            current_pressure_and_derivative_couplings_retained=a0[5,3]!=0 and a1[5,3]!=0,
            diagonal_preserves_parity=diag*p0==p0*diag,
            variable_diagonal_kernel_blocks_annihilate=a1*kernel*a1==s.zeros(6),
            differentiated_right_block_annihilates=a1*kernel*s.diff(a1,e)==s.zeros(6),
            odd_swirl_source_control_breaks_even_X_descent=wrong.subs(xi,-xi)!=-p0*wrong))

def complex_source_review():
    """Independent full symbolic source at complex points; diagnostics only."""
    import sympy as s
    import numpy as np
    parent=module('second_coefficient_source','source_bounds.py')
    X,e,fields,sources=parent.exact_source()
    regular_sources=tuple(s.cancel(g) for g in sources)
    # |X|<=1 and |eta|<=1.01. Sum each polynomial's absolute differentiated
    # coefficients, including eta^j at 1.01, for all raw derivatives required.
    def bound(expr,R,K):
        return max(sum(abs(F(str(v)))*F(101,100)**j for (i,j),v in s.Poly(s.diff(expr,X,r,e,k),X,e).terms()) for r in range(R+1) for k in range(K+1))
    H0=max(F(1),*(bound(f,3,3) for f in fields[:2]));H1=max(F(1),*(bound(f,3,2) for f in fields[2:]))
    op=parent.operator_and_source(1,1,H0,H1,F(1,200))
    source=s.lambdify((X,e),regular_sources[:3],'numpy',cse=True);rows=[]
    for xi in (complex(.3,.2),complex(.4,-.3),complex(-.2,.5)):
        for eta in (complex(-.8,.001),complex(0,-.001),complex(.9,.001)):
            P,T,Z=source(xi*xi,eta);L=1-.01*eta*eta
            values=(0,0,0,2*xi*P,2*T,2*(Z-2*eta*xi*xi*P/L))
            checks=dict(complex_full_source_caps=bool(all(abs(v)<=float(F(op['source'][i]['value'])) for i,v in enumerate(values))),
                no_branch_of_sqrt_in_regular_source=abs(xi)**2<=1 and abs(eta)<=1.01)
            rows.append(dict(xi=[xi.real,xi.imag],eta=[eta.real,eta.imag],checks=checks))
    # An exact axis pressure slope from the FULL radial previous viscosity.
    # Cancel removable average quotients BEFORE evaluating X=0. Direct
    # substitution into the unsimplified symbolic primitive gives nan.
    slope=regular_sources[0].subs({X:0,e:0})
    from types import SimpleNamespace
    original=module('background_summation','hierarchy.py');from residual import Poly
    def poly(expr):return Poly({(i,j):float(c) for (i,j),c in s.Poly(expr,X,e).terms()})
    profiles=[SimpleNamespace(F=poly(ff),U=poly(uu),average_U=poly(s.cancel(s.integrate(uu,(X,0,X))/X)),Pi=Poly({}))
        for ff,uu in ((fields[0],fields[1]),(fields[2],fields[3]))]
    axis=original.lower_sources(profiles,0.,0.,.005,1.)['pressure']
    return dict(rows=rows,manufactured_not_selected_PDE=True,not_interval_certification=True,
        H0=str(H0),H1=str(H1),complete_axis_pressure_slope=str(slope),original_axis_pressure_diagnostic=axis,
        checks=dict(all_nine_complex_full_source_controls=all(all(r['checks'].values()) for r in rows),
            complete_axis_pressure_source_not_set_zero=bool(slope.is_Rational and slope!=0),
            all_full_sources_have_no_radial_denominator=all(not s.fraction(g)[1].has(X) for g in regular_sources),
            removable_axis_quotient_matches_original_hierarchy=bool(abs(axis-float(slope))<1e-12)))

def review():
    core=nested_core(2,1,3,F(1,200));join=joined_jets(2,1,3,F(1,200),2,2)
    foundation=module('core_all_order_foundation','foundation.py');norm=256*2+9+F(2,3)
    radial=[foundation.mixed_bound(r,0,norm,1,F(5),F(1,8))*3**r<=core['leading_radial_three'] for r in range(4)]
    monomial=[]
    for n in range(16):
        for r in range(min(3,n)+1):
            val=F(factorial(n),factorial(n-r))*3**r*F(4)**(n-r)/F(9,2)**n
            monomial.append(val<=core['radial_factors'][r])
    cut=cut_rows(join['core_Ceta1_factors'],join['collar_affine'],3,(1,2,3,5))
    omit=cut_rows(join['core_Ceta1_factors'],join['collar_affine'],3,(1,2,3,0))
    checks=dict(leading_reference_raw_radial_three_enclosed=all(radial),
        leading_four_delta_inside_joint_reference_strip=4*core['delta']<F(1,8),
        first_profile_radial_disks_remain_inside_five=F(9,2)+F(1,4)<5,
        first_profile_eta_disks_remain_inside_first_solution=core['delta']/2+core['delta']/2<3*core['delta']/2,
        second_core_Cauchy_disks_remain_inside_nine_halves=4+F(1,4)<F(9,2),
        second_source_and_Picard_strips_nested=2*core['delta2']==core['source_outer_eta_strip'] and core['solution_eta_strip']<core['source_outer_eta_strip'],
        first_raw_second_eta_factor_includes_two_factorial=core['H1_M1core_factor']==6*13**3*(1+2/core['delta']+8/core['delta']**2),
        complex_xi_ray_length_bound=core['xi_path_upper']**2>=core['second_complex_radius'],
        analytic_monomials_obey_second_core_Cauchy=all(monomial),
        core_and_positive_collar_overlap=join['overlap'],
        separate_core_collar_source_exponentials_not_identified=core['M2core']!=join['collar']['M2'],
        scalar_source_envelope_and_nested_exponent_not_evaluated=isinstance(core['Tcore'],str) and isinstance(core['M2core'],str),
        third_cutoff_cost_still_present=cut[3][0]-omit[3][0]==5*join['core_Ceta1_factors'][0],
        all_join_cut_coefficients_nonnegative=all(min(r)>=0 for r in cut))
    sym=symbolic_review();complex_=complex_source_review();checks.update(sym['checks']);checks.update(complex_['checks'])
    return dict(nested_core_fixture=serial(core),inner_join_fixture=serial(join),
        cut_rows_fixture=serial(cut),symbolic=sym,complex_source_diagnostics=complex_,checks=checks)
