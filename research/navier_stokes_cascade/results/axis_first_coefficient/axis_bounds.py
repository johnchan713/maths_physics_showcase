"""Complex-core first-coefficient jets and a real activated-collar join."""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial,comb
from pathlib import Path
import importlib.util
import sys

HERE=Path(__file__).resolve().parent
def rat(x):
    if isinstance(x,(float,bool)):raise TypeError('Exact rational inputs required')
    return F(x)
def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x
@lru_cache(None)
def module(folder,file):
    p=HERE.parent/folder/file;prior=sys.path[:];sys.path.insert(0,str(p.parent))
    try:
        s=importlib.util.spec_from_file_location('axis_first_'+folder,p);m=importlib.util.module_from_spec(s)
        sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=prior

def core_envelope(S,rho,Lambda):
    """Analytic SEED on |Lambda X|<=5, not the activated continuation there."""
    S,rho,Lambda=map(rat,(S,rho,Lambda))
    if S<2 or not 0<rho<=1 or Lambda<=0:raise ValueError('Retained core norm, radius and scale required')
    delta=rho/(800*S);x=5/Lambda;b=1+x
    B=2*(256*S+9+S/Lambda)*(1+Lambda)**2
    H=B*(1+1/delta+2/delta**2)
    W=512*(1+x)*H*H;Z=128*(1+x)**2*H
    A0=128*(1+x+b)*H;A1=32*(1+x)*H;source=(b+3*x)*W+2*Z
    C=max(F(1),A0,A1,source);z=2*C*b;beta=6*z*z/delta
    return dict(delta=delta,Xa=4/Lambda,X_complex_radius=x,xi_path_upper=b,
        reference_radial_two=B,leading_mixed_two=H,full_w0=W,full_shifted_axial=Z,
        A0=A0,A1=A1,source=source,C=C,z=z,beta_upper=beta,
        Mcore='(z+beta_upper)*exp(beta_upper); beta_upper=3*z^2/(delta/2)',
        X_Cauchy_radius=1/(2*Lambda),
        radial_factors=tuple(F(factorial(r))*(2*Lambda)**r for r in range(4)),
        angular_Ceta2_factor=1+2/delta+4/delta**2,
        definition='each F1,U1 X jet<=3 bounded by radial_factor*Mcore on unchanged real core;factorial angular C2 costs angular_Ceta2_factor')

def joined_jets(Lambda,delta,Q,Qxi):
    """Triples (A,B,C): A*Mcore+B*Mcollar+C, after angular Ceta2 loss."""
    Lambda,delta,Q,Qxi=map(rat,(Lambda,delta,Q,Qxi))
    if Lambda<=0 or not 0<delta<=1 or Q<1 or Qxi<0:raise ValueError('Positive scale/strip and valid collar envelopes required')
    Xa=4/Lambda;a=min(F(1),Xa)
    collar=module('collar_third_coefficient','third_bounds.py').coefficient_third(Q,Qxi,a,delta)
    factor=1+2/delta+4/delta**2
    core=tuple(F(factorial(r))*(2*Lambda)**r for r in range(4))
    rows=tuple((factor*core[r],factor*collar['scalar_radial_jets'][r][0],factor*collar['scalar_radial_jets'][r][1]) for r in range(4))
    return dict(Xa=Xa,xi_positive_lower=a,overlap=a*a<=Xa,Ceta2_factor=factor,
        radial_Ceta2_bounds=rows,core_radial_factors=core,
        common_mixed_bound=tuple(max(row[i] for row in rows) for i in range(3)),
        common_raw_mixed_bound=tuple(2*max(row[i] for row in rows) for i in range(3)),
        collar_scalar_affine=collar['scalar_radial_jets'],
        domain='actual inner first coefficient on 0<=X<=Xa exp(tc/4);not yet multiplied by coefficient cutoff or repaired on I3')

def cut_product(rows,kappa):
    """Same eta-independent cutoff; triples remain coefficients of two M bounds."""
    rows=tuple(tuple(map(rat,row)) for row in rows);kappa=tuple(map(rat,kappa))
    if len(rows)!=4 or any(len(row)!=3 or min(row)<0 for row in rows) or len(kappa)!=4 or min(kappa)<0:
        raise ValueError('Four complete nonnegative triple rows and cutoff jets required')
    return tuple(tuple(sum(F(comb(r,j))*kappa[j]*rows[r-j][i] for j in range(r+1)) for i in range(3)) for r in range(4))

def symbolic_review():
    """Every matrix parity entry; variable-kernel derivative-block cancellation."""
    import sympy as s
    xi,e=s.symbols('xi e');P=s.diag(1,1,1,1,-1,-1);D=s.diag(0,0,2,0,3,1)
    # Each independent abstract coefficient is even in xi (a function of X).
    A0=s.zeros(6);A1=s.zeros(6)
    for i,j in ((0,4),(1,5),(2,5),(4,0),(4,1),(4,2),(5,0),(5,1),(5,2),(5,3)):
        A0[i,j]=s.Symbol(f'b{i}{j}')
    for i,j in ((3,0),(4,4),(5,5)):A0[i,j]=xi*s.Symbol(f'b{i}{j}')
    for i,j in ((4,0),(4,1),(4,2),(5,1),(5,2),(5,3)):A1[i,j]=(1+e)*s.Symbol(f'c{i}{j}')
    source=s.Matrix([0,0,0,-xi*s.Symbol('w0'),-2*s.Symbol('CF'),s.Symbol('pU')-2*s.Symbol('CU')])
    kernel=s.diag(*s.symbols('k0:6'));different=s.diag(*s.symbols('l0:6'))
    odd_control=source+s.Matrix([0,0,0,0,xi,0])
    kernel_checks=[]
    for exponent in (0,1,2,3):
        for n in range(7):
            inverse=xi**(n+1)/s.Integer(n+1+exponent)
            kernel_checks.append(s.expand(s.diff(inverse,xi)+exponent*inverse/xi-xi**n)==0)
    return dict(nonzero_A0=[[i,j] for i in range(6) for j in range(6) if A0[i,j]!=0],
        nonzero_A1=[[i,j] for i in range(6) for j in range(6) if A1[i,j]!=0],
        parity=[1,1,1,1,-1,-1],diagonal=[0,0,2,0,3,1],
        checks=dict(all_A0_entries_have_reflection_covariance=A0.subs(xi,-xi)==-P*A0*P,
            all_A1_entries_have_reflection_covariance=A1.subs(xi,-xi)==-P*A1*P,
            full_pressure_and_axial_source_has_reflection_covariance=source.subs(xi,-xi)==-P*source,
            diagonal_commutes_with_reflection=D*P==P*D,
            variable_diagonal_kernel_derivative_blocks_annihilate=A1*kernel*A1==s.zeros(6),
            derivative_hitting_right_block_also_annihilates=A1*kernel*s.diff(A1,e)==s.zeros(6),
            different_radii_kernels_do_not_spoil_annihilation=A1*kernel*s.diff(A1,e)*different==s.zeros(6),
            complex_ray_regular_Volterra_monomial_identities=all(kernel_checks),
            odd_swirl_source_breaks_needed_parity=odd_control.subs(xi,-xi)!=-P*odd_control))

def recurrence(terms=12,drop_pressure=False,odd_source=False):
    """Exact analytic manufactured six-system coefficients, not actual PDE data."""
    import sympy as s
    if type(terms) is not int or terms<1:raise ValueError('Positive finite recurrence length required')
    e=s.Symbol('e');D=(0,0,2,0,3,1);a0=[s.zeros(6) for _ in range(3)];a1=[s.zeros(6) for _ in range(3)]
    for i,j in ((0,4),(1,5),(2,5),(4,0),(4,1),(4,2),(5,0),(5,1),(5,2),(5,3)):
        a0[0][i,j]=s.Rational(i+j+1,19)*(1+e)
        a0[2][i,j]=s.Rational(i+j+1,37)*(2-e)
    for i,j,value in ((0,4,1),(1,5,1),(2,5,-1)):
        a0[0][i,j]=value;a0[2][i,j]=0
    for i,j in ((3,0),(4,4),(5,5)):a0[1][i,j]=s.Rational(i+j+1,23)*(1+e*e)
    for i,j in ((4,0),(4,1),(4,2),(5,1),(5,2),(5,3)):
        a1[0][i,j]=s.Rational(i+j+1,29)*(2-e)
        a1[2][i,j]=s.Rational(i+j+1,41)*(1+e)
    f=[s.Matrix([0,0,0,0,1+e,2-e]),s.Matrix([0,0,0,0 if drop_pressure else -1-e*e,1 if odd_source else 0,0]),s.Matrix([0,0,0,0,1+e*e,3-e])]
    out=[s.zeros(6,1)]
    for n in range(terms):
        rhs=f[n].copy() if n<len(f) else s.zeros(6,1)
        for p in range(min(n,2)+1):rhs+=a0[p]*out[n-p]+a1[p]*out[n-p].diff(e)
        out.append(s.Matrix([s.expand(rhs[i]/(n+1+D[i])) for i in range(6)]))
    return out,a0,a1,f,e

@lru_cache(None)
def recurrence_review():
    import sympy as s
    out,a0,a1,f,e=recurrence();no_pressure=recurrence(drop_pressure=True)[0];odd=recurrence(odd_source=True)[0]
    signs=(1,1,1,1,-1,-1);residual=[]
    for n in range(12):
        rhs=f[n].copy() if n<len(f) else s.zeros(6,1)
        for p in range(min(n,2)+1):rhs+=a0[p]*out[n-p]+a1[p]*out[n-p].diff(e)
        residual.extend(s.expand((n+1+(0,0,2,0,3,1)[i])*out[n+1][i]-rhs[i])==0 for i in range(6))
    return dict(terms=12,exact_scalar_recurrences=len(residual),
        first_axis_X_slopes=[str(out[2][i]) for i in (0,1,3)],
        manufactured_analytic_system_not_selected_PDE=True,
        checks=dict(all_six_component_recurrences=all(residual),
            all_component_parities=all(out[n][i]==0 or (-1)**n==signs[i] for n in range(13) for i in range(6)),
            zero_axis_data=out[0]==s.zeros(6,1),
            exact_regular_average_component=all(s.expand(out[n][2]+s.Rational(n,n+2)*out[n][1])==0 for n in range(13)),
            first_four_components_descend_to_holomorphic_X=all(out[n][i]==0 for n in range(1,13,2) for i in range(4)),
            last_two_are_odd_not_declared_holomorphic_X=out[1][4]!=0 and out[1][5]!=0,
            removing_pressure_source_changes_axial_coefficient=any(out[n][1]!=no_pressure[n][1] for n in range(13)),
            parity_breaking_source_creates_odd_regular_profile=any(odd[n][0]!=0 for n in range(1,13,2))))

def envelope_review():
    b=core_envelope(2,1,3);j=joined_jets(3,F(1,4),2,3)
    foundation=module('core_all_order_foundation','foundation.py')
    B=b['reference_radial_two'];norm=256*2+9+F(2,3)
    checks=[foundation.mixed_bound(r,0,norm,1,F(5),F(1,8))*3**r<=B for r in range(3)]
    monomials=[]
    for degree in range(13):
        for r in range(min(3,degree)+1):
            value=F(factorial(degree),factorial(degree-r))*3**r*F(4)**(degree-r)/5**degree
            monomials.append(value<=F(factorial(r))*6**r)
    cut=cut_product(j['radial_Ceta2_bounds'],(1,2,3,5));without=cut_product(j['radial_Ceta2_bounds'],(1,2,3,0))
    return dict(core=serial(b),joined=serial(j),cut_product_fixture=serial(cut),
        checks=dict(core_majorant_bounds_all_needed_radial_orders=all(checks),
            actual_four_delta_strip_inside_core_joint_majorant=4*b['delta']<F(1,8),
            Cauchy_disks_stay_strictly_inside_complex_core=4+F(1,2)<5,
            positive_collar_and_core_cover_without_axis_division=j['overlap'],
            separate_core_and_collar_M_constants_both_charged=all(a>0 and c>0 for a,c,d in j['radial_Ceta2_bounds']),
            all_join_constants_nonnegative=all(min(row)>=0 for row in j['radial_Ceta2_bounds']),
            second_angular_factor_retained=j['Ceta2_factor']>1+2/F(1,4),
            full_pressure_source_envelope_dominates_both_shifted_viscosities=b['source']>=2*b['full_shifted_axial'],
            complex_path_length_bound_valid=b['xi_path_upper']**2>=b['X_complex_radius'],
            analytic_monomials_obey_core_Cauchy_factors=all(monomials),
            third_cutoff_jet_not_discharged_by_axis_join=cut[3][0]-without[3][0]==5*j['radial_Ceta2_bounds'][0][0],
            no_selected_exponentials_materialized=isinstance(b['Mcore'],str)))

def legacy_review():
    """Compare stored structural patterns to every nonzero original matrix entry."""
    import numpy as np
    from types import SimpleNamespace
    original=module('first_correction_construction','system.py')
    from residual import Poly
    profile=SimpleNamespace(F=Poly({(0,0):1.,(1,0):.2,(0,1):.1}),
        U=Poly({(0,0):.3,(1,0):.4,(0,1):.2}),
        average_U=Poly({(0,0):.3,(1,0):.2,(0,1):.2}),Pi=Poly({}))
    expected=symbolic_review();rows=[]
    for xi,e in ((.3,-.6),(.7,.4),(1.1,.8)):
        A0,A1,f=original.matrices(profile,xi,e,.005,1.)
        patterns=[[[i,j] for i in range(6) for j in range(6) if a[i,j]!=0] for a in (A0,A1)]
        rows.append(dict(xi=xi,eta=e,checks=dict(A0_pattern_is_original=patterns[0]==expected['nonzero_A0'],
            A1_pattern_is_original=patterns[1]==expected['nonzero_A1'],
            original_pressure_and_two_axial_source_entries_nonzero=bool(np.all(f[3:]!=0)))))
    return dict(rows=rows,floating_manufactured_structure_only=True,
        checks=dict(original_system_structure_matches=all(all(r['checks'].values()) for r in rows)))

def review():
    a,b,c,d=symbolic_review(),recurrence_review(),envelope_review(),legacy_review()
    return dict(symbolic=a,recurrence=b,envelopes=c,original_structure=d,
        checks=dict(**a['checks'],**b['checks'],**c['checks'],**d['checks']))
