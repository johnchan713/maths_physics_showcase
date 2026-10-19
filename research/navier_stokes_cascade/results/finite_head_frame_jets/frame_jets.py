"""Finite-head Borel radial jets and spectral-frame C2 bounds.

All estimates are conditional on actual same-field coefficient/leading
inputs. A scalar countercontrol is not a PDE coefficient family.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb
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
def norms(v,length):
    v=tuple(map(rat,v))
    if len(v)!=length or min(v)<0:raise ValueError('Complete nonnegative jet list required')
    return v
@lru_cache(None)
def module(name,file):
    p=HERE.parent/name/file;s=importlib.util.spec_from_file_location('finite_head_'+name,p)
    m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m

def borel_majorant(order,head,z):
    """Head norms of E_n,U_n through the fixed order, z=q^h<=1.

    Tail n>=J=max(2,order) uses the ORIGINAL diagonal m<=n budget.
    No X derivative acts on chi(c_n q) at fixed q,eta.
    """
    if type(order) is not int or order<1:raise ValueError('Positive integer derivative order required')
    J=max(2,order);head=norms(head,J-1);z=rat(z)
    if not 0<z<=1:raise ValueError('0<q^h<=1 required')
    tail=(z/2)**J/(1-z/2)
    terms=tuple(C*z**(2*n) for n,C in enumerate(head,1))
    total=sum(terms)+tail
    return dict(first_tail_order=J,head_terms=terms,tail=tail,total=total,
        q2h_normalized=total/(z*z),whole_domain_constant=sum(C for C in head)+F(1,2**(J-1)))

def radius_jets(lo,hi):
    lo,hi=map(rat,(lo,hi))
    if lo<=0 or hi<lo:raise ValueError('Positive ordered normalized radii required')
    R=(hi,1/lo,1/lo**3,3/lo**5)
    inv=(1/lo,1/lo**3,3/lo**5,15/lo**7)
    return R,inv
def first_coefficient_norm(S,lo,hi):
    """E1=sqrt(2X)F1; U1 unchanged; input common F1,U1 jets<=3."""
    S=norms(S,4);R,_=radius_jets(lo,hi)
    E=tuple(sum(F(comb(m,j))*R[j]*S[m-j] for j in range(m+1)) for m in range(4))
    return dict(E=E,U=S,combined=max(F(1),*E,*S))

def profile_error_jets(B,lo,hi):
    """B*q^(2h) bounds EACH E,U X jet<=3 on this annulus.

    Return F=e/R and g=(2X F_X,R U_X) error factors through order two.
    Euclidean g norms are bounded by the component sum.
    """
    B=rat(B);lo,hi=map(rat,(lo,hi));R,inv=radius_jets(lo,hi)
    if B<0:raise ValueError('Nonnegative error factor required')
    f=tuple(B*sum(F(comb(m,j))*inv[j] for j in range(m+1)) for m in range(4))
    gt=tuple(hi**2*f[m+1]+2*m*f[m] for m in range(3))
    gz=tuple(B*sum(F(comb(m,j))*R[j] for j in range(m+1)) for m in range(3))
    return dict(F=f,g_theta=gt,g_z=gz,g=tuple(x+y for x,y in zip(gt,gz)))

def spectral_bounds(floor,fjets,gjets):
    """Current native frame N=g/|g|,K=rot N,c=-sqrt(D/P).

    f<=F<=Fmax,-gtheta>=t,|g|<=G,D>=d; fjets/ gjets are
    total first/second radial derivative bounds. Floors must already hold.
    """
    f,Fmax,t,G,d=map(rat,floor);F1,F2=norms(fjets,2);g1,g2=norms(gjets,2)
    if min(f,Fmax,t,G,d)<=0 or Fmax<f or G<t:raise ValueError('Positive consistent current spectral floors required')
    Pmin=2*f*t;Pmax=2*Fmax*G;Dmax=G*G+2*Fmax*G
    if d>Dmax:raise ValueError('Discriminant lower bound exceeds its enclosure')
    cmin=min(F(1),d/Pmax);cmax=max(F(1),Dmax/Pmin)
    N1=g1/t;N2=g2/t+3*g1*g1/(t*t)
    D1=2*G*g1+2*F1*G+2*Fmax*g1
    D2=2*g1*g1+2*G*g2+2*F2*G+4*F1*g1+2*Fmax*g2
    P1=2*(F1*G+Fmax*g1);P2=2*(F2*G+2*F1*g1+Fmax*g2)
    q1=D1/Pmin+Dmax*P1/Pmin**2
    q2=D2/Pmin+2*D1*P1/Pmin**2+Dmax*P2/Pmin**2+2*Dmax*P1**2/Pmin**3
    c1=q1/(2*cmin);c2=q2/(2*cmin)+q1*q1/(4*cmin**3)
    lam1=2*(F1*cmax+Fmax*N1*cmax+Fmax*c1)
    lam2=2*(F2*cmax+Fmax*N2*cmax+Fmax*c2+2*F1*N1*cmax+2*F1*c1+2*Fmax*N1*c1)
    return dict(Pmin=Pmin,Pmax=Pmax,Dmax=Dmax,D1=D1,D2=D2,P1=P1,P2=P2,
        squared_c_derivatives=(q1,q2),cmin=cmin,cmax=cmax,
        N=(N1,N2),K=(N1,N2),c=(c1,c2),growth=(lam1,lam2))

def radial_physical(jets,rhi):
    """Coefficients AFTER extracting q^-1/2,q^-1; not loss-free jets."""
    j1,j2=norms(jets,2);rhi=rat(rhi)
    if rhi<=0:raise ValueError('Positive normalized radius required')
    return (rhi*j1,j1+rhi*rhi*j2)

def second_coefficient_control(n,z=F(1,16)):
    if type(n) is not int or n<1:raise ValueError('Positive integer oscillation frequency required')
    z=rat(z)
    if not 0<z<=F(1,2):raise ValueError('Profile diagonal n=2 budget needs z^2<=1/4')
    # psi is the SAME smooth plateau from the previous analytical control;
    # its derivative supports are disjoint, and psi=1 near X=1/2.
    p1,p2,p3=F(256),F(77824),F(112640000);a=1/(1+2*p1+p2)
    C2=F(1);C3=max(C2,a*(n+3*p1+3*p2/n+p3/F(n*n)))
    actual_center_third=a*n*z**4
    corrected=borel_majorant(3,(0,C3),z)
    naive=z*z/2
    return dict(n=n,z=z,low_order_stage_norm=C2,third_order_coefficient_norm=C3,
        value_first_second_uniform_bound=C2,third_perturbation_at_center=actual_center_third,
        naive_first_coefficient_plus_half_bound=naive,corrected_bound=corrected['total'],
        checks=dict(diagonal_order2_profile_budget=C2*z*z<=F(1,4),
            first_two_derivative_norms_below_one=a*(1+2*p1+p2)<=1,
            corrected_finite_head_bound_encloses_center=actual_center_third<=corrected['total']),
        naive_bound_fails=actual_center_third>naive)

def symbolic_review():
    import sympy as s
    X,z=s.symbols('X z',positive=True);e=s.Function('e')(X);u=s.Function('u')(X)
    R=s.sqrt(2*X);f=e/R
    checks=dict(tail3_geometric=s.simplify((z/2)**3/(1-z/2)-z**3/(8*(1-z/2)))==0,
        third_order_normalized_finite_head=s.simplify((z**2*s.Symbol('C1')+z**4*s.Symbol('C2')+(z/2)**3/(1-z/2))/z**2-(s.Symbol('C1')+z**2*s.Symbol('C2')+z/(8*(1-z/2))))==0)
    for m in range(3):
        checks['exact_shear_product_'+str(m)]=s.simplify(s.diff(2*X*s.diff(f,X),X,m)-(2*X*s.diff(f,X,m+1)+2*m*s.diff(f,X,m)))==0
    for m,v in enumerate((1/R,-R**-3,3*R**-5,-15*R**-7)):
        checks['inverse_radius_jet_'+str(m)]=s.simplify(s.diff(1/R,X,m)-v)==0
    D=s.Function('D')(X);P=s.Function('P')(X);q=D/P
    q2=s.diff(D,X,2)/P-2*s.diff(D,X)*s.diff(P,X)/P**2-D*s.diff(P,X,2)/P**2+2*D*s.diff(P,X)**2/P**3
    checks['squared_c_second_quotient']=s.simplify(s.diff(q,X,2)-q2)==0
    c=s.sqrt(D);checks['square_root_second']=s.simplify(s.diff(c,X,2)-(s.diff(D,X,2)/(2*c)-s.diff(D,X)**2/(4*c**3)))==0
    # Independent two-component normal identity; g has nonconstant magnitude.
    g=s.Matrix([X+1,X*X+2]);G=s.sqrt(g.dot(g));N=g/G;proj=s.eye(2)-N*N.T
    rhs=proj*g.diff(X,2)/G-2*N.diff(X)*s.diff(G,X)/G-N*(N.diff(X).dot(N.diff(X)))
    checks['normal_second_projected_identity']=all(s.simplify(v)==0 for v in N.diff(X,2)-rhs)
    r,qscale=s.symbols('r qscale',positive=True);v=s.Function('v')(X)
    checks['physical_second_radial_scale']=s.simplify(s.diff(v.subs(X,r*r/(2*qscale)),r,2)-(s.diff(v,X)+2*X*s.diff(v,X,2)).subs(X,r*r/(2*qscale))/qscale)==0
    return dict(checks=checks)

def spectral_diagnostics():
    """Exact symbolic test profiles; numerical evaluation is diagnostic only."""
    import sympy as s
    X=s.symbols('X',real=True);f=1+X/10;g=s.Matrix([-4-X/5,1+X/7]);G=s.sqrt(g.dot(g))
    D=g.dot(g)+2*f*g[0];P=-2*f*g[0];c=s.sqrt(D/P);N=g/G;lam=2*f*(-N[0])*c
    caps=spectral_bounds((1,2,4,6,4),(F(1,10),0),(F(1,2),0))
    rows=[]
    for x in (F(0),F(1,4),F(1,2),F(3,4),F(1)):
        sub={X:s.Rational(x.numerator,x.denominator)}
        n1=s.sqrt(N.diff(X).dot(N.diff(X)));n2=s.sqrt(N.diff(X,2).dot(N.diff(X,2)))
        actual={'N':(n1,n2),'c':(s.diff(c,X),s.diff(c,X,2)),'growth':(s.diff(lam,X),s.diff(lam,X,2))}
        checks={name:all(abs(float(v.subs(sub)))<=float(bound) for v,bound in zip(values,caps[name])) for name,values in actual.items()}
        checks['profile_floor_hypotheses']=bool(f.subs(sub)>=1 and f.subs(sub)<=2 and -g[0].subs(sub)>=4 and G.subs(sub)<=6 and D.subs(sub)>=4)
        rows.append(dict(X=str(x),checks=checks))
    return dict(manufactured_not_actual_profile=True,rows=rows,
        checks=dict(all_five_spectral_jet_controls=all(all(v['checks'].values()) for v in rows)))

def review():
    inherited=module('cutoff_higher_jet_interface','jet_bounds.py')
    # Inputs here are auxiliary finite budgets. Actual C2,3 remains an input.
    first=first_coefficient_norm((2,3,5,7),1,4)
    borel=borel_majorant(3,(first['combined'],11),F(1,4));err=profile_error_jets(borel['q2h_normalized'],1,4)
    caps=spectral_bounds((1,2,4,6,4),(3,5),(7,11))
    controls=[second_coefficient_control(n) for n in (1,1000,10**12,10**20)]
    sym=symbolic_review();diag=spectral_diagnostics()
    low=borel_majorant(1,(first['combined'],),1)
    checks=dict(symbolic_all=all(sym['checks'].values()),spectral_diagnostics=all(diag['checks'].values()),
        order3_tail_starts_at_three_and_second_is_in_head=borel['first_tail_order']==3 and len(borel['head_terms'])==2,
        old_value_half_tail_recovered=low['q2h_normalized']==first['combined']+F(1,2),
        all_abstract_second_coefficient_controls_meet_diagonal_profile_budget=all(all(v['checks'].values()) for v in controls),
        first_coefficient_alone_does_not_bound_third_jet=controls[-1]['naive_bound_fails'],
        inverse_radius_and_lower_shear_term_charged=err['F'][3]>err['F'][0] and err['g_theta'][2]>16*err['F'][3],
        strict_current_frame_denominators_retained=caps['cmin']>0 and caps['Pmin']>0,
        physical_factors_match_independent_general_chain=radial_physical(caps['N'],4)==tuple(inherited.physical_bound((0,*caps['N']),4,m) for m in (1,2)),
        physical_scale_losses_not_erased=radial_physical(caps['N'],4)[1]>caps['N'][1])
    return dict(first_coefficient_example=serial(first),finite_head_borel_example=serial(borel),
        native_profile_error_example=serial(err),current_spectral_jet_example=serial(caps),
        physical_N_derivative_factors=serial(radial_physical(caps['N'],4)),
        second_coefficient_profile_controls=serial(controls),symbolic=sym,spectral_diagnostics=diag,
        checks=checks)
