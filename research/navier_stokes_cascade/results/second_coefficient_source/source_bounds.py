"""Complete second inner source and conditional third jets; no global closure."""
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
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x
@lru_cache(None)
def module(folder,file):
    p=HERE.parent/folder/file;prior=sys.path[:];sys.path.insert(0,str(p.parent))
    try:
        s=importlib.util.spec_from_file_location('second_source_'+folder,p);m=importlib.util.module_from_spec(s)
        sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=prior

Pair=module('collar_third_coefficient','third_bounds.py').Pair

def radial_velocity_cost(h,lam):
    """Raw eta0,1,2 costs given U/Ubar eta<=3; every X derivative identical."""
    h,lam=map(rat,(h,lam))
    if not 0<h<=F(1,100) or lam not in (0,2*h):raise ValueError('Leading or first lower exponent required')
    e=F(101,100);d=F(2021,1000);l=F(97,100);D=F(1,2)-h
    r0=1/l;r1=4*h*e/l**2;r2=4*h/l**2+32*h*h*e*e/l**3
    n0=2*e+2*(D+lam)*e+d
    n1=2+2*e+2*(D+lam)*(1+e)+2*e+d
    n2=4+2*e+2*(D+lam)*(2+e)+2+4*e+d
    return (r0*n0,r0*n1+r1*n0,r0*n2+2*r1*n1+r2*n0)

def leading_angular_upgrade(delta):
    delta=rat(delta)
    if not 0<delta<=1:raise ValueError('Positive retained angular radius required')
    return sum(F(factorial(k))/delta**k for k in range(4))

def shifted_second(g,X):
    """Both coefficient derivatives included; |exponent|,|shift|<=2, nu=1."""
    aa=3*X;bb=Pair(3);cc=Pair(6);ax=Pair(3);ae=3*X;be=Pair(3);ce=Pair(6)
    return (aa*aa*g+2*aa*bb*g+bb*bb*g+
        (aa*ax+bb*ae+cc*aa+aa*cc)*g+
        (bb*be+cc*bb+bb*cc)*g+(bb*ce+cc*cc)*g)

def operator_and_source(x,b,H0,H1,h):
    """Pairs (value,xi gradient), leading X<=3 eta<=3; first X<=3 eta<=2."""
    x,b,H0,H1,h=map(rat,(x,b,H0,H1,h))
    if min(x,b)<0 or b*b<x or min(H0,H1)<1 or not 0<h<=F(1,100):raise ValueError('Enclosing radii, H0/H1>=1 and retained h required')
    e=F(101,100);d=F(2021,1000);L=F(97,100);A=F(1,2)+h;D=F(1,2)-h
    lam1=2*h;lam2=4*h;be=-1-h;ce=-A;kp=-2*A+lam2
    X=Pair(x,2*b);xi=Pair(b,1);g0=Pair(H0,2*b*H0);g1=Pair(H1,2*b*H1)
    v0=20*g0;v1=12*g1
    radial_axial=shifted_second(v0,X)
    # Every density term is charged with its ORIGINAL lower exponent.
    w=((abs(1-lam1)*v1+X*v1+D*e*v1)/L+2*(2*v1+X*v1)
        +v0*(v1/2+X*v1)+v1*(v0/2+X*v0)
        +g0*(2*abs(lam1-1)*e*v1+2*e*X*v1+d*v1)/L
        +g1*(2*e*v0+2*e*X*v0+d*v0)/L+radial_axial)
    pressure=g1*g1+w/2
    Zg1=(3*X+9)*g1;Z2g1=shifted_second(g1,X)
    theta=v1*(X+1)*g1+g1*Zg1+Z2g1
    axial=X*v1*g1+g1*Zg1+Z2g1
    source=[Pair(0)]*3+[2*xi*pressure,2*theta,2*axial+4*e*X*pressure/L]
    # CURRENT-order matrix lambda2; leading ZF0/ZU0 still use lambda0.
    jf=(1+X)*g0;ju=X*g0;zg=(3*X+9)*g0
    hc=D*e+d*g0;transport=xi*(F(1,2)+e*g0)/L+xi*v0/2
    zero=Pair(0);a0=[[zero for _ in range(6)] for _ in range(6)];a1=[[zero for _ in range(6)] for _ in range(6)]
    a0[0][4]=a0[1][5]=a0[2][5]=Pair(1);a0[3][0]=4*xi*g0
    a0[4][0]=2*(abs(be+lam2)*(2*e*g0+1)/L+v0)
    a0[4][1]=2*(2*abs(A-lam2)*e*jf/L+zg)
    a0[4][2]=4*abs(D+lam2)*e*jf/L;a0[4][4]=2*transport
    a0[5][0]=8*e*X*g0/L
    a0[5][1]=2*(abs(ce+lam2)*(2*e*g0+1)/L+zg+2*abs(A-lam2)*e*ju/L)
    a0[5][2]=4*abs(D+lam2)*e*ju/L;a0[5][3]=Pair(4*abs(kp)*e/L);a0[5][5]=2*transport
    a1[4][0]=2*hc/L;a1[4][1]=a1[4][2]=2*d*jf/L
    a1[5][1]=2*(hc+d*ju)/L;a1[5][2]=2*d*ju/L;a1[5][3]=Pair(2*d/L)
    rows0=[sum(row,zero) for row in a0];rows1=[sum(row,zero) for row in a1];allpairs=rows0+rows1+source
    return dict(Q=max(F(1),*(p.v for p in allpairs)),Qxi=max(p.d for p in allpairs),
        A0=[p.record() for p in rows0],A1=[p.record() for p in rows1],source=[p.record() for p in source],
        full_w1=w.record(),pressure_P2=pressure.record(),radial_axial=radial_axial.record(),
        theta_source=theta.record(),axial_source=axial.record(),
        exponents=dict(current=lam2,lower=lam1,radial_axial_a=F(-1),radial_axial_shift=-1-D,
            theta_axial_lower=be+lam1,axial_lower=ce+lam1,current_theta=be+lam2,current_pressure=kp))

def second_collar(x,b,h,a,delta):
    """T>=max(1,H0,H1), outer strip delta/2; M2 is NOT evaluated at T=1."""
    x,b,h,a,delta=map(rat,(x,b,h,a,delta))
    if a<=0 or a*a>x or not 0<delta<=F(1,1600):raise ValueError('Positive covered lower radius and retained delta<=1/1600 required')
    op=operator_and_source(x,b,1,1,h);d2=delta/4
    jets=module('collar_third_coefficient','third_bounds.py').coefficient_third(op['Q'],op['Qxi'],a,d2)
    zfactor=2*b*op['Q'];betafactor=6*zfactor*zfactor/d2
    return dict(q=op['Q'],qxi=op['Qxi'],delta2=d2,outer_source_strip=delta/2,
        z2_factor=zfactor,beta2_factor=betafactor,
        radial_Ceta1_affine=jets['C_eta1_radial_jets'],T_powers=(0,0,2,4),
        M2='z2=2*b*q*T^2;B2=6*z2^2/delta2;M2=(z2+B2)*exp(B2), on 3delta2/2',
        interpretation='J2_r<=T^power[r]*(A_r*M2+B_r) on xi>=a;not a second-axis or finalized/global H23 bound')

def lower_cut_interface(Lambda,delta,Q1,Q1xi,kappa):
    m=module('axis_first_coefficient','axis_bounds.py');j=m.joined_jets(Lambda,delta,Q1,Q1xi)
    kappa=tuple(map(rat,kappa))
    if len(kappa)!=4 or kappa[0]!=1:raise ValueError('All actual cutoff jets through three, with kappa0=1, required')
    cut=m.cut_product(j['radial_Ceta2_bounds'],kappa)
    raw=tuple(tuple(2*v for v in row) for row in cut)
    return dict(raw_cut_rows=raw,H1_triple=tuple(max(row[i] for row in raw) for i in range(3)),
        meaning='H1<=A*Mcore+B*Mcollar+C on delta/2, including axis;I3 repair is disjoint and does not feed backward into inner averages')

def exact_source(wrong_lower=False):
    """Independent symbolic full source on manufactured profiles, fixed nu=1."""
    import sympy as s
    X,e=s.symbols('X e');h=s.Rational(1,200);D=s.Rational(1,2)-h;A=1-D;L=1-2*h*e*e;d=1-e*e
    F0=1+X/3+e/5+X*e*e/7;U0=s.Rational(1,5)+e/7+X/11+X*e**3/13
    F1=X*(s.Rational(1,11)+e/13+X/17);U1=X*(e/17+X/19+e*e/23)
    ell=4*h if wrong_lower else 2*h
    def Z(k,f):return (2*k*e*f-2*e*X*s.diff(f,X)+d*s.diff(f,e))/L
    def second(k,f):return Z(k-D,Z(k,f))
    def v(U,lam):
        avg=s.integrate(U,(X,0,X))/X
        return (2*e*U-2*(D+lam)*e*avg-d*s.diff(avg,e))/L
    v0=v(U0,0);v1=v(U1,ell)
    w=((1-ell)*v1+X*s.diff(v1,X)+D*e*s.diff(v1,e))/L-2*(2*s.diff(v1,X)+X*s.diff(v1,X,2))
    w+=v0*(v1/2+X*s.diff(v1,X))+v1*(v0/2+X*s.diff(v0,X))
    w+=U0*(2*(ell-1)*e*v1-2*e*X*s.diff(v1,X)+d*s.diff(v1,e))/L
    w+=U1*(-2*e*v0-2*e*X*s.diff(v0,X)+d*s.diff(v0,e))/L-second(-1,v0)
    P=F1*F1-w/2;theta=v1*(X*s.diff(F1,X)+F1)+U1*Z(-1-h+ell,F1)-second(-1-h+ell,F1)
    axial=X*v1*s.diff(U1,X)+U1*Z(-A+ell,U1)-second(-A+ell,U1)
    return X,e,(F0,U0,F1,U1),(P,theta,axial,w)

def symbolic_review():
    import sympy as s
    X,e,h=s.symbols('X e h');D=s.Rational(1,2)-h;L=1-2*h*e*e;d=1-e*e
    def Z(k,f):return (2*k*e*f-2*e*X*s.diff(f,X)+d*s.diff(f,e))/L
    U=e**3;v=(2*e*U-2*D*e*U-d*s.diff(U,e))/L
    control=s.simplify((-Z(-1-D,Z(-1,v))).subs(e,0))
    q,r,u,eta,M=s.symbols('q r u eta M',positive=True)
    B=q+3*u;V1=B*M+q*2*eta*M+q
    V1e=q*eta*M+(B+q*eta)*2*eta*M+q*8*eta**2*M+q*eta
    V2=(r+3*u*u)*M+B*V1+r*2*eta*M+q*V1e+r
    J3=s.expand((1+2*eta)*(V2*u**3+3*V1*u**4+3*M*u**5)/8)
    f=s.Function('f')(X,e);identity=s.simplify(Z(s.Symbol('k'),X*f)-X*Z(s.Symbol('k')-1,f))
    return dict(pointwise_third_angular_control='U0=eta^3,F1=U1=0; lower eta jets vanish at eta0 but w1=6',
        control_w1=str(control),J3_degree_in_q_qxi=s.Poly(J3,q,r).total_degree(),
        checks=dict(regular_axial_identity_preserves_shift_minus_one=identity==0,
            third_leading_angular_derivative_term_is_detected=control==6,
            lower_three_pointwise_angular_traces_do_not_supply_third=all(s.diff(U,e,k).subs(e,0)==0 for k in range(3)) and s.diff(U,e,3)==6,
            third_collar_affine_has_degree_at_most_two=s.Poly(J3,q,r).total_degree()==2,
            third_collar_affine_positive_polynomial=all(c>0 for c in s.Poly(J3,q,r,u,eta,M).coeffs())))

def legacy_review():
    import sympy as s
    import numpy as np
    from types import SimpleNamespace
    original=module('background_summation','hierarchy.py');from residual import Poly
    X,e,fields,sources=exact_source();_,_,_,wrong=exact_source(True)
    def poly(expr):return Poly({(i,j):float(c) for (i,j),c in s.Poly(expr,X,e).terms()})
    profiles=[]
    for FF,UU in ((fields[0],fields[1]),(fields[2],fields[3])):
        profiles.append(SimpleNamespace(F=poly(FF),U=poly(UU),average_U=poly(s.integrate(UU,(X,0,X))/X),Pi=Poly({})))
    cost=lambda z,k:sum(abs(F(str(c)))*(i+1)**3*(j+1)**k for (i,j),c in s.Poly(z,X,e).terms())
    H0=max(F(1),*(cost(z,3) for z in fields[:2]));H1=max(F(1),*(cost(z,2) for z in fields[2:]))
    op=operator_and_source(1,1,H0,H1,F(1,200))
    values=s.lambdify((X,e),[*sources[:3],*[s.diff(z,X) for z in sources[:3]],sources[3],wrong[3]],'numpy',cse=True)
    rows=[]
    for xi,ee in ((.2,-.7),(.2,0),(.2,.8),(.6,-.7),(.6,0),(.6,.8),(1.,-.7),(1.,0),(1.,.8)):
        a0,a1,f=original.matrices(profiles,xi,ee,.005,1.);v=values(xi*xi,ee);P,theta,axial,PX,thetaX,axialX,w,wrongw=v
        expected=np.array([0.,0.,0.,2*xi*P,2*theta,2*(axial-2*ee*xi*xi*P/(1-.01*ee*ee))])
        derivative=np.array([0.,0.,0.,2*P+4*xi*xi*PX,4*xi*thetaX,
            4*xi*axialX-8*ee*xi*(P+xi*xi*PX)/(1-.01*ee*ee)])
        eps=1e-6;plus=original.matrices(profiles,xi+eps,ee,.005,1.);minus=original.matrices(profiles,xi-eps,ee,.005,1.)
        gradients=[(p-m)/(2*eps) for p,m in zip(plus,minus)]
        checks=dict(independent_exact_full_source_matches=bool(np.max(np.abs(f-expected)/np.maximum(1.,np.abs(expected)))<2e-10),
            analytic_source_xi_gradient_matches=bool(np.max(np.abs(gradients[2]-derivative)/np.maximum(1.,np.abs(derivative)))<2e-7))
        for key,mat,grad in (('A0',a0,gradients[0]),('A1',a1,gradients[1])):
            checks[key+'_row_value_caps']=bool(all(sum(abs(mat[i]))<=float(F(op[key][i]['value'])) for i in range(6)))
            checks[key+'_row_gradient_caps']=bool(all(sum(abs(grad[i]))<=float(F(op[key][i]['xi']))*(1+1e-7)+1e-7 for i in range(6)))
        checks['all_source_value_caps']=bool(all(abs(f[i])<=float(F(op['source'][i]['value'])) for i in range(6)))
        checks['all_source_gradient_caps']=bool(all(abs(derivative[i])<=float(F(op['source'][i]['xi'])) for i in range(6)))
        omitted=original.lower_sources(profiles,xi*xi,ee,.005,1.,True)['omega_over_X']
        rows.append(dict(xi=xi,eta=ee,checks=checks,omission_detected=bool(abs(w-omitted)>1e-9),wrong_lower_exponent_detected=bool(abs(w-wrongw)>1e-9)))
    return dict(rows=rows,manufactured_not_selected_PDE=True,finite_difference_tolerance='2e-7,not interval certification',
        checks=dict(original_matrix_and_source_caps_all_pass=all(all(r['checks'].values()) for r in rows),
            dropping_previous_axial_viscosity_detected=any(r['omission_detected'] for r in rows),
            blanket_current_exponent_substitution_detected=any(r['wrong_lower_exponent_detected'] for r in rows)))

def review():
    op=operator_and_source(1,1,1,1,F(1,200));scaled=operator_and_source(1,1,3,7,F(1,200))
    jets=second_collar(1,1,F(1,200),1,F(1,1600));old=module('collar_third_coefficient','third_bounds.py')
    big=old.coefficient_third(op['Q']*49,op['Qxi']*49,1,F(1,6400))['C_eta1_radial_jets']
    lower=lower_cut_interface(3,F(1,1600),2,3,(1,2,3,5))
    # h<=.01 and D+lambda<=.51: this SAME endpoint first-order majorant
    # dominates every numerator and reciprocal cost, including leading eta2.
    uniform=radial_velocity_cost(F(1,100),F(1,50))
    checks=dict(leading_velocity_eta2_cost_below20=max(uniform)<20,
        first_velocity_eta1_cost_below12=max(uniform[:2])<12,
        source_common_T_squared_factor=scaled['Q']<=op['Q']*49 and scaled['Qxi']<=op['Qxi']*49,
        second_collar_T_power_four_keeps_true_M2_input=all(big[r][i]<=jets['radial_Ceta1_affine'][r][i]*7**jets['T_powers'][r] for r in range(4) for i in range(2)),
        second_outer_strip_matches_first_complex_strip=jets['outer_source_strip']==F(1,3200) and jets['delta2']==F(1,6400),
        leading_Cauchy_upgrade_charges_third_eta_jet=leading_angular_upgrade(F(1,4))==421,
        first_cutoff_raw_input_retains_factorial_conversion=all(lower['H1_triple'][i]>=row[i] for row in lower['raw_cut_rows'] for i in range(3)))
    symbolic=symbolic_review();legacy=legacy_review();checks.update(symbolic['checks']);checks.update(legacy['checks'])
    return dict(operator_source_fixture=serial(op),second_collar_fixture=serial(jets),first_cut_interface=serial(lower),
        symbolic=symbolic,original_hierarchy_diagnostics=legacy,checks=checks)
