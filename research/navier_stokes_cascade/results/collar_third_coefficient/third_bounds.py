"""Actual activated leading third jets and a first-coefficient J3 bound.

Exponentials stay symbolic/separated. Rational fixtures are not the selected
enormous core. Matrix envelopes retain all regular pressure/source terms.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
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
@lru_cache(None)
def module(name,file):
    p=HERE.parent/name/file;old=sys.path[:];sys.path.insert(0,str(p.parent))
    try:
        s=importlib.util.spec_from_file_location('third_collar_'+name,p)
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=old

def leading_third(S,rho,Lambda,t,ell):
    S,rho,Lambda,t,ell=map(rat,(S,rho,Lambda,t,ell))
    if S<2 or not 0<rho<=1 or Lambda<=0 or not 0<t<=F(1,100) or not 0<ell<t/20:
        raise ValueError('Retained positive core, strip and actual collar ranges required')
    Xa=4/Lambda;x=F(41,10)/Lambda;b=1+x;delta=rho/(800*S)
    B=2*(256*S+9+S/Lambda)*(1+Lambda)**3
    j1=25*S;j2=1600*S/t+1000*S*S
    j3=30400*S/t**2+128000*S*S/t+50000*S**3
    expjets=(F(1),j1,j2+j1*j1,j3+3*j1*j2+j1**3)
    rlog=(B,x*B,(x+x*x)*B,(x+3*x*x+x**3)*B)
    flog=tuple(sum(F(comb(m,j))*rlog[m-j]*expjets[j] for j in range(m+1)) for m in range(4))
    # U_a,y=(1-e_a) U_r,y, including both derivatives of e_a.
    ulog=(B*(1+ell*x),x*B,(x+x*x+64*x/t)*B,
        (x+3*x*x+x**3+128*(x+x*x)/t+1216*x/t**2)*B)
    c=module('cutoff_higher_jet_interface','jet_bounds.py')
    def radial(logs,m):
        if m==0:return logs[0]
        return sum(abs(v)*logs[j] for j,v in enumerate(c.log_coefficients(m)))/Xa**m
    # A SUM, not a sampled supremum: it bounds all four radial orders of
    # both fields when multiplied by E>=1. It also encloses the core side.
    K=B+sum(radial(flog,m)+radial(ulog,m) for m in range(4))
    q=50*ell*S;power=(q.numerator+q.denominator-1)//q.denominator
    H=K*(1+1/delta+2/delta**2)
    return dict(Xa=Xa,X_upper=x,xi_upper=b,xi_lower=min(F(1),Xa),delta=delta,
        reference_third_jet=B,J_y=j1,J_yy=j2,J_yyy=j3,
        swirl_log_factors=flog,axial_log_bounds=ulog,
        exponential_power=power,radial_third_factor=K,leading_mixed_third_factor=H,
        convention='actual K3<=factor*2^power,H3<=factor*2^power;radial<=3,eta<=2 on2delta')

@dataclass(frozen=True)
class Pair:
    """Nonnegative bounds for a function and its xi derivative; no cancellation."""
    v:F
    d:F=F(0)
    def __post_init__(self):
        object.__setattr__(self,'v',rat(self.v));object.__setattr__(self,'d',rat(self.d))
        if min(self.v,self.d)<0:raise ValueError('Pair bounds must be nonnegative')
    def __add__(self,o):
        o=o if isinstance(o,Pair) else Pair(o);return Pair(self.v+o.v,self.d+o.d)
    __radd__=__add__
    def __mul__(self,o):
        o=o if isinstance(o,Pair) else Pair(o);return Pair(self.v*o.v,self.d*o.v+self.v*o.d)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=rat(o)
        if o<=0:raise ValueError('Only a positive constant divisor is allowed')
        return Pair(self.v/o,self.d/o)
    def __pow__(self,n):
        if type(n) is not int or n<0:raise ValueError('Nonnegative integer pair power required')
        out=Pair(1)
        for _ in range(n):out=out*self
        return out
    def record(self):return dict(value=str(self.v),xi=str(self.d))

def operator_pairs(x,b,H,h):
    """Direct positive-product majorants of the EXACT system.py rows, nu=1.

    H bounds leading X<=3,eta<=2 on2delta, including the radial average.
    x bounds X, b bounds xi. Signed terms are individually charged.
    """
    x,b,H,h=map(rat,(x,b,H,h))
    if min(x,b)<0 or b*b<x or H<1 or not 0<h<=F(1,100):raise ValueError('Valid enclosing radii,H>=1 and retained h required')
    e=F(101,100);d=F(2021,1000);L=F(97,100);A=F(1,2)+h;D=F(1,2)-h;lam=2*h
    be=-1-h;ce=-A;kp=-1
    X=Pair(x,2*b);xi=Pair(b,1);g=Pair(H,2*b*H);v=12*g
    # Any needed v,v_X,v_XX,v_eta has this pair, via the regular average.
    w=(v+X*v+D*e*v)/L+v*(v/2+X*v)+g*(2*e*(v+X*v)+d*v)/L+2*(2*v+X*v)
    aa=3*X;bb=Pair(3);cc=Pair(6);ax=Pair(3);ae=3*X;be_eta=Pair(3);ce_eta=Pair(6)
    Z2=(aa*aa*g+2*aa*bb*g+bb*bb*g+
        (aa*ax+bb*ae+cc*aa+aa*cc)*g+
        (bb*be_eta+cc*bb+bb*cc)*g+(bb*ce_eta+cc*cc)*g)
    jf=(1+X)*g;ju=X*g;zg=(3*X+9)*g
    hc=D*e+d*g;transport=xi*(F(1,2)+e*g)/L+xi*v/2
    zero=Pair(0);a0=[[zero for _ in range(6)] for _ in range(6)]
    a1=[[zero for _ in range(6)] for _ in range(6)]
    a0[0][4]=a0[1][5]=a0[2][5]=Pair(1)
    a0[3][0]=4*xi*g
    a0[4][0]=2*(abs(be+lam)*(2*e*g+1)/L+v)
    a0[4][1]=2*(2*abs(A-lam)*e*jf/L+zg)
    a0[4][2]=4*abs(D+lam)*e*jf/L;a0[4][4]=2*transport
    a0[5][0]=8*e*X*g/L
    a0[5][1]=2*(abs(ce+lam)*(2*e*g+1)/L+zg+2*abs(A-lam)*e*ju/L)
    a0[5][2]=4*abs(D+lam)*e*ju/L;a0[5][3]=Pair(4*abs(kp)*e/L);a0[5][5]=2*transport
    a1[4][0]=2*hc/L;a1[4][1]=a1[4][2]=2*d*jf/L
    a1[5][1]=2*(hc+d*ju)/L;a1[5][2]=2*d*ju/L;a1[5][3]=Pair(2*d/L)
    source=[zero,zero,zero,xi*w,2*Z2,2*e*X*w/L+2*Z2]
    rows0=[sum(row,zero) for row in a0];rows1=[sum(row,zero) for row in a1]
    allpairs=rows0+rows1+source
    return dict(Q=max(F(1),*(p.v for p in allpairs)),Qxi=max(p.d for p in allpairs),
        A0=[p.record() for p in rows0],A1=[p.record() for p in rows1],source=[p.record() for p in source],
        full_w0=w.record(),full_shifted_axial=Z2.record())

def coefficient_third(Q,Qxi,a,delta):
    """Return affine coefficients (A,B) of bounds A*M+B, never exp(M)."""
    Q,Qxi,a,delta=map(rat,(Q,Qxi,a,delta))
    if Q<1 or Qxi<0 or a<=0 or not 0<delta<=1:raise ValueError('Q>=1,Qxi>=0,a>0,0<delta<=1 required')
    def add(*xs):return tuple(sum(v[i] for v in xs) for i in range(2))
    def mul(v,c):return tuple(c*x for x in v)
    M=(F(1),F(0));M1=(2/delta,F(0));M2=(8/delta**2,F(0))
    B=Q+3/a;Bx=Qxi+3/a**2
    wx=add(mul(M,B),mul(M1,Q),(F(0),Q))
    wxe=add(mul(M,Q/delta),mul(M1,B+Q/delta),mul(M2,Q),(F(0),Q/delta))
    wxx=add(mul(M,Bx),mul(wx,B),mul(M1,Qxi),mul(wxe,Q),(F(0),Qxi))
    J=[M,mul(M,1/(2*a)),add(mul(wx,1/(4*a*a)),mul(M,1/(4*a**3))),
       add(mul(wxx,1/(8*a**3)),mul(wx,3/(8*a**4)),mul(M,3/(8*a**5)))]
    ceta=1+2/delta
    return dict(W_xi=wx,W_xi_eta=wxe,W_xixi=wxx,
        scalar_radial_jets=J,C_eta1_radial_jets=[mul(v,ceta) for v in J],
        C_eta1_loss=ceta,strip_Picard_loss=delta/2,
        domain='sqrt(X)>=a>0 on cutoff collar;not an axis bound or whole background/pulse bound')

def symbolic_review():
    import sympy as s
    xi=s.symbols('xi',positive=True);P=s.Function('P')(xi)
    f1=P/(2*xi);f2=s.diff(f1,xi)/(2*xi);f3=s.diff(f2,xi)/(2*xi)
    expected=s.diff(P,xi,2)/(8*xi**3)-3*s.diff(P,xi)/(8*xi**4)+3*P/(8*xi**5)
    y=s.symbols('y');f=s.Function('f')(y);J=s.Function('J')(y)
    third=s.exp(J)*(s.diff(f,y,3)+3*s.diff(f,y,2)*s.diff(J,y)+3*s.diff(f,y)*(s.diff(J,y,2)+s.diff(J,y)**2)+f*(s.diff(J,y,3)+3*s.diff(J,y)*s.diff(J,y,2)+s.diff(J,y)**3))
    e=s.Function('e')(y);u=s.Function('u')(y)
    checks=dict(third_X_from_regular_P_Q=s.simplify(f3-expected)==0,
        missing_inverse_radius_terms_detected=s.simplify(f3-s.diff(P,xi,2)/(8*xi**3))!=0,
        complete_swirl_third_product=s.simplify(s.diff(f*s.exp(J),y,3)-third)==0,
        axial_activation_third=s.expand(s.diff((1-e)*s.diff(u,y),y,2)-((1-e)*s.diff(u,y,3)-2*s.diff(e,y)*s.diff(u,y,2)-s.diff(e,y,2)*s.diff(u,y)))==0,
        Phi_third_majorant_below_S=F(6,20**3)*F(8,5)**4<1,
        core_third_field_factor_below_two=F(6,20**3)*F(8,5)**4<2,
        log_shear_second_bound=F(75700)<100000,
        actual_sigma_second_from_previous_bound=module('cutoff_higher_jet_interface','jet_bounds.py').sigma_bound(2)==1216,
        explicit_Picard_diagnostic_bound=(8+1536<2**11 and 11+3072<4000))
    return dict(checks=checks)

def legacy_review():
    """Independent old system values/finite differences; diagnostics, not proof."""
    import numpy as np
    system=module('first_correction_construction','system.py')
    fixtures=module('physical_residual_budget','review.py');base,_=fixtures.fixtures();rows=[]
    # Polynomial absolute-coefficient majorants on X<=4,|eta|<=1.01.
    H=F(100);op=operator_pairs(4,2,H,F(1,100));step=1e-5
    for xi,eta,h in product((.2,.8,1.2),(-1.,0.,1.),(.005,.01)):
        value=system.matrices(base,xi,eta,h,1.)
        plus=system.matrices(base,xi+step,eta,h,1.);minus=system.matrices(base,xi-step,eta,h,1.)
        grad=[(a-b)/(2*step) for a,b in zip(plus,minus)]
        # h<=.01; use its own exact envelope, not the .01 coefficients.
        op=operator_pairs(4,2,H,F(str(h)))
        def checks(arr,tag):
            if tag=='source':return [np.max(np.abs(arr))]
            return np.sum(np.abs(arr),axis=1)
        ok=[]
        for i,tag in enumerate(('A0','A1','source')):
            bound=op[tag]
            if tag=='source':
                ok.extend([all(abs(v)<=float(F(b['value'])) for v,b in zip(value[i],bound)),
                           all(abs(v)<=float(F(b['xi']))+1e-7 for v,b in zip(grad[i],bound))])
            else:
                ok.extend([all(v<=float(F(b['value'])) for v,b in zip(checks(value[i],tag),bound)),
                           all(v<=float(F(b['xi']))+1e-7 for v,b in zip(checks(grad[i],tag),bound))])
        rows.append(dict(xi=str(xi),eta=str(eta),h=str(h),checks=dict(value_rows_and_full_sources=bool(all(ok[::2])),xi_derivative_diagnostics=bool(all(ok[1::2])))))
    A0,A1,f=system.matrices(base,.8,.4,.005,1.)
    return dict(manufactured_not_selected_core=True,derivative_method='centered difference h=1e-5;not interval certification',
        leading_H='100 bounds the given degree-two polynomial fixture and its X<=3,eta<=2 and average jets on X<=4,|eta|<=1.01',
        rows=rows,checks=dict(all_18_original_system_value_and_gradient_controls=all(all(r['checks'].values()) for r in rows),
            pressure_and_both_shifted_sources_nonzero=bool(min(abs(f[3]),abs(f[4]),abs(f[5]))>1e-5),
            both_pressure_column_couplings_retained=bool(A0[5,3]!=0 and A1[5,3]!=0)))

def review():
    lead=leading_third(2,F(1,8),4,F(1,100),F(1,10000));refined=leading_third(2,F(1,8),4,F(1,1000),F(1,100000))
    # Every positive operator expression has degree <=2 in H; E>=1.
    op=operator_pairs(lead['X_upper'],lead['xi_upper'],lead['leading_mixed_third_factor'],F(1,200))
    jets=coefficient_third(2,3,1,F(1,4));M=2**4000
    J=tuple(A*M+B for A,B in jets['C_eta1_radial_jets'])
    old=module('cutoff_higher_jet_interface','jet_bounds.py')
    k=old.radial_cutoff((1,33,1217,10**9),F(1,64),1);B=old.repair_bump_bounds(3)
    repaired=old.extended_jets(k,J,B,17,9,3)
    sym=symbolic_review();legacy=legacy_review()
    zero_grad=coefficient_third(2,0,1,F(1,4));wide=coefficient_third(2,3,1,F(1,2))
    checks=dict(symbolic_all=all(sym['checks'].values()),legacy_all=all(legacy['checks'].values()),
        true_activation_width_cost_retained=refined['radial_third_factor']>lead['radial_third_factor'],
        no_huge_exponential_materialized=lead['exponential_power']==1,
        operator_derivative_and_pressure_sources_retained=op['Qxi']>0 and F(op['source'][3]['xi'])>0,
        missing_operator_derivative_underestimates_J3=jets['scalar_radial_jets'][3][0]>zero_grad['scalar_radial_jets'][3][0],
        eta_strip_loss_needed=jets['scalar_radial_jets'][3][0]>wide['scalar_radial_jets'][3][0],
        interface_keeps_positive_unknown_third_cutoff_cost=repaired[3]>J[3],
        physical_order3_bound_positive=old.physical_bound(repaired,4,3)>0)
    return dict(actual_leading_factor_fixture=serial(lead),actual_operator_factor_fixture=serial(op),
        separated_operator_power=2*lead['exponential_power'],
        conditional_coefficient_jet_fixture=serial(jets),
        auxiliary_M_upper='2^4000,valid only for Q=2,b=2,Picard loss=1/8',
        conditional_repaired_order3_bound=str(repaired[3]),symbolic=sym,retained_system_diagnostic=legacy,checks=checks)
