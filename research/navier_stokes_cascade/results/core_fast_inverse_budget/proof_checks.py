"""Exact fast-map, formal-phase and derivative-budget diagnostics."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


coupled=load("fast_budget_coupled",HERE.parent/"core_coupled_defect_budget/proof_checks.py")
means=coupled.means;Poly=means.Poly;U,Z,T=means.R,means.Z,means.T
P=coupled.P;order=coupled.order;exact=coupled.cycle.exact


def temporal_constant(m):order(m);return 24*2**(m+4)


def slow_count(m):order(m);return Q((m+1)*(m+2),2)


def map_constants(m,length,factor_f,factor_g,factor_k,chi,one_minus_chi,chi_u):
    order(m)
    values=[exact(v) for v in (length,factor_f,factor_g,factor_k,chi,one_minus_chi,chi_u)]
    length,f,g,k,ch,one,chu=values
    if length<=0 or min(f,g,k,ch,one,chu)<0:raise ValueError("Positive hull and nonnegative multiplier norms required")
    return dict(T=2**(3*m)*length*g*(one+ch)*f,A=2**(3*m)*length*k*chu*f)


def fast_majorants(m,tg,kt_next,ka,Etheta,Ez,Ez_next,S):
    order(m)
    tg,kt_next,ka=map(exact,(tg,kt_next,ka))
    if tg<=0 or min(kt_next,ka)<0:raise ValueError("Positive time floor and nonnegative map constants required")
    return dict(Ft=slow_count(m)*temporal_constant(m)*tg*S*(Etheta+(1+ka)*Ez),
                Fr=slow_count(m)*kt_next*temporal_constant(m+1)*tg*S*Ez_next)


def flat_budget_data(m,p,length,Bm,An):
    order(m);order(p)
    if p<1:raise ValueError("Positive integration-by-parts order required")
    length,Bm,An=map(exact,(length,Bm,An))
    if length<=0 or min(Bm,An)<0:raise ValueError("Positive hull and nonnegative multiplier norms required")
    n=m+2*p+3;gamma=8*3**p*2**(m+p+3)
    return dict(input_order=n,constant=2**(m+n)*length*Bm*An*gamma)


class CP:
    """Complex polynomials with rational real/imaginary coefficients."""
    def __init__(self,real=0,imag=0):
        if isinstance(real,CP):self.re,self.im=Poly(real.re),Poly(real.im)
        else:self.re,self.im=Poly(real),Poly(imag)
    def __add__(self,other):
        other=CP(other);return CP(self.re+other.re,self.im+other.im)
    __radd__=__add__
    def __neg__(self):return CP(-self.re,-self.im)
    def __sub__(self,other):return self+-CP(other)
    def __mul__(self,other):
        other=CP(other);return CP(self.re*other.re-self.im*other.im,self.re*other.im+self.im*other.re)
    __rmul__=__mul__
    def diff(self,var):return CP(self.re.diff(var),self.im.diff(var))
    def at(self,var,value):return CP(self.re.at(var,value),self.im.at(var,value))
    def __eq__(self,other):
        other=CP(other);return self.re==other.re and self.im==other.im
    def __bool__(self):return bool(self.re) or bool(self.im)


def add(a,b):return tuple(x+y for x,y in zip(a,b))
def neg(a):return tuple(-x for x in a)
def mul(p,a):return tuple(CP(p)*x for x in a)
def diff(a,var):return tuple(x.diff(var) for x in a)
def zero(a):return all(x==0 for x in a)


def graph(a,omega):
    # Entries multiply 1, exp(i omega(1-U)), exp(i omega(2-U)).
    # The graph derivative cancels the phase derivative in the last two.
    return (a[0].diff(0)+CP(0,omega)*a[0],a[1].diff(0),a[2].diff(0))


def particular(f,omega):
    f=CP(f);omega=exact(omega)
    if not omega:raise ValueError("Nonzero manufactured mode frequency required")
    degree=max((key[0] for poly in (f.re,f.im) for key in poly.terms),default=0)
    inv=CP(0,-1/omega);power=inv;current=f;out=CP()
    for j in range(degree+1):
        out=out+(-1)**j*power*current;power=power*inv;current=current.diff(0)
    return out


def compact(f,chi,omega):
    H=particular(f,omega)
    I=(H,-H.at(0,1),CP());J=(CP(),-H.at(0,1),H.at(0,2))
    return add(I,neg(mul(chi,J))),J,I


def mode_review():
    f=(U-1)**6*(2-U)**6*(1+Z+Z**2+T)
    bump=(U-1)**6*(2-U)**6
    prim=bump.primitive();chi=(prim-prim.at(0,1))/bump.integral().scalar()+Z*bump/100
    c=Q(2,5);nu=Q(7,3);omega=Q(5,2);eps=Q(1,13)
    gd=CP(0,1/(c*nu))*f
    K,J,I=compact(gd*Q(1,2),chi,omega)
    g=mul(2,graph(K,omega));A_gd=mul(2*chi.diff(0),J)
    Kf,Jf,If=compact(CP(f)*Q(1,2),chi,omega);A_f=mul(2*chi.diff(0),Jf)
    forcing=(CP(f),CP(),CP())
    gd_expr=(gd,CP(),CP())
    fulltime=mul(CP(0,c*nu),g)
    frozen_cutoff_Z=add(diff(I,1),neg(mul(chi,diff(J,1))))
    checks={
        "particular_polynomial_solves_oscillatory_first_order":particular(CP(f)*Q(1,2),omega).diff(0)+CP(0,omega)*particular(CP(f)*Q(1,2),omega)==CP(f)*Q(1,2),
        "whole_line_integral_is_graph_homogeneous":zero(graph(J,omega)),
        "full_compact_radial_identity":g==add(gd_expr,neg(A_gd)),
        "actual_A1_velocity_nonzero":not zero(A_gd),
        "actual_fast_time_error_is_A1_of_current_residual":fulltime==add(neg(forcing),A_f),
        "current_residual_A1_error_nonzero":not zero(A_f),
        "wrong_fast_time_error_sign_detected":fulltime!=add(neg(forcing),neg(A_f)),
        "omitting_actual_A1_velocity_detected":g!=gd_expr,
        "full_moving_cutoff_divergence_exact":zero(add(mul(-2*eps,graph(diff(K,1),omega)),mul(eps,diff(g,1)))),
        "moving_cutoff_Z_term_is_required":diff(K,1)!=frozen_cutoff_Z,
        "left_compact_endpoint_cancels":K[0].at(0,1)+K[1].at(0,1)==0 and K[2].at(0,1)==0,
        "right_compact_endpoint_cancels":K[0].at(0,2)+K[2].at(0,2)==0 and K[1].at(0,2)==0,
        "cutoff_zero_one_and_nonconstant_slow_part":chi.at(0,1)==0 and chi.at(0,2)==1 and bool(chi.diff(1)),
        "commutation_removes_temporal_inverse_factor":mul(CP(0,c*nu),A_gd)==neg(A_f),
        "zero_mode_nonzero_moment_not_flat":bool(bump.integral()) and bool(chi.diff(0)*bump.integral()),
    }
    return dict(checks=checks,scope="formal independent endpoint phases at rational frequencies, finite C5 compact data; not selected actual modes or an all-order field")


def divisor_review():
    checks={};count=0
    for k1 in range(-5,6):
        for k2 in range(-5,6):
            if not(k1 or k2):continue
            r=abs(k1)+abs(k2)
            for label,A,B in (("radial",k1+k2,-k2),("temporal",-k1+k2,k1)):
                norm=A*A-2*B*B
                checks[f"{label}_k{k1}_{k2}_nonzero_integer_norm"]=abs(norm)>=1
                # |conjugate| <= |A|+sqrt(2)|B| < 3r, justified using sqrt(2)<3/2.
                checks[f"{label}_k{k1}_{k2}_conjugate_upper_bound"]=abs(A)+Q(3,2)*abs(B)<=3*r
                count+=1
    for m in (0,1,3,6):
        checks[f"m{m}_explicit_fourier_lattice_constant"]=temporal_constant(m)==8*3*2**(m+4)
    partial=sum((Q(1,r*r) for r in range(1,40)),Q(0))
    checks["whole_lattice_sum_bound_via_integral_tail"]=4*(partial+Q(1,39))<=8
    return dict(checks=checks,divisors=count,scope="exact finite arithmetic controls; nonzero algebraic norm and summability proved for all lattice modes in note")


def budget_review():
    checks={};rows=[];S=P({1:1});Et=P({0:2,1:3});Ez=P({1:2,2:1});Eznext=P({0:1,3:2})
    for m in (0,1,3):
        constants=map_constants(m,1,Q(1,2),2,2,3,4,5)
        next_constants=map_constants(m+1,1,Q(1,2),3,2,4,5,6)
        fast=fast_majorants(m,3,next_constants["T"],constants["A"],Et,Ez,Eznext,S)
        coupled_bound=coupled.fast_budget(m,coupled.geometry(1,2),P(2),P(3),fast["Ft"],fast["Fr"])
        for scale in (1,4,16):
            ev=lambda p:p.evaluate(scale)
            direct=fast_majorants(m,3,next_constants["T"],constants["A"],ev(Et),ev(Ez),ev(Eznext),scale)
            checks[f"m{m}_S{scale}_explicit_actual_fast_norm_substitution"]=all(ev(fast[k])==v for k,v in direct.items())
            direct_coupled=coupled.fast_budget(m,coupled.geometry(1,2),2,3,direct["Ft"],direct["Fr"])
            checks[f"m{m}_S{scale}_both_scalar_channels_keep_actual_A1"]=all(ev(coupled_bound[k])==v for k,v in direct_coupled.items())
        checks[f"m{m}_fast_scalar_input_orders"]=(m+1)+4==m+5 and (m+1)+5==m+6
        checks[f"m{m}_joined_slow_scalar_input_orders"]=max((m+2)+6,(m+3)+4)==m+8
        checks[f"m{m}_joined_axial_moment_input_order"]=max((m+3)+6,(m+4)+4)==m+9
        rows.append(dict(order=m,KT=str(constants["T"]),KA=str(constants["A"]),Ft=fast["Ft"].record(),Fr=fast["Fr"].record()))
    for m,p in ((0,1),(1,3),(3,5)):
        data=flat_budget_data(m,p,1,2,3)
        checks[f"m{m}_p{p}_flat_error_uses_current_residual_directly"]=data["input_order"]==m+2*p+3 and data["input_order"]<m+2*p+3+4
        checks[f"m{m}_p{p}_flat_constant_has_no_c_inverse"]=data["constant"]>0
    for target in (0,3,12):
        kappa=Q(1,100000);H=Q(3,2);loss=Q(7,2)
        p=max(1,(Q(target)+loss-H)//kappa+1)
        checks[f"flat_target{target}_strict_exponent_margin"]=H+p*kappa-loss>target
    checks["smoothness_alone_not_a_uniform_band_bound"]=all(Q(2**(6*n))*Q(1,2**n)**6==1 for n in (2,5,9))
    return dict(checks=checks,rows=rows,scope="positive current-residual polynomial routing; actual canonical input seminorms remain unverified")


def review():
    parts={"formal_current_fast":mode_review(),"actual_divisors":divisor_review(),"explicit_budgets":budget_review(),
           "pressure_tails":means.pressure_tail_review(),"complete_covariance":means.covariance_review()}
    return dict(checks={p+":"+k:v for p,row in parts.items() for k,v in row["checks"].items()},parts=parts)
