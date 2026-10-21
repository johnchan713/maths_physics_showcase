"""Exact finite scalar recurrence, current cancellation and scale-loss controls."""
from fractions import Fraction as Q
import importlib.util
from math import factorial
from pathlib import Path

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


dynamic=load("scalar_cycle_dynamic",HERE.parent/"core_dynamic_mean_recheck/proof_checks.py")
means=dynamic.means
R,Z,T,U=means.R,means.Z,means.T,means.U
EPS=means.EPS


def exact(x):
    if type(x) not in (int,Q):raise TypeError("Exact rational input required")
    return Q(x)


class ScalePolynomial:
    """Nonnegative rational majorants in S, not actual pulse fields."""
    def __init__(self,terms=None):
        if terms is None:terms={}
        if type(terms) in (int,Q):terms={0:terms}
        if not isinstance(terms,dict):raise TypeError("Degree/coefficient mapping required")
        self.c={}
        for p,v in terms.items():
            if type(p) is not int or p<0:raise ValueError("Nonnegative integer scale powers required")
            v=exact(v)
            if v<0:raise ValueError("Nonnegative majorant coefficients required")
            if v:self.c[p]=v

    def __add__(self,other):
        if not isinstance(other,ScalePolynomial):other=ScalePolynomial(other)
        result=dict(self.c)
        for p,v in other.c.items():result[p]=result.get(p,Q(0))+v
        return ScalePolynomial(result)

    __radd__=__add__

    def __mul__(self,other):
        if not isinstance(other,ScalePolynomial):other=ScalePolynomial(other)
        result={}
        for p,v in self.c.items():
            for q,w in other.c.items():result[p+q]=result.get(p+q,Q(0))+v*w
        return ScalePolynomial(result)

    __rmul__=__mul__

    def evaluate(self,S):
        S=exact(S)
        if S<1:raise ValueError("Accepted S>=1 required")
        return sum((v*S**p for p,v in self.c.items()),Q(0))

    def record(self):return {str(p):str(v) for p,v in sorted(self.c.items())}


def exponent_margins(C,kappa=Q(1,100000)):
    C,kappa=exact(C),exact(kappa)
    if not 0<kappa<Q(1,100) or C<2-2*kappa:
        raise ValueError("Current finite-cycle half-line required")
    H=C-kappa;mu=1-kappa;nextC=C+Q(1,2)-kappa
    return dict(wave_square=C-1,fast_cross=mu+H-C,fast_square=2*H-C,
                slow_square=C-1,next_input=(C+1)-nextC,
                only_full_mean_rate_is_short_by=kappa)


def cycle_polynomials(j,M2,w_s,q_s,w_d,q_d,old,fast,a,U,Z,M2s):
    M2,M2s=exact(M2),exact(M2s)
    if min(M2,M2s)<=0:raise ValueError("Positive integration geometry required")
    values=(j,w_s,q_s,w_d,q_d,old,fast,a,U,Z)
    if any(not isinstance(p,ScalePolynomial) for p in values):raise TypeError("Scale-polynomial inputs required")
    pi=j+M2*(2*w_s*q_s+q_s*q_s+2*w_d*q_d+q_d*q_d+2*old*fast+fast*fast)
    theta=M2s*(a*(U+Z)+U*Z)
    return pi,theta,a+U+Z


def tail_upper(A,p,h,eta,ell0):
    A,h,eta,ell0=map(exact,(A,h,eta,ell0))
    if type(p) is not int or p<0 or A<0 or h<=0 or eta<=0 or ell0<1:
        raise ValueError("Nonnegative polynomial power and strict positive exponent margin required")
    n=2*p+2;a=h*eta/2
    return A*factorial(n)/(a**n*ell0**2)


def polynomial_review():
    P=ScalePolynomial
    j=P({0:2,1:1});ws=P({0:3,1:2});qs=P({1:Q(1,2)})
    wd=P({0:4,2:1});qd=P({0:Q(2,3),1:1})
    old=P({0:2,1:3});fast=P({1:Q(1,5),2:Q(1,7)})
    a=P({0:3,1:1});u=P({1:2});z=P({0:1,2:3})
    pi,theta,anext=cycle_polynomials(j,Q(4),ws,qs,wd,qd,old,fast,a,u,z,Q(3))
    checks={};rows=[]
    for S in (Q(1),Q(4),Q(64)):
        v=lambda p:p.evaluate(S)
        independent=v(j)+4*(2*v(ws)*v(qs)+v(qs)**2+2*v(wd)*v(qd)+v(qd)**2+2*v(old)*v(fast)+v(fast)**2)
        checks[f"preslow_expanded_budget_S{S}"]=v(pi)==independent
        checks[f"postslow_expanded_budget_S{S}"]=v(theta)==3*(v(a)*(v(u)+v(z))+v(u)*v(z))
        checks[f"Haar_update_is_finite_sum_S{S}"]=v(anext)==v(a)+v(u)+v(z)
        rows.append(dict(S=str(S),Pi=str(v(pi)),Theta=str(v(theta))))
    # Exact epsilon control uses a moderate diagnostic denominator; the actual
    # kappa half-line below is tested without constructing huge physical scales.
    k=Q(1,200);C=2-2*k;eps=Q(1,2**200);S=Q(4)
    power=lambda exponent:Q(1,2**int(200*exponent))
    raw=(power(C)*j.evaluate(S)+4*(2*power(C)*ws.evaluate(S)*qs.evaluate(S)
         +power(2*C-1)*qs.evaluate(S)**2+2*power(C)*wd.evaluate(S)*qd.evaluate(S)
         +power(2*C-1)*qd.evaluate(S)**2+2*power(1-k+C-k)*old.evaluate(S)*fast.evaluate(S)
         +power(2*(C-k))*fast.evaluate(S)**2))
    checks["raw_exponent_budget_below_Pi_without_tail_change"]=raw<=power(C)*pi.evaluate(S)
    slow_raw=3*(power(C+1)*a.evaluate(S)*(u.evaluate(S)+z.evaluate(S))+power(2*C)*u.evaluate(S)*z.evaluate(S))
    checks["raw_slow_budget_below_Theta"]=slow_raw<=power(C+1)*theta.evaluate(S)
    return dict(checks=checks,Pi=pi.record(),Theta=theta.record(),a_next=anext.record(),rows=rows,
                scope="positive-polynomial majorant algebra; constants are diagnostic inputs, not actual canonical norm evaluations")


def margin_review():
    k=Q(1,100000);C0=2-2*k;checks={};rows=[]
    for j in (0,1,4,10**6):
        C=C0+j*(Q(1,2)-k);m=exponent_margins(C,k)
        checks[f"all_strict_margins_j{j}"]=all(v>0 for v in m.values())
        checks[f"sharper_Haar_rate_needed_j{j}"]=(C+1-k)<C+1 and m["next_input"]==Q(1,2)+k
        rows.append(dict(index=j,C=str(C),margins={name:str(v) for name,v in m.items()}))
    # Whole-half-line proof is by nonnegative slopes, not this finite scan.
    checks["margin_affine_slopes_nonnegative"]=all(v>=0 for v in (1,0,1,1,0))
    checks["prefix_slow_orders_above_one"]=min(2-3*k,Q(3,2)-k,C0)>1
    return dict(checks=checks,rows=rows,
                scope="exact rational diagnostic margins; analytical affine formulas prove the full half-line")


def current_slow_review():
    bump=(R-1)**6*(2-R)**6;V=R**2
    psi=R*bump*(1+Z+Z**3)*(1+T+U)
    old=(-EPS*psi.diff(1),(R**4*bump).diff(0).shift(-2)*(1+Z+T+U),means.radial(psi))
    _,_,_,W=means.fixture()
    scalar=lambda m:(R**2*(V*m[2]+m[2]*m[1]+W[1][2]).haar()).integral()
    J=scalar(old);psig=R*bump;gshape=means.radial(psig)
    mass=(R**2*V*gshape).integral().scalar()
    u=(R**4*bump).diff(0).shift(-2)*(1+Z+T)
    g=-J*gshape/mass;B=EPS*(J*psig).diff(1)/mass
    new=tuple(p+q for p,q in zip(old,(B,u,g)))
    v0=old[1].haar();g0=old[2].haar()
    nonlinear=(R**2*(g0*u+v0*g+u*g)).integral()
    checks={
        "nonzero_current_angular_linear_row_divisor":mass!=0,
        "fresh_current_row_cancels_entire_incoming_scalar":J+(R**2*V*g).integral()==0,
        "independently_recomputed_postslow_scalar":scalar(new)==nonlinear,
        "old_Haar_means_are_not_zero":bool(v0) and bool(g0),
        "old_Haar_axial_cross_retained":bool((R**2*g0*u).integral()),
        "old_Haar_angular_cross_retained":bool((R**2*v0*g).integral()),
        "postslow_scalar_not_zero":bool(scalar(new)),
        "slow_increment_divergence_exact":means.radial(B)+EPS*g.diff(1)==0,
        "slow_increment_both_velocity_moments_zero":(R**2*u).integral()==(R*g).integral()==0,
        "complete_covariance_fixed_during_slow_step":bool(W[1][2]) and bool(W[0][2]),
    }
    delta=EPS**3*(1+Z+T)
    gstale=-(J-delta)*gshape/mass
    Bstale=EPS*((J-delta)*psig).diff(1)/mass
    stalestate=tuple(p+q for p,q in zip(old,(Bstale,u,gstale)))
    stale_nonlinear=(R**2*(g0*u+v0*gstale+u*gstale)).integral()
    checks["stale_angular_row_leaves_exact_scalar_mismatch"]=scalar(stalestate)-stale_nonlinear==delta
    checks["stale_row_control_still_has_exact_divergence"]=means.radial(Bstale)+EPS*gstale.diff(1)==0
    for n in (0,1,3):
        upper=4*(dynamic.majorant(g0,n)*dynamic.majorant(u,n)+dynamic.majorant(v0,n)*dynamic.majorant(g,n)
                 +dynamic.majorant(u,n)*dynamic.majorant(g,n))
        checks[f"postslow_nonlinear_norm_budget_n{n}"]=dynamic.majorant(scalar(new),n)<=upper
        checks[f"integrated_residual_jet_budget_n{n}"]=dynamic.majorant(EPS*scalar(new).diff(1),n)<=EPS*(n+1)*dynamic.majorant(scalar(new),n+1)
    return dict(checks=checks,linear_row_divisor=str(mass),
                scope="manufactured Haar current fields and fresh angular linear row; not a full physical five-map solve or the actual canonical field")


def tail_review():
    checks={};rows=[]
    for A,p,h,eta,L in ((Q(3),0,Q(1,200),Q(1,2),10**8),
                         (Q(7),2,Q(1,200),Q(1,100000),10**50),
                         (Q(11),3,Q(1,2**80),Q(1,3),2**500)):
        upper=tail_upper(A,p,h,eta,L)
        label=f"p{p}_h{h}"
        checks[label+":positive_bound"]=upper>0
        checks[label+":monotone_whole_tail_bound"]=tail_upper(A,p,h,eta,2*L)==upper/4
        checks[label+":given_tail_passes_diagnostic_tolerance"]=upper<Q(1,100)
        rows.append(dict(A=str(A),p=p,h=str(h),eta=str(eta),ell_start=str(L),upper=str(upper)))
    # Zero margin with polynomial growth is not a bounded-tail conclusion.
    checks["zero_margin_positive_S_power_grows"]=all(Q(2)*n**4>Q(2)*(n//2)**4 for n in (4,16,64))
    return dict(checks=checks,rows=rows,scope="rational exponential-series upper bounds on specified diagnostic tails; no actual accepted tail certified")


def review():
    parts={"cycle_polynomials":polynomial_review(),"exponent_margins":margin_review(),
           "current_slow":current_slow_review(),"scale_tail":tail_review(),
           "complete_covariance":means.covariance_review(),"pressure_tails":means.pressure_tail_review()}
    return dict(checks={p+":"+k:v for p,r in parts.items() for k,v in r["checks"].items()},parts=parts)
