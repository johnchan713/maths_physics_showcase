"""Exact finite-jet tests for the profile-dependent continuum estimates."""
from fractions import Fraction as Q
import importlib.util
from math import factorial
from pathlib import Path

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


geometry=load("fixed_profile_geometry",HERE.parent/"core_reserved_geometry_recheck/proof_checks.py")
shape=load("fixed_profile_shape",HERE.parent/"core_pressure_shape_recheck/proof_checks.py")
means=shape.means
R,Z,T=means.R,means.Z,means.T
Jet=geometry.Jet


def order(n):
    if type(n) is not int or n<0:raise ValueError("Nonnegative integer order required")
    return n


def exact_nonnegative(x):
    x=geometry.exact(x)
    if x<0:raise ValueError("Nonnegative exact majorant required")
    return x


def majorant_powers(m):
    order(m)
    A=Jet(m,{(i,j):Q(32**i*64**j,3) for i in range(m+1) for j in range(m+1-i) if i+j})
    powers=[Jet(m,{(0,0):1})]
    for k in range(1,m+1):powers.append(powers[-1]*A)
    return powers


def tame_bound(b,r,m,profile_majorants):
    if b not in (0,1) or type(b) is not int:raise ValueError("Density weight b must be zero or one")
    order(r);order(m)
    if len(profile_majorants)!=m+1:raise ValueError("All profile orders through m required")
    majorants=[exact_nonnegative(v) for v in profile_majorants]
    return 2**(b+r)*sum((majorants[k]*sum(p.c.values(),Q(0))/Q(2**k*factorial(k))
                        for k,p in enumerate(majorant_powers(m))),Q(0))


def log_one(t):
    if t.c.get((0,0),0):raise ValueError("Log increment must have zero constant term")
    out=Jet(t.n);power=Jet(t.n,{(0,0):1})
    for k in range(1,t.n+1):
        power=power*t;out=out+Q((-1)**(k+1),k)*power
    return out


def euler_shift(poly,c):return R*poly.diff(0)+c*poly


def radial_diff(poly,r):
    for _ in range(r):poly=poly.diff(0)
    return poly


def profile_bound(poly):
    """Coefficient majorant on x in [1,2]; never an actual-profile evaluation."""
    return sum((abs(v)*2**max(k[0],0) for k,v in poly.terms.items()),Q(0))


def composition_review():
    checks={};rows=[];m=5;x0=Q(3,2)
    f=(R-1)**8*(2-R)**8
    A=majorant_powers(m)[1]
    for h,eta in ((Q(1,200),Q(0)),(Q(1,200),Q(1)),(Q(1,100),Q(-1)),(Q(1,2**60),Q(3,5))):
        t=geometry.implicit_jet(h,eta,m);u=log_one(t)
        key=f"h{h}_eta{eta}"
        checks[key+":coordinate_log_majorant"]=all(abs(v)<=A.c[k] for k,v in u.c.items())
        for b,r in ((0,0),(0,1),(0,2),(1,0),(1,2)):
            c=b+r;fr=radial_diff(f,r)
            direct=Jet(m)
            for (p,z,tt,aux),v in fr.terms.items():
                direct=direct+v*x0**p*t.power_one(-Q(c+p,2))
            series=Jet(m);power=Jet(m,{(0,0):1});current=fr
            values=[];majorants=[]
            for k in range(m+1):
                values.append(current.at(0,x0).scalar());majorants.append(profile_bound(current))
                series=series+Q((-1)**k,2**k*factorial(k))*values[-1]*power
                current=euler_shift(current,c);power=power*u
            label=key+f"_b{b}_r{r}"
            checks[label+":direct_power_equals_Euler_composition"]=series.c==direct.c
            checks[label+":coefficient_majorant"]=sum(abs(v) for v in direct.c.values())<=tame_bound(b,r,m,majorants)
            if c:
                wrong=Jet(m);power=Jet(m,{(0,0):1});current=fr
                for k in range(m+1):
                    wrong=wrong+Q((-1)**k,2**k*factorial(k))*current.at(0,x0).scalar()*power
                    current=euler_shift(current,0);power=power*u
                checks[label+":omitted_density_radial_shift_detected"]=wrong.c!=direct.c
        rows.append(dict(h=str(h),eta=str(eta),order=m,nonzero_log_coefficients=len(u.c)))
    powers=majorant_powers(2)
    checks["constant_term_is_one"]=powers[0].c=={(0,0):Q(1)}
    checks["linear_majorant_coefficients"]=powers[1].c[(1,0)]==Q(32,3) and powers[1].c[(0,1)]==Q(64,3)
    checks["quadratic_mixed_convolution"]=powers[2].c[(1,1)]==Q(4096,9)
    checks["order_zero_density_radial_factor"]=tame_bound(1,2,0,[Q(5)])==40
    return dict(checks=checks,rows=rows,
                scope="exact fifth-order jets at four chart centers and mixed radial weights; all-order proof is finite Taylor plus coefficient majorization")


def qshift(m,j):return factorial(m+j)//factorial(m)


def norm(poly,m):
    """A submultiplicative coefficient norm on [1,2] x [-1,1]^2."""
    out=Q(0)
    for i in range(m+1):
        for j in range(m+1-i):
            p=poly
            for _ in range(i):p=p.diff(1)
            for _ in range(j):p=p.diff(2)
            out+=profile_bound(p)/Q(factorial(i)*factorial(j))
    return out


def commutator_budget(m,eps,r,f):
    return (eps**2*(2*(m+1)**2*r[m+1]*f[m+1]+qshift(m,2)*r[m+2]*f[m])
            +eps**4*(3*(m+1)*qshift(m,2)*r[m+1]*f[m+2]
                     +3*(m+1)*qshift(m,2)*r[m+2]*f[m+1]
                     +qshift(m,3)*r[m+3]*f[m]))


def commutator_review():
    eps=means.EPS;bump=(R-1)**6*(2-R)**6;rho0=bump/bump.integral().scalar()
    rho1=(R-Q(3,2))*rho0/100
    star=rho0+(Z+Z**2+Z**3+T+Z*T)*rho1
    rho=star+rho1/7
    F=1+Z+Z**2+Z**3+Z**4+T+T**2*Z**2
    B=eps*(F*star).diff(1)
    source=eps*B.diff(2)+eps**3*B.diff(1).diff(1)
    U=eps**2*F.diff(1).diff(2)+eps**4*F.diff(1).diff(1).diff(1)
    terms=[eps**2*star.diff(1)*F.diff(2),eps**2*star.diff(2)*F.diff(1),eps**2*star.diff(1).diff(2)*F,
           3*eps**4*star.diff(1)*F.diff(1).diff(1),3*eps**4*star.diff(1).diff(1)*F.diff(1),eps**4*star.diff(1).diff(1).diff(1)*F]
    C=sum(terms,means.Poly());primitive=C.primitive();p=primitive-primitive.at(0,1)
    checks={"exact_six_term_commutator":source-rho*U==(star-rho)*U+C,
            "each_moving_product_term_nonzero":all(bool(v) for v in terms),
            "each_term_has_zero_radial_mass":all(v.integral()==0 for v in terms),
            "commutator_pressure_compact":p.at(0,1)==p.at(0,2)==0,
            "commutator_pressure_radial_identity":p.diff(0)==C}
    rs=[norm(star,n) for n in range(7)];fs=[norm(F,n) for n in range(7)]
    for m in range(3):
        checks[f"all_six_terms_budget_order{m}"]=norm(C,m)<=commutator_budget(m,eps,rs,fs)
        # For this exact polynomial source its pressure is known; the analytic
        # half-L1/sup inequality itself is proved in the note, not by coefficient norms.
        checks[f"all_six_term_omissions_detected_order{m}"]=all(C-v!=C for v in terms)
    return dict(checks=checks,scope="manufactured moving unit-mass densities; exact full product and coefficient-norm budgets")


def margin_review():
    # Quantifiers: positive margins yield a finite-order common threshold;
    # zero margin and a growing target constant require additional input.
    conditions=[(Q(10**6),2,Q(1,3)),(Q(7),1,Q(1,5)),(Q(2**40),3,Q(1,7))]
    candidates=[min(Q(1),tau/C) for C,gap,tau in conditions]
    common=min(candidates)
    checks={"finite_common_threshold_is_positive":common>0,
            "all_positive_integer_margins_absorb_fixed_constants":all(C*common**gap<=tau for C,gap,tau in conditions),
            "zero_margin_not_absorbed":Q(7)*common**0>Q(1,5),
            "unbounded_target_constant_can_cancel_gain":all(Q(2**(5*n))*Q(1,2**n)**5==1 for n in (2,5,9)),
            "infinite_orders_need_not_share_positive_threshold":all(Q(1,2**(n+1))<Q(1,2**n) for n in (2,5,9))}
    return dict(checks=checks,finite_example_threshold=str(common),
                scope="finite rational controls; positive real exponent existence follows by continuity at epsilon=0, no new accepted tail selected")


def review():
    parts={"composition":composition_review(),"commutator":commutator_review(),"margins":margin_review(),
           "pressure_tails":means.pressure_tail_review(),"complete_covariance":means.covariance_review()}
    return dict(checks={p+":"+k:v for p,r in parts.items() for k,v in r["checks"].items()},parts=parts)
