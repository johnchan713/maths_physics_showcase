"""Exact piecewise potential moments and high target-jet response controls."""
from fractions import Fraction as Q
import importlib.util
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("slow_jet_mean_polynomials",HERE.parent/"complete_covariance_initial_means/means.py")
means=importlib.util.module_from_spec(spec);spec.loader.exec_module(means)


def logs(one=0,g=0,b=0):return {k:Q(v) for k,v in dict(one=one,g=g,b=b).items() if v}
def ladd(a,b):return {k:a.get(k,Q(0))+b.get(k,Q(0)) for k in set(a)|set(b) if a.get(k,Q(0))+b.get(k,Q(0))}
def lscale(a,c):return {k:v*c for k,v in a.items() if v*c}


def integrate(f,a,b):
    p=f.primitive()
    return p.at(0,b)-p.at(0,a)


def integrate_laurent(f,a,b,la,lb):
    out={}
    for (r,z,t,u),c in f.terms.items():
        if z or t or u:raise ValueError("Radial polynomial required")
        if r==-1:out=ladd(out,lscale(ladd(lb,lscale(la,-1)),c))
        else:out=ladd(out,logs(one=c*(b**(r+1)-a**(r+1))/Q(r+1)))
    return out


def integrate_xlog(f,a,b,la,lb):
    out={}
    for (r,z,t,u),c in f.terms.items():
        if z or t or u or r<-1:raise ValueError("Supported radial power required")
        n=r+2
        right=ladd(lscale(lb,b**n/Q(n)),logs(one=-b**n/Q(n*n)))
        left=ladd(lscale(la,a**n/Q(n)),logs(one=-a**n/Q(n*n)))
        out=ladd(out,lscale(ladd(right,lscale(left,-1)),c))
    return out


def radial_fixture():
    R=means.R;g,b=Q(9,4),Q(25,16)
    eta=(R-1)**4*(b-R)**4
    phi1=means.Poly({k:v/g**(k[0]+1) for k,v in eta.terms.items()})
    mu1=integrate(R*eta,Q(1),b).scalar()
    mu3=integrate(R**3*eta,Q(1),b).scalar()
    # lambda=1/4: x^(1-2lambda)=sqrt(x), endpoints have rational square roots.
    def sqrt_moment(poly,lower_root,upper_root):
        return sum((v*(upper_root**(2*r+3)-lower_root**(2*r+3))/Q(r+Q(3,2))
                    for (r,z,t,u),v in poly.terms.items()),Q(0))
    muh=sqrt_moment(eta,Q(1),Q(5,4))
    q=Q(2,3);a=Q(7,5)
    c1=1/(a*muh*g*(1-q))
    sources=(-g*eta,means.Poly(),phi1)
    bounds=((Q(1),b),(b,g),(g,g*b))
    ends=((logs(),logs(b=1)),(logs(b=1),logs(g=1)),(logs(g=1),logs(g=1,b=1)))
    pieces=[];history=means.Poly()
    for f,(left,right),(ll,lr) in zip(sources,bounds,ends):
        raw=(R*f).primitive();A=raw-raw.at(0,left)+history
        pieces.append((left,right,ll,lr,c1*f,c1*A.shift(-1)))
        history=A.at(0,right)
    moment={};weighted={};logrow={}
    for left,right,ll,lr,f,psi in pieces:
        moment=ladd(moment,integrate_laurent(psi,left,right,ll,lr))
        weighted=ladd(weighted,integrate_laurent(R**2*psi,left,right,ll,lr))
        logrow=ladd(logrow,integrate_xlog(f,left,right,ll,lr))
    Hlog=mu1/(a*muh*(1-q))
    H2=mu3*(g*g-1)/(2*a*muh*(1-q))
    checks=dict(
        axial_zero_flux=history==0,
        fractional_moment_dilation=sqrt_moment(phi1,Q(3,2),Q(15,8))==g*q*muh,
        exact_current_angular_target=a*c1*(sqrt_moment(phi1,Q(3,2),Q(15,8))-g*muh)==-1,
        logarithmic_potential_moment=moment==logs(g=-Hlog),
        weighted_potential_moment=weighted==logs(one=-H2),
        independent_log_source_row=logrow==logs(g=Hlog),
        no_spurious_log_endpoint_or_constant="b" not in moment and "one" not in moment,
        gap_potential_is_nonzero=bool(pieces[1][5]),
        both_potential_endpoints_zero=pieces[0][5].at(0,1)==pieces[-1][5].at(0,g*b)==0,
        positive_log_coefficient=Hlog>0 and H2>0,
    )
    gap=integrate_laurent(pieces[1][5],b,g,logs(b=1),logs(g=1))
    checks["omitting_gap_changes_moment"]=bool(gap)
    return dict(checks=checks,pieces=pieces,Hlog=Hlog,H2=H2,
                mu1=mu1,mu3=mu3,muh=muh,g=g,q=q,a=a,
                controls=dict(potential_moment={k:str(v) for k,v in moment.items()},
                              weighted_moment={k:str(v) for k,v in weighted.items()},
                              gap_moment={k:str(v) for k,v in gap.items()},
                              mu1=str(mu1),mu3=str(mu3),mu_half=str(muh)))


def operator_review(radial):
    R,Z,T,eps=means.R,means.Z,means.T,means.EPS
    _,_,_,W=means.fixture();totalP={};totalQ={};checks={}
    for degree in (3,4):
        J=Z**degree/(6 if degree==3 else 24)
        P,Qg,Jz={},{},{}
        for left,right,ll,lr,f,psi in radial["pieces"]:
            gamma=f*J;beta=-eps*psi*J.diff(1)
            # Nonzero old mean curl and radial base, kept in both evaluations.
            bump=(R-left)**4*(right-R)**4
            oldpsi=R*bump*(1+Z+T)
            old=(-eps*oldpsi.diff(1),(R**4*bump).diff(0).shift(-2)*(1+Z),means.radial(oldpsi))
            base=(eps*(1+R+Z),R**2*(2+Z),means.Poly())
            plus=tuple(v+d for v,d in zip(old,(beta,means.Poly(),gamma)))
            minus=tuple(v-d for v,d in zip(old,(beta,means.Poly(),gamma)))
            gp=means.averaged_gr(base,plus,W[0][0],W[1][1],W[0][2])
            gm=means.averaged_gr(base,minus,W[0][0],W[1][1],W[0][2])
            linear=(gp-gm)/2
            at=linear.at(1,0) if degree==3 else linear.diff(1).at(1,0)
            zp=(2*base[2]*plus[2]+plus[2]**2+W[2][2]).haar()
            zm=(2*base[2]*minus[2]+minus[2]**2+W[2][2]).haar()
            full_axial=R*(zp-zm)/2-R**2*linear/2
            axial_at=full_axial.at(1,0) if degree==3 else full_axial.diff(1).at(1,0)
            # Highest target jet isolates axial viscosity; all retained lower terms vanish.
            checks[f"full_operator_principal_piece_{degree}_{left}"]=at==-eps**4*psi
            P=ladd(P,integrate_laurent(at,left,right,ll,lr))
            Qg=ladd(Qg,integrate_laurent(R**2*at,left,right,ll,lr))
            Jz=ladd(Jz,integrate_laurent(axial_at,left,right,ll,lr))
            no_axial=linear-eps**3*beta.diff(1).diff(1)
            omitted=no_axial.at(1,0) if degree==3 else no_axial.diff(1).at(1,0)
            checks[f"axial_viscosity_omission_detected_{degree}_{left}"]=omitted==0 and bool(at)
        checks[f"pressure_principal_coefficient_degree{degree}"]=P==logs(g=eps**4*radial["Hlog"])
        checks[f"axial_defect_principal_coefficient_degree{degree}"]=lscale(Qg,-Q(1,2))==logs(one=-eps**4*radial["H2"]/2)
        checks[f"complete_axial_moment_odd_part_degree{degree}"]=Jz==logs(one=-eps**4*radial["H2"]/2)
        if degree==3:totalP,totalQ=P,Qg
    # Physical normalization of both moments, exact exponents without evaluating s^h.
    for h,lam in ((Q(1,200),Q(1,4)),(Q(1,2**80),Q(1,1000))):
        checks[f"Q_log_moment_exponent_h{h}"]=Q(1,2)-(lam-h)-(1-2*lam)/2==h
        checks[f"Q_weighted_moment_exponent_h{h}"]=Q(3,2)-(lam-h)-(1-2*lam)/2==1+h
    return dict(checks=checks,principal_P={k:str(v) for k,v in totalP.items()},
                principal_Q={k:str(v) for k,v in totalQ.items()},
                scope="exact odd part of full quadratic radial operator; polynomial target-jet isolators, not prescribed incoming states")


def pressure_review(radial):
    K,K2=radial["Hlog"],radial["H2"]
    # c_hat*log(g) is the quantity in the exact cancellation gate.
    tuned_log_product=K2/(2*K)
    exact_gate=radial["mu3"]*(radial["g"]**2-1)/(4*radial["mu1"])
    Z=means.Z;F=Z**3/6;c=2+3*Z+Z**2
    lhs=c*F.diff(1).diff(1).diff(1)-(c*F).diff(1).diff(1).diff(1)
    rhs=-3*c.diff(1)*F.diff(1).diff(1)-3*c.diff(1).diff(1)*F.diff(1)-c.diff(1).diff(1).diff(1)*F
    return dict(checks=dict(
        exact_fourth_order_cancellation_gate=tuned_log_product==exact_gate,
        algebraically_tuned_principal_coefficient_zero=K*tuned_log_product-K2/2==0,
        pressure_third_coefficient_not_cancelled=K>0,
        full_moving_pressure_commutator=lhs==rhs,
        moving_pressure_commutator_survives=lhs.diff(1).at(1,0)==-9,
    ),c_hat_times_log_g_required=str(tuned_log_product),
       scope="algebraic gate only; actual compensator moment has not been shown to satisfy it")


def rate_review():
    rows=[]
    for n in (2,16,1024):
        rows.append(dict(frequency=n,normalized_third_response=n**3,normalized_fourth_response=n**4))
    checks={"same_norm_third_response_unbounded_control":rows[-1]["normalized_third_response"]>10**9,
            "frequency_ratios_increase":all(a["normalized_third_response"]<b["normalized_third_response"] for a,b in zip(rows,rows[1:]))}
    for q in (Q(1,2),Q(15,16),1-Q(1,2**80)):
        scaled=1/(1-q)
        checks[f"inverse_tilt_denominator_q{q}"]=(1-q)*scaled==1
    return dict(checks=checks,rows=rows,scope="finite exact controls; unboundedness argument is the smooth oscillatory family in README.md")


def review():
    r=radial_fixture()
    parts={"radial_moments":dict(checks=r["checks"],controls=r["controls"]),
           "full_operator":operator_review(r),"pressure_gate":pressure_review(r),
           "rate_control":rate_review(),"pressure_tails":means.pressure_tail_review(),
           "complete_covariance":means.covariance_review()}
    return dict(checks={p+":"+k:value for p,part in parts.items() for k,value in part["checks"].items()},parts=parts)
