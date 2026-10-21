"""Exact pointwise pressure-shape diagnostics; no actual compensator replacement."""
from fractions import Fraction as Q
import importlib.util
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("pressure_shape_potential_moments",HERE.parent/"core_slow_jet_loss_recheck/proof_checks.py")
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
means=old.means


def expr(one=0,g=0,b=0):return {(i,j):Q(v) for (i,j),v in (((0,0),one),((1,0),g),((0,1),b)) if v}
def eadd(a,b):return {k:a.get(k,Q(0))+b.get(k,Q(0)) for k in set(a)|set(b) if a.get(k,Q(0))+b.get(k,Q(0))}
def escale(a,c):return {k:v*c for k,v in a.items() if v*c}
def gshift(a,n):return {(i+n,j):v for (i,j),v in a.items()}
def adapt(a):return eadd(expr(one=a.get("one",0)),eadd(expr(g=a.get("g",0)),expr(b=a.get("b",0))))
def display(a):return {f"logg^{i}*logb^{j}":str(v) for (i,j),v in sorted(a.items())}


def fixture():
    R=means.R;r=old.radial_fixture();mass={};second={}
    for left,right,ll,lr,f,psi in r["pieces"]:
        numerator=-psi/r["Hlog"]  # rho_* = numerator/log(g).
        mass=eadd(mass,adapt(old.integrate_laurent(numerator,left,right,ll,lr)))
        second=eadd(second,adapt(old.integrate_laurent(R**2*numerator,left,right,ll,lr)))
    mass=gshift(mass,-1);second=escale(gshift(second,-1),Q(1,2))
    lo,hi=Q(7,4),Q(2)
    w=(R-lo)**8*(hi-R)**8
    h=w.diff(0).diff(0).diff(0)
    M=sum((abs(v)*max(lo**k[0],hi**k[0]) for k,v in h.terms.items()),Q(0))
    delta=1/(2*r["g"]*(r["g"]-1)*max(Q(1),M))
    D=-delta*w.diff(0).diff(0)
    raw=(R*D).primitive();target=(raw-raw.at(0,lo)).shift(-1)
    checks=dict(
        reference_mass_one=mass==expr(one=1),
        reference_second_moment_gate=second=={(-1,0):r["H2"]/(2*r["Hlog"])},
        reference_gap_density_exact=(-r["pieces"][1][5]/r["Hlog"])==R.shift(-2),
        perturbation_support_inside_gap=r["pieces"][1][0]<lo<hi<r["pieces"][1][1],
        perturbation_zero_mass=old.integrate(h,lo,hi)==0,
        perturbation_zero_second_moment=old.integrate(R**2*h,lo,hi)==0,
        exact_positivity_budget=0<delta*M<=1/(2*r["g"]*(r["g"]-1)),
        compact_pressure_primitive=D.at(0,lo)==D.at(0,hi)==0,
        pressure_primitive_derivative=D.diff(0)==-delta*h,
        integrated_fourth_coefficient_zero=old.integrate(R*D,lo,hi)==0,
        pointwise_fourth_coefficient_nonzero=D.at(0,(lo+hi)/2).scalar()>0,
        compact_signed_target_primitive=target.at(0,lo)==target.at(0,hi)==0,
        signed_target_derivative=means.radial(target)==D,
        signed_target_pointwise_nonzero=bool(target.at(0,(lo+hi)/2)),
    )
    return dict(checks=checks,radial=r,lo=lo,hi=hi,w=w,h=h,delta=delta,D=D,target=target,
                controls=dict(reference_mass=display(mass),reference_second_moment=display(second),
                              positivity_delta=str(delta),perturbation_majorant=str(M),
                              pressure_midpoint=str(D.at(0,(lo+hi)/2).scalar()),
                              signed_target_midpoint=str(target.at(0,(lo+hi)/2).scalar())))


def operator_review(f):
    R,Z,T,eps=means.R,means.Z,means.T,means.EPS
    r=f["radial"];_,_,_,W=means.fixture()
    checks={};pressure_history=means.Poly();radial_integral={};radial_log_integral=means.Poly();axial_moment=means.Poly()
    # Split the gap at the perturbation's actual support; no polynomial is used off support.
    first,gap,last=r["pieces"]
    pieces=[first,
            (gap[0],f["lo"],gap[2],{},gap[4],gap[5]),
            (f["lo"],f["hi"],{},{},gap[4],gap[5]),
            (f["hi"],gap[1],{},gap[3],gap[4],gap[5]),last]
    for index,(left,right,ll,lr,gamma_profile,psi) in enumerate(pieces):
        J=Z**4/24;gamma=gamma_profile*J;B=-eps*psi*J.diff(1)
        bump=(R-left)**4*(right-R)**4;oldpsi=R*bump*(1+Z+T)
        m=(-eps*oldpsi.diff(1),(R**4*bump).diff(0).shift(-2)*(1+Z),means.radial(oldpsi))
        plus=tuple(a+b for a,b in zip(m,(B,means.Poly(),gamma)))
        minus=tuple(a-b for a,b in zip(m,(B,means.Poly(),gamma)))
        base=(eps*(1+R+Z),R**2*(2+Z),means.Poly())
        odd=(means.averaged_gr(base,plus,W[0][0],W[1][1],W[0][2])-means.averaged_gr(base,minus,W[0][0],W[1][1],W[0][2]))/2
        source_jet=odd.diff(1).at(1,0)
        numerator=-psi/r["Hlog"]
        perturb=f["delta"]*f["h"] if index==2 else means.Poly()
        # P_v principal coefficient is eps^4 Hlog*log(g); reference part cancels exactly.
        constant_part=source_jet-eps**4*r["Hlog"]*numerator
        log_part=-eps**4*r["Hlog"]*perturb
        checks[f"complete_source_reference_cancellation_piece{index}"]=constant_part==0
        primitive=log_part.primitive();pressure=primitive-primitive.at(0,left)+pressure_history
        pressure_history=pressure.at(0,right)
        expected=eps**4*r["Hlog"]*f["D"] if index==2 else means.Poly()
        checks[f"pointwise_compensated_pressure_piece{index}"]=pressure==expected
        checks[f"all_nonpressure_axial_velocity_jets_lower_piece{index}"]=all(v.at(1,0)==0 for v in (gamma,gamma.diff(1),gamma.diff(1).diff(1),B,B.diff(1)))
        axial_moment += old.integrate(R*eps*pressure,left,right)
        # Radial residual coefficient = derivative of compensated pressure - complete source.
        force_constant=-source_jet
        force_log=pressure.diff(0)
        radial_log_integral += old.integrate(force_log,left,right)
        checks[f"radial_compensation_identity_piece{index}"]=force_constant==-eps**4*r["Hlog"]*numerator and force_log==-eps**4*r["Hlog"]*perturb
        # Only the original three pieces need endpoint log symbols for the constant force.
        if index in (0,4):
            radial_integral=eadd(radial_integral,adapt(old.integrate_laurent(force_constant,left,right,ll,lr)))
    # Gap contribution is integrated on its original endpoints, avoiding new log symbols.
    gap_force=eps**4*gap[5]
    radial_integral=eadd(radial_integral,adapt(old.integrate_laurent(gap_force,gap[0],gap[1],gap[2],gap[3])))
    radial_integral=eadd(radial_integral,expr(g=radial_log_integral.scalar()))
    checks["radial_high_jet_total_nonzero"]=radial_integral==expr(g=-eps**4*r["Hlog"])
    checks["complete_pressure_right_endpoint_zero"]=pressure_history==0
    checks["integrated_axial_high_jet_zero_but_pointwise_nonzero"]=axial_moment==0 and bool(eps**5*r["Hlog"]*f["D"])
    return dict(checks=checks,radial_integral=display(radial_integral),
                pointwise_axial_coefficient_per_logg_at_midpoint=str((eps**5*r["Hlog"]*f["D"]).at(0,(f["lo"]+f["hi"])/2).scalar()),
                scope="exact full averaged-source response; same-moment comparison density, not an actual compensator substitution")


def moving_review():
    R,Z,T,eps=means.R,means.Z,means.T,means.EPS
    bump=(R-1)**6*(2-R)**6;rho0=bump/bump.integral().scalar()
    rho1=(R-Q(3,2))*rho0/100
    star=rho0+(Z+Z**2+Z**3+T+Z*T)*rho1
    rho=star+rho1/7
    F=1+Z+Z**2+Z**3+Z**4+T+T**2*Z**2
    B=eps*(F*star).diff(1)
    source=eps*B.diff(2)+eps**3*B.diff(1).diff(1)
    P=eps**2*F.diff(2).diff(1)+eps**4*F.diff(1).diff(1).diff(1)
    terms=[eps**2*star.diff(1)*F.diff(2),eps**2*star.diff(2)*F.diff(1),eps**2*star.diff(1).diff(2)*F,
           3*eps**4*star.diff(1)*F.diff(1).diff(1),3*eps**4*star.diff(1).diff(1)*F.diff(1),eps**4*star.diff(1).diff(1).diff(1)*F]
    compensated=source-rho*P;commutator=sum(terms,means.Poly())
    primitive=compensated.primitive();p=primitive-primitive.at(0,1)
    matched=source-star*P
    mp=matched.primitive();matched_pressure=mp-mp.at(0,1)
    cr=(R**2*rho).integral()/2
    return dict(checks=dict(
        moving_reference_and_actual_mass_one=star.integral()==rho.integral()==1,
        complete_uncompensated_source_total=source.integral()==P,
        full_shape_commutator_identity=compensated==(star-rho)*P+commutator,
        every_moving_shape_term_nonzero=all(bool(t) for t in terms),
        full_compensated_source_zero_total=compensated.integral()==0,
        full_pressure_compact=p.at(0,1)==p.at(0,2)==0,
        full_pressure_radial_identity=p.diff(0)-source==-rho*P,
        entire_moving_pressure_moment=eps*(R*p).integral().diff(1)==eps*(cr*P-(R**2*source).integral()/2).diff(1),
        shape_match_leaves_commutator=matched==commutator and bool(matched_pressure.diff(1)),
        moving_compensator_product_nonzero=bool(cr.diff(1)*P),
    ),scope="manufactured moving unit-mass polynomial densities for the exact product identity; not actual physical reference profiles")


def primitive_bound_review():
    # Exact signed masses of two separated nonnegative bumps saturate the half-L1 bound.
    R=means.R;b=(R-1)**4*(2-R)**4;mass=b.integral().scalar()
    translated=(R-4)**4*(5-R)**4
    positive=old.integrate(Q(3,7)*b/mass,Q(1),Q(2)).scalar()
    negative=old.integrate(-Q(3,7)*translated/mass,Q(4),Q(5)).scalar()
    total_variation=positive-negative
    return dict(checks=dict(zero_total=positive+negative==0,
                            half_L1_bound_sharp_on_gap=positive==total_variation/2,
                            reference_radial_mass_cannot_vanish=mass>0),
                gap_primitive=str(positive),L1=str(total_variation),
                scope="disjoint normalized smooth-bump mass control; the arbitrary-source estimate is proved analytically")


def review():
    f=fixture()
    parts={"reference_and_countercontrol":dict(checks=f["checks"],controls=f["controls"]),
           "pointwise_operator":operator_review(f),"moving_density":moving_review(),
           "primitive_bound":primitive_bound_review(),"pressure_tails":means.pressure_tail_review(),
           "complete_covariance":means.covariance_review()}
    return dict(checks={p+":"+k:v for p,part in parts.items() for k,v in part["checks"].items()},parts=parts)
