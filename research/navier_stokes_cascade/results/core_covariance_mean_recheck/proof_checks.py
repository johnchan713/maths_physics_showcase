"""Exact manufactured controls for the analytical projected mean-defect bounds."""
from fractions import Fraction as Q
import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


means = load("mean_defect_polynomials", HERE.parent / "complete_covariance_initial_means/means.py")


def dot(x, y, weights):
    if not len(x) == len(y) == len(weights) or not x:
        raise ValueError("Nonempty equal-length vectors required")
    if any(type(v) not in (int, Q) for v in (*x, *y, *weights)):
        raise TypeError("Exact rational inputs required")
    if any(w <= 0 for w in weights):
        raise ValueError("Positive weights required")
    return sum((w*a*b for w, a, b in zip(weights, x, y)), Q(0))


def project(target, null, weights):
    d = dot(null, null, weights)
    if not d:
        raise ValueError("Nonzero null vector required")
    coefficient = dot(target, null, weights)/d
    p = [t-coefficient*n for t, n in zip(target, null)]
    return p, dot(p, p, weights)


def dynamic_review():
    R, Z, T, U, eps = means.R, means.Z, means.T, means.U, means.EPS
    bump = (R-1)**4*(2-R)**4
    psi = R*bump*(1+Z+Z**3)*(1+T+T**2)
    B, g = -eps*psi.diff(1), means.radial(psi)
    u = (R**4*bump).diff(0).shift(-2)*(1+Z+T)
    b, V, G = eps*R*bump*(1+Z+T), R**2*(2+Z), means.Poly()
    fast = (eps*R*bump*U, R*bump*U, R**2*bump*U)
    _, _, _, W = means.fixture()
    new = tuple(f+s for f, s in zip(fast, (B, u, g)))
    old_gr = means.averaged_gr((b, V, G), fast, W[0][0], W[1][1], W[0][2])
    new_gr = means.averaged_gr((b, V, G), new, W[0][0], W[1][1], W[0][2])
    direct = new_gr-old_gr-(2*V*u).shift(-1)
    S, H = 2*b*B+B**2, b*g+B*g
    radial_viscosity = B.diff(0).diff(0)+B.diff(0).shift(-1)-B.shift(-2)
    axial_viscosity = eps**3*B.diff(1).diff(1)
    expected = eps*B.diff(2)-means.radial(S)-eps*H.diff(1)+(u**2).shift(-1)+eps*radial_viscosity+axial_viscosity
    tau = lambda f: eps*f.diff(2)+eps**3*f.diff(1).diff(1)
    cp_time, cp_axvis = eps*B.integral().diff(2), axial_viscosity.integral()
    cp_axial = -eps*H.integral().diff(1)
    cp = tau(B.integral())-S.shift(-1).integral()+cp_axial
    cq = tau((R**2*B).integral())+(R*S).integral()-eps*(R**2*H).integral().diff(1)
    eu, fu, eg = (u**2).shift(-1).integral(), (R*u**2).integral(), (R*g**2).integral()
    # Independently reconstruct the entire compensated pressure primitive.
    mass=bump.integral().scalar()
    rho=bump/mass+Z*(R-Q(3,2))*bump/mass
    crho=(R**2*rho).integral()/2
    primitive=(direct-rho*direct.integral()).primitive()
    pm=primitive-primitive.at(0,1)
    pressure_moment=(R*pm).integral()
    # Independent complete theta/axial defect differences, before linear rows.
    jt = lambda m: (R**2*(V*m[2]+m[1]*m[2]).haar()).integral()
    jz = lambda m, gr: (R*(m[2]**2).haar()).integral()-(R**2*gr).integral()/2
    later = (fast[0], fast[1]+R*bump, fast[2])
    later_new = tuple(f+s for f, s in zip(later, (B, u, g)))
    bad = means.averaged_gr((b,V,G),later_new,W[0][0],W[1][1],W[0][2])-means.averaged_gr((b,V,G),later,W[0][0],W[1][1],W[0][2])-(2*V*u).shift(-1)
    return dict(checks={
        "compact_endpoint_second_jets": all(f.diff(0).diff(0).at(0,e)==0 for f in (B,u,g) for e in (1,2)),
        "exact_slow_divergence": means.radial(B)+eps*g.diff(1)==0,
        "both_velocity_moments_zero": (R**2*u).integral()==(R*g).integral()==0,
        "zero_fast_haar_means": all(f.haar()==0 for f in fast),
        "complete_operator_remainder": direct==expected,
        "plain_pressure_identity": direct.integral()==eu+cp,
        "weighted_pressure_identity": (R**2*direct).integral()==fu+cq,
        "compensator_unit_mass": rho.integral()==1,
        "full_pressure_primitive_compact": pm.at(0,1)==pm.at(0,2)==0,
        "whole_axial_pressure_moment": (eg+pressure_moment).diff(1)==(eg-(fu+cq)/2+crho*(eu+cp)).diff(1),
        "moving_compensator_derivative_omission_detected": bool(crho.diff(1)*(eu+cp)),
        "radial_vector_viscosity_plain_cancels": radial_viscosity.integral()==0,
        "radial_vector_viscosity_weighted_cancels": (R**2*radial_viscosity).integral()==0,
        "radial_vector_viscosity_pointwise_nonzero": bool(radial_viscosity),
        "time_omission_detected": bool(cp_time),
        "axial_viscosity_omission_detected": bool(cp_axvis),
        "axial_transport_omission_detected": bool(cp_axial),
        "angular_nonlinear_identity": jt(new)-jt(fast)-(R**2*V*g).integral()==(R**2*u*g).integral(),
        "axial_nonlinear_identity": jz(new,new_gr)-jz(fast,old_gr)+(R*V*u).integral()==eg-(fu+cq)/2,
        "later_nonzero_old_mean_breaks_initial_formula": bad!=expected,
        "Haar_cross_cancels_only_after_average": bool(fast[1]*u) and (fast[1]*u).haar()==0,
        "Haar_square_retained": bool((fast[1]**2).haar()),
        "stationary_positive_angular_energy": eu.at(1,0).at(2,0).scalar()>0,
    }, controls={"E_u_at_origin":str(eu.at(1,0).at(2,0).scalar()),
                 "E_g_at_origin":str(eg.at(1,0).at(2,0).scalar())},
       fixture="finite-regularity compact polynomial, arbitrary manufactured V; not actual five-map targets")


def projection_review():
    radii = [Q(1),Q(2),Q(3)]
    weights = [1/r for r in radii]
    V = [1/r for r in radii]
    null = [r**3 for r in radii]
    p, dp = project([2*v for v in V],null,weights)
    pressure = Q(7,5)
    u = [-pressure*v/dp for v in p]
    eu = dot(u,u,weights)
    raw = dot([2*v for v in V],[2*v for v in V],weights)
    checks = {
        "angular_null_moment": dot(u,null,weights)==0,
        "angular_target": dot(u,[2*v for v in V],weights)==-pressure,
        "projected_angular_bound_saturated": eu==pressure**2/dp,
        "unprojected_bound_strictly_weaker": eu>pressure**2/raw,
    }
    controls = []
    for m in (1,4,20,80):
        q = 1-Q(1,2**m)
        f, w, one = [Q(1),q], [Q(1),Q(2)], [Q(1),Q(1)]
        pz, dz = project(f,one,w)
        target = Q(3,7)
        g = [-target*v/dz for v in pz]
        pair = w[0]*w[1]/sum(w)*(1-q)**2
        checks[f"axial_pairwise_variance_m{m}"] = dz==pair
        checks[f"axial_zero_flux_m{m}"] = dot(g,one,w)==0
        checks[f"axial_target_m{m}"] = dot(g,f,w)==-target
        checks[f"axial_inverse_square_energy_m{m}"] = dot(g,g,w)==target**2/dz
        controls.append({"m":m,"D_z":str(dz),"energy":str(dot(g,g,w))})
    pzero,dzero=project([Q(1),Q(1)],[Q(1),Q(1)],[Q(1),Q(2)])
    checks["zero_tilt_constant_target_projection_zero"] = dzero==0 and pzero==[0,0]
    return dict(checks=checks, small_denominators=controls,
                scope="finite rational Hilbert-space controls, not continuum quadrature")


def continuum_moments():
    # lambda=1/4, a_eff=1, I=[1,16]; the fractional powers have rational endpoints.
    i0, mixed, norm = Q(255,2),Q(42),Q(15)
    dz = norm-mixed**2/i0
    angular_norm = Q(1- Q(1,16**3),3)
    radial_sixth = Q(16**6-1,6)
    dp = 4*(angular_norm-mixed**2/radial_sixth)
    return dict(checks={"exact_continuum_axial_variance":dz==Q(99,85),
                        "both_continuum_denominators_positive":dz>0 and dp>0,
                        "projected_axial_bound_stronger_than_unprojected":1/dz>1/norm},
                D_z=str(dz),D_P=str(dp),lambda_value="1/4",interval=[1,16])


def joint_review():
    # Rational samples of V=R^(-3/2), not a continuum quadrature.
    r=[Q(1),Q(4),Q(9)]
    V=[Q(1),Q(1,8),Q(1,27)]
    wa,wz=[1/x for x in r],r
    angular=[2*v for v in V]
    pa,dp=project(angular,[x**3 for x in r],wa)
    pz,dz=project([x*v for x,v in zip(r,V)],[1,1,1],wz)
    p,j=Q(2,5),Q(3,7)
    u,g=[-p*x/dp for x in pa],[-j*x/dz for x in pz]
    eu,fu,eg=dot(u,u,wa),dot(u,u,wz),dot(g,g,wz)
    # Retained terms can cancel these TWO defects. This is not a PDE solution.
    cp,cq=-eu,2*eg-fu
    pnew,jznew=eu+cp,eg-(fu+cq)/2
    eliminated=jznew+(max(r)**2*(pnew-cp)+cq)/2
    return dict(checks={
        "compact_hull_energy_comparison":fu<=max(r)**2*eu,
        "two_defects_can_cancel_with_dynamic_terms":pnew==jznew==0,
        "joint_eliminated_budget":eliminated>=j**2/dz,
        "pressure_cancellation_requires_negative_C_P":-cp>=p**2/dp>0,
        "joint_zero_defect_gate":cq-max(r)**2*cp>=2*j**2/dz,
        "stationary_zero_dynamic_terms_fail_joint_gate":Q(0)<2*j**2/dz,
    },scope="finite-space necessity control, no dynamic realizability or third-defect closure")


def review():
    parts = {"dynamic":dynamic_review(), "projection":projection_review(),
             "continuum":continuum_moments(), "joint_budget":joint_review(),
             "complete_covariance":means.covariance_review(),
             "pressure_tails":means.pressure_tail_review()}
    checks = {name+":"+key:value for name,part in parts.items() for key,value in part["checks"].items()}
    return {"checks":checks,"parts":parts,
            "analytical_bounds_in_readme_not_inferred_from_fixtures":True}
