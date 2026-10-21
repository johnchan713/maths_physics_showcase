"""Exact current-state and explicit norm diagnostics; no actual pulse sampling."""
from fractions import Fraction as Q
import importlib.util
from math import comb
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("dynamic_mean_polynomials",HERE.parent/"complete_covariance_initial_means/means.py")
means=importlib.util.module_from_spec(spec)
spec.loader.exec_module(means)
KEYS=("b","G","old_B","old_v","old_g","B","u","g","Drr","Drt","Drz","Dtt","Dtz","Dzz","c_rho")


def budget(m,eps,a,d,norms):
    """Bounds for order m from field norms through m+2; exact rational only."""
    if type(m) is not int or m<0:
        raise ValueError("Nonnegative integer order required")
    if set(norms)!=set(KEYS):
        raise ValueError("All base, old means, increments and six covariance entries required")
    if any(type(x) not in (int,Q) for x in (eps,a,d,*norms.values())):
        raise TypeError("Exact rational bounds required")
    if not 0<eps<=1 or not 0<a<d or any(x<0 for x in norms.values()):
        raise ValueError("Positive shell/epsilon and nonnegative norm bounds required")
    eps,a,d=Q(eps),Q(a),Q(d)
    n={k:Q(v) for k,v in norms.items()}
    x,y,z=n["B"],n["u"],n["g"]
    A,C=n["b"]+n["old_B"],n["G"]+n["old_g"]
    ss=2*A*x+x*x+n["Drr"]
    hh=A*z+C*x+x*z+n["Drz"]
    tt=2*n["old_v"]*y+y*y+n["Dtt"]
    kk=n["old_g"]*y+n["old_v"]*z+y*z+n["Dtz"]
    zz=2*n["old_g"]*z+z*z+n["Dzz"]
    M={j:(d-a)*max(a**j,d**j) for j in (-1,0,1,2)}
    tm=eps*(m+1)+eps**3*(m+1)*(m+2)
    p=M[0]*tm*x+M[-1]*(ss+tt)+eps*(m+1)*M[0]*hh
    q=M[2]*tm*x+M[1]*(ss+tt)+eps*(m+1)*M[2]*hh
    jt,jz=M[2]*kk,M[1]*zz+q/2
    return dict(order=m,input_order=m+2,P=p,Qrem=q,Jtheta=jt,Jz=jz,
                residual_order=m-1 if m else None,
                Mz=eps*m*(jz+n["c_rho"]*p) if m else None,
                Mtheta=eps*m*jt if m else None,
                Drt_kept_for_full_target=n["Drt"],
                terms=dict(S=ss,H=hh,U=tt,K=kk,Zf=zz,tau=tm))


def majorant(f,n,a=Q(1),d=Q(2)):
    """Provable coefficient l1 bound on normalized Z,T jets, |Z|,|T|,|U|<=1."""
    if type(n) is not int or n<0:
        raise ValueError("Nonnegative derivative order required")
    if any(z<0 or t<0 or u<0 for (r,z,t,u) in f.terms):
        raise ValueError("Only polynomial slow/auxiliary variables allowed")
    return sum((abs(c)*comb(z,i)*comb(t,j)*max(a**r,d**r)
        for (r,z,t,u),c in f.terms.items() for i in range(min(z,n)+1)
        for j in range(min(t,n-i)+1)),Q(0))


def fixture():
    R,Z,T,U,eps=means.R,means.Z,means.T,means.U,means.EPS
    f=(R-1)**4*(2-R)**4
    psi0=R*f*(1+Z+Z**3)*(1+T+U)
    old=(-eps*psi0.diff(1),(R**4*f).diff(0).shift(-2)*(1+Z+T+U),means.radial(psi0))
    psi=R*f*(1+2*Z+Z**3)*(1+T+T**2)
    B,u,g=-eps*psi.diff(1),(R**4*f).diff(0).shift(-2)*(1+2*Z+T),means.radial(psi)
    b,V,G=eps*R*f*(1+Z+T),R**2*(2+Z),(1+Z)*R*f
    _,_,_,W=means.fixture()
    ratio=Q(2,7)
    Wnew=[[(1+ratio)**2*w for w in row] for row in W]
    D=[[ (Wnew[i][j]-W[i][j]).haar() for j in range(3)] for i in range(3)]
    B0,v0,g0=(f.haar() for f in old)
    new=tuple(f+h for f,h in zip(old,(B,u,g)))
    gr0=means.averaged_gr((b,V,G),old,W[0][0],W[1][1],W[0][2])
    gr1=means.averaged_gr((b,V,G),new,Wnew[0][0],Wnew[1][1],Wnew[0][2])
    rem=gr1-gr0-(2*V*u).shift(-1)
    S=2*(b+B0)*B+B**2+D[0][0]
    H=(b+B0)*g+(G+g0)*B+B*g+D[0][2]
    UU=2*v0*u+u**2+D[1][1]
    K=g0*u+v0*g+u*g+D[1][2]
    Zf=2*g0*g+g**2+D[2][2]
    tau=lambda x:eps*x.diff(2)+eps**3*x.diff(1).diff(1)
    expected=tau(B)-means.radial(S)-eps*H.diff(1)+UU.shift(-1)+eps*(B.diff(0).diff(0)+B.diff(0).shift(-1)-B.shift(-2))
    P=tau(B.integral())+(UU-S).shift(-1).integral()-eps*H.integral().diff(1)
    Qrem=tau((R**2*B).integral())+(R*(S+UU)).integral()-eps*(R**2*H).integral().diff(1)
    Jt=(R**2*K).integral(); Jz=(R*Zf).integral()-Qrem/2
    jt=lambda v,w:(R**2*(G*v[1]+V*v[2]+v[1]*v[2]+w[1][2]).haar()).integral()
    jz=lambda v,w,gr:(R*(2*G*v[2]+v[2]**2+w[2][2]).haar()).integral()-(R**2*gr).integral()/2
    rho=f/f.integral().scalar()+Z*(R-Q(3,2))*f/f.integral().scalar()
    cr=(R**2*rho).integral()/2
    # Compensate the nonlinear remainder and, separately, the FULL updated source.
    primitive=(rem-rho*P).primitive();pm=primitive-primitive.at(0,1)
    fullP=gr1.integral()
    fullprimitive=(gr1-rho*fullP).primitive()
    fullpressure=fullprimitive-fullprimitive.at(0,1)
    fullaxial=(2*G*new[2]+new[2]**2+Wnew[2][2]).haar()
    fields=dict(zip(KEYS,(b,G,B0,v0,g0,B,u,g,D[0][0],D[0][1],D[0][2],D[1][1],D[1][2],D[2][2],cr)))
    checks=dict(
        old_and_increment_divergence=means.radial(old[0])+eps*old[2].diff(1)==means.radial(B)+eps*g.diff(1)==0,
        old_and_new_velocity_moments=all((R**2*v[1]).integral().haar()==(R*v[2]).integral().haar()==0 for v in (old,new)),
        old_Haar_means_nonzero=all(bool(x) for x in (B0,v0,g0)),
        full_current_radial_difference=rem==expected,
        plain_moment_identity=P==rem.integral(),
        weighted_moment_identity=Qrem==(R**2*rem).integral(),
        theta_current_difference=jt(new,Wnew)-jt(old,W)-(R**2*(G*u+V*g)).integral()==Jt,
        axial_current_difference=jz(new,Wnew,gr1)-jz(old,W,gr0)-(2*R*G*g-R*V*u).integral()==Jz,
        complete_pressure_both_ends=pm.at(0,1)==pm.at(0,2)==0,
        whole_pressure_axial_moment=(R*Zf+R*pm).integral().diff(1)==(Jz+cr*P).diff(1),
        full_updated_pressure_both_ends=fullpressure.at(0,1)==fullpressure.at(0,2)==0,
        full_updated_pressure_moving_moment=(R*(fullaxial+fullpressure)).integral().diff(1)==(jz(new,Wnew,gr1)+cr*fullP).diff(1),
        moving_compensator_term_nonzero=bool(cr.diff(1)*P),
        omitted_old_mean_pressure_detected=bool((2*v0*u).shift(-1).integral()),
        omitted_old_axial_cross_detected=bool((2*R*g0*g).integral()),
        all_six_covariance_changes_nonzero=all(bool(D[i][j]) for i in range(3) for j in range(i,3)),
        complete_wave_rescaling_identity=all(Wnew[i][j]-W[i][j]==(2*ratio+ratio**2)*W[i][j] for i in range(3) for j in range(3)),
        Drt_changes_target_flux=bool(D[0][1]),
    )
    rows=[]
    for m in (0,1,3):
        ns={k:majorant(v,m+2) for k,v in fields.items()}
        upper=budget(m,eps,Q(1),Q(2),ns)
        actual={k:majorant(v,m) for k,v in dict(P=P,Qrem=Qrem,Jtheta=Jt,Jz=Jz).items()}
        for key,v in actual.items():
            checks[f"norm_bound_{key}_m{m}"]=v<=upper[key]
        if m:
            checks[f"whole_axial_residual_bound_m{m-1}"]=majorant(eps*(Jz+cr*P).diff(1),m-1)<=upper["Mz"]
            checks[f"angular_residual_bound_m{m-1}"]=majorant(eps*Jt.diff(1),m-1)<=upper["Mtheta"]
        rows.append(dict(order=m,input_order=m+2,bounds={k:str(upper[k]) for k in ("P","Qrem","Jtheta","Jz","Mz","Mtheta")},
                         coefficient_majorants={k:str(v) for k,v in actual.items()}))
    return dict(checks=checks,norm_rows=rows,scope="finite compact polynomial box; no actual five-map or all-scale suprema")


def correlation_controls():
    rows=[];checks={}
    for k in (4,12,40):
        eps=Q(1,2**k)
        # mu=1, alpha=2; energy divided by epsilon^(2alpha).
        cross=2*eps**3+eps**4
        ratio=cross/eps**4
        checks[f"old_mean_inverse_epsilon_rate_k{k}"]=ratio==2/eps+1
        # nu=1/2, chi=3/2, use epsilon=t^2 for rational amplitudes.
        t=Q(1,2**k)
        old,new=t,t**3
        delta=2*old*new+new**2
        checks[f"wave_cross_not_square_k{k}"]=delta/new**2==2/t**2+1
        checks[f"same_phase_covariance_cross_retained_k{k}"]=delta>new**2
        rows.append(dict(k=k,old_mean_quadratic_ratio=str(ratio),wave_square_ratio=str(delta/new**2)))
    checks["ratios_strictly_increase"]=all(Q(a["old_mean_quadratic_ratio"])<Q(b["old_mean_quadratic_ratio"]) for a,b in zip(rows,rows[1:]))
    # Exact circle averages: cos^2=sin^2=1/2, cos sin=0.
    w,q=Q(3,5),Q(2,7)
    bilinear=lambda x,y:x[0]*y[0]+x[1]*y[1]
    same_cross=bilinear((w,Q(0)),(q,Q(0)))
    orthogonal_cross=bilinear((w,Q(0)),(Q(0),q))
    checks["orthogonal_phase_special_cross_zero"]=orthogonal_cross==0
    checks["same_phase_generic_cross_nonzero"]=same_cross!=0
    return dict(checks=checks,rate_rows=rows,scope="exact rational sequences corroborate the smooth continuum countercontrols")


def review():
    parts={"current_update":fixture(),"correlations":correlation_controls(),"pressure_tails":means.pressure_tail_review()}
    return dict(checks={p+":"+k:v for p,r in parts.items() for k,v in r["checks"].items()},parts=parts)
