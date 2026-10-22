"""Exact coupled moment budgets, complete covariance and Haar correlations."""
from fractions import Fraction as Q
import importlib.util
from pathlib import Path

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


cycle=load("coupled_scalar_cycle",HERE.parent/"core_scalar_cycle_budget/proof_checks.py")
dynamic=cycle.dynamic;means=dynamic.means
R,Z,T,U=means.R,means.Z,means.T,means.U
EPS=means.EPS;P=cycle.ScalePolynomial


def geometry(a,d):
    a,d=cycle.exact(a),cycle.exact(d)
    if not 0<a<d:raise ValueError("Positive ordered integration hull required")
    return {j:(d-a)*max(a**j,d**j) for j in (-1,0,1,2)}


def order(m):
    if type(m) is not int or m<0:raise ValueError("Nonnegative integer order required")
    return m


def wave_budget(m,M,d0,d1):
    order(m)
    return dict(P=2*M[-1]*d0+(m+1)*M[0]*d1,
                Jz=2*M[1]*d0+Q(m+1,2)*M[2]*d1)


def fast_budget(m,M,Ot,Or,Ft,Fr):
    order(m)
    s=2*Or*Fr+Fr*Fr;h=Or*Ft+Ot*Fr+Fr*Ft
    u=2*Ot*Ft+Ft*Ft;z=u
    p=M[-1]*(s+u)+(m+1)*M[0]*h
    q=M[1]*(s+u)+(m+1)*M[2]*h
    return dict(P=p,Qrem=q,Jz=M[1]*z+q*Q(1,2))


def slow_budget(m,M,b,r,a,X,U,Z):
    order(m)
    s=2*(b+r)*X+X*X;h=(b+r)*Z+a*X+X*Z
    t=2*a*U+U*U;z=2*a*Z+Z*Z
    tm=(m+1)+(m+1)*(m+2)
    p=M[0]*tm*X+M[-1]*(s+t)+(m+1)*M[0]*h
    q=M[2]*tm*X+M[1]*(s+t)+(m+1)*M[2]*h
    jt=M[2]*(a*(U+Z)+U*Z)
    return dict(P=p,Qrem=q,Jtheta=jt,Jz=M[1]*z+q*Q(1,2))


def phi(poly):
    p=poly.primitive()
    return p-p.at(0,1)


def fields():
    f=(R-1)**6*(2-R)**6
    psi0=R*f*(1+Z+Z**3)*(1+T)*(1+U)
    old=(-EPS*psi0.diff(1),(R**4*f).diff(0).shift(-2)*(1+Z+T)*(1+U),means.radial(psi0))
    base=(EPS*R*f*(1+Z+T),R**2*(2+Z),means.Poly())
    _,_,_,W=means.fixture()
    return f,old,base,W


def scalar_z(base,mean,W,gr):
    G=base[2];g=mean[2]
    return (R*(2*G*g+g*g+W[2][2]).haar()).integral()-(R**2*gr).integral()/2


def wave_review():
    f,old,base,W=fields();ratio=Q(2,7)
    Wnew=[[(1+ratio)**2*x for x in row] for row in W]
    D=[[(Wnew[i][j]-W[i][j]).haar() for j in range(3)] for i in range(3)]
    gr0=means.averaged_gr(base,old,W[0][0],W[1][1],W[0][2])
    gr1=means.averaged_gr(base,old,Wnew[0][0],Wnew[1][1],Wnew[0][2])
    delta=gr1-gr0
    dP=(D[1][1]-D[0][0]).shift(-1).integral()-EPS*D[0][2].integral().diff(1)
    dJz=(R*(D[2][2]-(D[0][0]+D[1][1])/2)).integral()+EPS*(R**2*D[0][2]).integral().diff(1)/2
    checks={
        "full_current_wave_source_difference":delta==-means.radial(D[0][0])-EPS*D[0][2].diff(1)+D[1][1].shift(-1),
        "current_wave_plain_pressure_moment":delta.integral()==dP,
        "current_wave_axial_scalar_moment":scalar_z(base,old,Wnew,gr1)-scalar_z(base,old,W,gr0)==dJz,
        "all_six_covariance_changes_nonzero":all(bool(D[i][j]) for i in range(3) for j in range(i,3)),
        "full_covariance_cross_plus_square":all(D[i][j]==(2*ratio+ratio**2)*W[i][j].haar() for i in range(3) for j in range(3)),
        "axial_velocity_covariance_enters_Jz":bool((R*D[2][2]).integral()),
        "radial_axial_derivative_enters_moments":bool(D[0][2].integral().diff(1)) or bool((R**2*D[0][2]).integral().diff(1)),
    }
    M=geometry(1,2)
    for m in (0,1,3):
        d0=max(dynamic.majorant(D[i][j],m) for i in range(3) for j in range(i,3))
        d1=max(dynamic.majorant(D[i][j],m+1) for i in range(3) for j in range(i,3))
        upper=wave_budget(m,M,d0,d1)
        checks[f"wave_pressure_budget_m{m}"]=dynamic.majorant(dP,m)<=upper["P"]
        checks[f"wave_axial_scalar_budget_m{m}"]=dynamic.majorant(dJz,m)<=upper["Jz"]
    return dict(checks=checks,scope="independent current covariance rescaling, not a canonical wave inverse realization")


def fast_review():
    f,old,base,W=fields()
    psif=R*f*(1+Z**2+Z**4)*(1+T+T**2)*U
    fast=(-EPS*psif.diff(1),(R**4*f).diff(0).shift(-2)*(1+Z+T)*U,means.radial(psif))
    new=tuple(x+y for x,y in zip(old,fast))
    beta,v,gamma=old;B,u,g=fast
    sf=(2*beta*B+B*B).haar();hf=(beta*g+gamma*B+B*g).haar()
    uf=(2*v*u+u*u).haar();zf=(2*gamma*g+g*g).haar()
    gr0=means.averaged_gr(base,old,W[0][0],W[1][1],W[0][2])
    gr1=means.averaged_gr(base,new,W[0][0],W[1][1],W[0][2]);delta=gr1-gr0
    dP=(uf-sf).shift(-1).integral()-EPS*hf.integral().diff(1)
    dQ=(R*(sf+uf)).integral()-EPS*(R**2*hf).integral().diff(1)
    dJz=(R*zf).integral()-dQ/2
    rho=f/f.integral().scalar()+Z*(R-Q(3,2))*f/f.integral().scalar()
    cr=(R**2*rho).integral()/2;p=phi(delta-rho*dP)
    checks={
        "all_actual_fast_components_zero_Haar":all(x.haar()==0 for x in fast),
        "full_Haar_fast_source_is_complete_nonlinear_change":delta==-means.radial(sf)-EPS*hf.diff(1)+uf.shift(-1),
        "fast_plain_moment_exact":delta.integral()==dP,
        "fast_weighted_source_moment_exact":(R**2*delta).integral()==dQ,
        "fast_current_axial_scalar_exact":scalar_z(base,new,W,gr1)-scalar_z(base,old,W,gr0)==dJz,
        "old_new_radial_Haar_correlation_nonzero":bool((beta*B).haar()),
        "old_new_axial_radial_Haar_correlation_nonzero":bool((gamma*B).haar()) and bool((beta*g).haar()),
        "products_of_Haar_means_would_erase_crosses":beta.haar()*B.haar()==gamma.haar()*B.haar()==0,
        "fast_compensated_pressure_compact":p.at(0,1)==p.at(0,2)==0,
        "whole_fast_moving_pressure_moment":EPS*(R*(zf+p)).integral().diff(1)==EPS*(dJz+cr*dP).diff(1),
        "fast_moving_compensator_derivative_nonzero":bool(cr.diff(1)*dP),
    }
    # Independent coefficient norms supply sufficient bounds for the exact products.
    M=geometry(1,2)
    for m in (0,1,2):
        n=m+1
        Ot=max(dynamic.majorant(v,n),dynamic.majorant(gamma,n));Or=dynamic.majorant(beta,n)
        Ft=max(dynamic.majorant(u,n),dynamic.majorant(g,n));Fr=dynamic.majorant(B,n)
        upper=fast_budget(m,M,Ot,Or,Ft,Fr)
        checks[f"fast_pressure_norm_bound_m{m}"]=dynamic.majorant(dP,m)<=upper["P"]
        checks[f"fast_axial_norm_bound_m{m}"]=dynamic.majorant(dJz,m)<=upper["Jz"]
    return dict(checks=checks,scope="manufactured Haar coefficient fields and zero-Haar curl algebra, not the full physical fast inverse")


def slow_review():
    f,old,base,W=fields();psi=R*f*(1+2*Z+Z**3)*(1+T+T**2)
    B=-EPS*psi.diff(1);u=(R**4*f).diff(0).shift(-2)*(1+2*Z+T);g=means.radial(psi)
    new=tuple(x+y for x,y in zip(old,(B,u,g)))
    b,V,G=base;beta,v,gamma=old;B0=beta.haar();v0=v.haar();g0=gamma.haar()
    S=2*(b+B0)*B+B*B;H=(b+B0)*g+g0*B+B*g
    Uu=2*v0*u+u*u;Zf=2*g0*g+g*g
    gr0=means.averaged_gr(base,old,W[0][0],W[1][1],W[0][2])
    gr1=means.averaged_gr(base,new,W[0][0],W[1][1],W[0][2])
    rem=gr1-gr0-(2*V*u).shift(-1)
    tau=lambda p:EPS*p.diff(2)+EPS**3*p.diff(1).diff(1)
    pm=tau(B.integral())+(Uu-S).shift(-1).integral()-EPS*H.integral().diff(1)
    qm=tau((R**2*B).integral())+(R*(S+Uu)).integral()-EPS*(R**2*H).integral().diff(1)
    jz=(R*Zf).integral()-qm/2
    rho=f/f.integral().scalar()+Z*(R-Q(3,2))*f/f.integral().scalar()
    cr=(R**2*rho).integral()/2;p=phi(rem-rho*pm)
    checks={
        "slow_full_source_before_moment_integration":rem==tau(B)-means.radial(S)-EPS*H.diff(1)+Uu.shift(-1)+EPS*(B.diff(0).diff(0)+B.diff(0).shift(-1)-B.shift(-2)),
        "slow_plain_remainder_moment":pm==rem.integral(),
        "slow_weighted_remainder_moment":qm==(R**2*rem).integral(),
        "slow_axial_scalar_difference_after_linear_row":scalar_z(base,new,W,gr1)-scalar_z(base,old,W,gr0)+(R*V*u).integral()==jz,
        "slow_divergence_exact":means.radial(B)+EPS*g.diff(1)==0,
        "slow_both_velocity_moments_zero":(R**2*u).integral()==(R*g).integral()==0,
        "radial_viscosity_plain_moment_zero":(B.diff(0).diff(0)+B.diff(0).shift(-1)-B.shift(-2)).integral()==0,
        "radial_viscosity_weighted_moment_zero":(R**2*(B.diff(0).diff(0)+B.diff(0).shift(-1)-B.shift(-2))).integral()==0,
        "axial_viscosity_retained_nonzero":bool(EPS**3*B.integral().diff(1).diff(1)),
        "slow_pressure_compact_at_both_ends":p.at(0,1)==p.at(0,2)==0,
        "whole_slow_pressure_moment":EPS*(R*(Zf+p)).integral().diff(1)==EPS*(jz+cr*pm).diff(1),
        "whole_moving_product_derivative_required":bool(cr.diff(1)*pm),
    }
    M=geometry(1,2);C=2
    for m in (0,1,2):
        n=m+2;major=lambda x:dynamic.majorant(x,n)
        inputs=(major(b)/EPS,major(B0)/EPS**2,max(major(v0),major(g0))/EPS,
                major(B)/EPS**(C+1),major(u)/EPS**C,major(g)/EPS**C)
        upper=slow_budget(m,M,*inputs)
        checks[f"reduced_postslow_pressure_bound_m{m}"]=dynamic.majorant(pm,m)<=EPS**(C+1)*upper["P"]
        checks[f"reduced_postslow_axial_bound_m{m}"]=dynamic.majorant(jz,m)<=EPS**(C+1)*upper["Jz"]
        if m:
            moment=EPS*(jz+cr*pm).diff(1)
            checks[f"whole_reduced_axial_residual_bound_m{m-1}"]=dynamic.majorant(moment,m-1)<=EPS**(C+2)*m*(upper["Jz"]+dynamic.majorant(cr,m)*upper["P"])
    return dict(checks=checks,scope="exact remainder after current linear-row subtraction; not a numerical physical five-map solve or actual asymptotic norm certificate")


def polynomial_review():
    M=geometry(Q(1),Q(2));checks={};rows=[]
    b=P({0:2,1:1});r=P({0:3});a=P({1:Q(1,2)})
    X=P({0:1,2:2});u=P({1:3});z=P({0:2,1:1})
    # Compare the reduced polynomial with the independently inherited epsilon budget.
    for m in (0,2):
        reduced=slow_budget(m,M,b,r,a,X,u,z)
        for Svalue in (1,4,32):
            values=[p.evaluate(Svalue) for p in (b,r,a,X,u,z)]
            for eps in (Q(1,2),Q(1,8)):
                C=2;bv,rv,av,xv,uv,zv=values
                norms={k:Q(0) for k in dynamic.KEYS}
                norms.update(b=eps*bv,old_B=eps**2*rv,old_v=eps*av,old_g=eps*av,
                             B=eps**(C+1)*xv,u=eps**C*uv,g=eps**C*zv)
                raw=dynamic.budget(m,eps,Q(1),Q(2),norms)
                for key in ("P","Qrem","Jtheta","Jz"):
                    checks[f"unreduced_below_polynomial_{key}_m{m}_S{Svalue}_e{eps}"]=raw[key]<=eps**(C+1)*reduced[key].evaluate(Svalue)
            rows.append(dict(order=m,S=Svalue,ThetaP=str(reduced["P"].evaluate(Svalue)),ThetaZ=str(reduced["Jz"].evaluate(Svalue))))
    # In the combined ledger scalar output m uses incoming angular m+3;
    # the extra axial moment derivative increases that to m+4.
    for m in (0,2,5):
        scalar_fields=m+2;moment_fields=(m+1)+2
        checks[f"scalar_input_target_order_m{m}"]=scalar_fields+1==m+3
        checks[f"axial_moment_input_target_order_m{m}"]=moment_fields+1==m+4
    margins=cycle.exponent_margins(2-Q(2,100000))
    checks["all_coupled_outputs_supply_next_cycle_order"]=margins["next_input"]==Q(1,2)+Q(1,100000)>0
    kappa=Q(1,100000);C0=2-2*kappa;mu=1-kappa
    # Every affine exponent margin is minimized at C0 on the whole half-line.
    for label,margin,slope in (("wave_square",C0-1,1),
                              ("fast_old_cross",mu-kappa,0),
                              ("fast_square",C0-2*kappa,1),
                              ("slow_square",C0-1,1)):
        checks["whole_half_line_exponent_"+label]=margin>0 and slope>=0
    first=wave_budget(1,M,P({0:2,1:1}),P({0:3,2:2}))
    second=wave_budget(1,M,P({0:4}),P({1:2}))
    fast=fast_budget(1,M,a,r,u,X)
    incoming_p=P({0:1});incoming_z=P({1:2})
    preslow_p=incoming_p+first["P"]+second["P"]+fast["P"]
    preslow_z=incoming_z+first["Jz"]+second["Jz"]+fast["Jz"]
    for Svalue in (1,4,32):
        ev=lambda x:x.evaluate(Svalue)
        checks[f"both_wave_and_fast_pressure_budgets_accumulate_S{Svalue}"]=ev(preslow_p)==ev(incoming_p)+ev(first["P"])+ev(second["P"])+ev(fast["P"])
        checks[f"both_wave_and_fast_axial_budgets_accumulate_S{Svalue}"]=ev(preslow_z)==ev(incoming_z)+ev(first["Jz"])+ev(second["Jz"])+ev(fast["Jz"])
    return dict(checks=checks,rows=rows,PiP=preslow_p.record(),PiZ=preslow_z.record(),scope="positive polynomial majorants and explicit input-order ledger, not verification of the retained inverse constants")


def review():
    parts={"complete_wave":wave_review(),"fast_correlations":fast_review(),"slow_remainder":slow_review(),
           "polynomial_budget":polynomial_review(),"pressure_tails":means.pressure_tail_review(),
           "complete_covariance":means.covariance_review()}
    return dict(checks={p+":"+k:v for p,r in parts.items() for k,v in r["checks"].items()},parts=parts)
