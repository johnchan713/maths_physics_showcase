"""Exact scalar/projected target identities and higher-jet pressure routing."""
from fractions import Fraction as Q
import importlib.util
from math import comb, factorial
from pathlib import Path

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


dynamic=load("scalar_join_dynamic",HERE.parent/"core_dynamic_mean_recheck/proof_checks.py")
means=dynamic.means
R,Z,T,U=means.R,means.Z,means.T,means.U
EPS=means.EPS


def diff(poly,axis,n):
    for _ in range(n):poly=poly.diff(axis)
    return poly


def primitive(poly,e=0):
    p=(R**e*poly).primitive()
    return (p-p.at(0,1)).shift(-e)


def pure_t_norm(poly,n):
    return sum((abs(c)*comb(t,j)*max(Q(1)**r,Q(2)**r)
                for (r,z,t,u),c in poly.terms.items() for j in range(min(t,n)+1)),Q(0))


def norm(poly,n):return dynamic.majorant(poly,n)


def qshift(m,j):return factorial(m+j)//factorial(m)


def current_fixture():
    bump=(R-1)**6*(2-R)**6
    basepsi=R*bump*(1+Z+T)
    base=(-EPS*basepsi.diff(1),R**2*(2+Z+Z**2+T),means.radial(basepsi))
    psi=R*bump*(1+Z+Z**3)*(1+T+U)
    mean=(-EPS*psi.diff(1),(R**4*bump).diff(0).shift(-2)*(1+Z+T+U),means.radial(psi))
    b,V,G=base;beta,v,gamma=mean
    _,_,_,W=means.fixture()
    stress=R**2*bump*(1+Z+T)
    H=(G*v+V*gamma+gamma*v+W[1][2]).haar()
    J=(R**2*H).integral()
    flux=(b*v+beta*V+beta*v+W[0][1]-stress).haar()
    v0=v.haar();g0=gamma.haar()
    # Evaluate the complete angular equation before projecting its scalar moment.
    E=(EPS*v0.diff(2)+means.radial(flux,2)+EPS*H.diff(1)
       -EPS*(means.lap(v0)-v0.shift(-2)))
    M=(R**2*E).integral()
    mass=(R**2*bump).integral().scalar()
    center=(R**3*bump).integral().scalar()/mass
    sigma=bump*(1+(Z/13+T/17)*(R-center))/mass
    QE=E-sigma*M;Xi=-primitive(QE,2)
    checks={
        "current_Haar_mean_exact_divergence":means.radial(beta.haar())+EPS*gamma.haar().diff(1)==0,
        "current_base_exact_divergence":means.radial(b)+EPS*G.diff(1)==0,
        "both_current_velocity_moments_zero":(R**2*v0).integral()==(R*g0).integral()==0,
        "current_scalar_is_integrated_complete_angular_equation":M==EPS*J.diff(1),
        "projector_weighted_mass_one":(R**2*sigma).integral()==1,
        "projected_source_zero_weighted_mass":(R**2*QE).integral()==0,
        "projected_target_compact":Xi.at(0,1)==Xi.at(0,2)==0,
        "two_channels_reconstruct_current_equation":-means.radial(Xi,2)+sigma*M==E,
        "both_current_channels_nonzero":bool(Xi) and bool(M),
        "averaging_products_after_multiplication_required":(gamma*v).haar()!=gamma.haar()*v.haar(),
        "complete_covariance_scalar_flux_retained":bool((R**2*W[1][2]).integral()),
        "complete_covariance_radial_flux_retained":bool(W[0][1]),
    }
    scalar_inputs=[]
    for n in (1,3,5):
        # On this fixture I=[1,2], M2=(2-1)*max(1,4)=4.
        F=4*(pure_t_norm(G,n)*pure_t_norm(v0,n)+pure_t_norm(V,n)*pure_t_norm(g0,n)
             +pure_t_norm(gamma,n)*pure_t_norm(v,n)+pure_t_norm(W[1][2],n))
        checks[f"pure_T_current_flux_bound_n{n}"]=pure_t_norm(J,n)<=F
        checks[f"mixed_scalar_jet_split_n{n}"]=norm(J,n)<=pure_t_norm(J,n)+norm(M,n-1)/EPS
        checks[f"current_data_scalar_jet_bound_n{n}"]=norm(J,n)<=F+norm(M,n-1)/EPS
        scalar_inputs.append(dict(order=n,pure_T_upper=str(F),mixed_upper=str(F+norm(M,n-1)/EPS)))
        for i in range(1,n+1):
            for j in range(n+1-i):
                checks[f"exact_scalar_jet_n{n}_i{i}_j{j}"]=diff(diff(J,1,i),2,j)==diff(diff(M,1,i-1),2,j)/EPS
    return dict(checks=checks,scalar_inputs=scalar_inputs,J=J,M=M,sigma=sigma,
                scope="complete Haar current angular equation for manufactured moment-preserving fields; no actual graph-field realization or pulse supremum")


def pressure_routing_review(J,M):
    K=1+Z+Z**2+Z**3+Z**4+T+Z*T
    D=(R-1)**6*(2-R)**6*(1+Z+T)
    checks={};rows=[]
    for p in (3,4):
        routed=diff(K,1,p)*J
        for k in range(1,p+1):routed+=comb(p,k)*diff(K,1,p-k)*diff(M,1,k-1)/EPS
        checks[f"exact_product_scalar_moment_routing_p{p}"]=diff(K*J,1,p)==routed
        checks[f"whole_K_dependence_retained_p{p}"]=bool(diff(K,1,p)*J)
    viscous=EPS**5*(D*diff(K*J,1,3)).diff(1)
    timed=EPS**3*(D*(K*J).diff(1).diff(2)).diff(1)
    for m in (0,1):
        # Use A=<J> instead of its looser raw-field upper bound F.
        n=m+4;d=norm(D,m+1);k=norm(K,n)
        upper_v=qshift(m,4)*d*k*(EPS**5*pure_t_norm(J,n)+EPS**4*norm(M,n-1))
        ntime=m+3
        upper_t=qshift(m,3)*d*norm(K,ntime)*(EPS**3*pure_t_norm(J,ntime)+EPS**2*norm(M,ntime-1))
        checks[f"axial_viscous_current_data_budget_m{m}"]=norm(viscous,m)<=upper_v
        checks[f"slow_time_current_data_budget_m{m}"]=norm(timed,m)<=upper_t
        rows.append(dict(order=m,axial_viscous_upper=str(upper_v),slow_time_upper=str(upper_t)))
    # M=0 only removes Z derivatives of J, not the variable K product.
    constant_J=1+T+T**2
    checks["scalar_M_zero_does_not_erase_K_pressure_response"]=bool((D*diff(K*constant_J,1,3)).diff(1))
    return dict(checks=checks,rows=rows,scope="exact product identities and sufficient coefficient-norm pressure budgets")


def projection_control_review(sigma):
    J=(1+T)*(Z+Z**4)
    M=EPS*J.diff(1)
    checks={}
    for n in (1,4,9):
        E=n*sigma*M;moment=(R**2*E).integral()
        Xi=-primitive(E-sigma*moment,2)
        checks[f"projected_target_zero_with_scalar_N{n}"]=Xi==0 and moment==n*M and bool(moment)
    checks["moving_projection_bump_retained"]=bool(sigma.diff(1)) and (R**2*sigma.diff(1)).integral()==0
    return dict(checks=checks,scope="operator sources E=sigma*M, not a canonical field realization or cycle counterexample")


def fast_scalar_review():
    bump=(R-1)**6*(2-R)**6;psi=R*bump;f=means.radial(psi)
    a,d=Q(3,7),Q(2,5)
    A=(1+Z**4)*(1+T);B=(1+2*Z**4)*(1+T)
    carrier=Q(3)
    # Pairs are exact cosine/sine coefficients for the commuting graph
    # derivatives D_R=partial_R+carrier*partial_y and D_Z=EPS*partial_Z.
    zero=means.Poly()
    add=lambda x,y:tuple(p+q for p,q in zip(x,y))
    dz=lambda x:tuple(EPS*p.diff(1) for p in x)
    dr1=lambda x:(means.radial(x[0])+carrier*x[1],means.radial(x[1])-carrier*x[0])
    haar_linear=lambda x:Q(0)*x[0]+Q(0)*x[1]  # Exact circle means of cos and sin.
    haar_product=lambda x,y:(x[0]*y[0]+x[1]*y[1])/2
    def curl(amplitude):
        gamma=(f*amplitude,-carrier*psi*amplitude)
        return ((-EPS*psi*amplitude.diff(1),zero),tuple(p.shift(-1) for p in gamma),gamma)
    old=curl(a*A);increment=curl(d*B);new=tuple(add(x,y) for x,y in zip(old,increment))
    # The base-linear terms average to zero at both states; W is held fixed.
    J=lambda state:(R**2*haar_product(state[2],state[1])).integral()
    delta=J(new)-J(old)
    energy=(R*(f**2+carrier**2*psi**2)).integral().scalar()
    predicted=(a*d*A*B+d*d*B*B/2)*energy
    _,dv,dg=increment;_,v,gamma=old
    exact_change=(R**2*(haar_product(gamma,dv)+haar_product(dg,v)+haar_product(dg,dv))).integral()
    checks={
        "old_and_fast_Haar_means_all_zero":all(haar_linear(x)==0 for s in (old,increment) for x in s),
        "old_and_fast_graph_divergence_exact":all(all(p==0 for p in add(dr1(s[0]),dz(s[2]))) for s in (old,increment,new)),
        "old_fast_new_both_Haar_moments_zero":all((R**2*haar_linear(s[1])).integral()==(R*haar_linear(s[2])).integral()==0 for s in (old,increment,new)),
        "graph_sine_curl_component_retained":bool(old[2][1]) and bool(increment[2][1]),
        "fast_scalar_update_exact_with_all_three_products":delta==exact_change==predicted,
        "fast_correlation_energy_positive":energy>0,
        "zero_Haar_does_not_preserve_scalar_defect":bool(delta),
        "fourth_Z_scalar_jet_survives":bool(diff(delta,1,4).at(1,0)),
        "old_new_cross_cannot_be_dropped":delta!=(R**2*haar_product(dg,dv)).integral(),
    }
    for k in (4,9,15):
        eps=Q(1,2**k);old_amp=eps;fast_amp=eps**2
        scalar=old_amp*fast_amp+fast_amp**2/2
        checks[f"no_extra_epsilon_for_scalar_change_k{k}"]=scalar/eps**3==1+eps/2
        checks[f"moment_response_keeps_one_epsilon_k{k}"]=(eps*scalar)/eps**4==1+eps/2
    return dict(checks=checks,energy=str(energy),graph_carrier=str(carrier),scope="complete single-phase graph curl; finite endpoint polynomial fixture, not the prescribed fast inverse")


def review():
    current=current_fixture()
    parts={"current_scalar":{k:v for k,v in current.items() if k not in ("J","M","sigma")},
           "pressure_routing":pressure_routing_review(current["J"],current["M"]),
           "projection_control":projection_control_review(current["sigma"]),
           "fast_scalar":fast_scalar_review(),
           "complete_covariance":means.covariance_review(),"pressure_tails":means.pressure_tail_review()}
    return dict(checks={p+":"+k:v for p,r in parts.items() for k,v in r["checks"].items()},parts=parts)
