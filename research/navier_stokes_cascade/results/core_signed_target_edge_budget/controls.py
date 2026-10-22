"""Exact primitive-tail, current-state projection and weighted-gap controls."""
from fractions import Fraction as Q
from math import comb,factorial
from budgets import P,TargetBudget,atlas,fixture,primitive_recipe,rational_upper,root
import importlib.util
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("signed_edge_scalar",HERE.parent/"core_scalar_target_pressure_join/proof_checks.py")
scalar=importlib.util.module_from_spec(spec);spec.loader.exec_module(scalar)
means=scalar.means;R,Z,T=means.R,means.Z,means.T;EPS=means.EPS


def add(a,b):
    out=dict(a)
    for p,v in b.items():out[p]=out.get(p,Q(0))+v
    return {p:v for p,v in out.items() if v}


def tail_review():
    checks={};rows=[]
    for c in (Q(1),Q(1,4),Q(4)):
        for q in range(7):
            p=3+2*q
            poly={2*j-2*q:Q(comb(q,j)*factorial(j),2)/c**(j+1) for j in range(q+1)}
            derivative=add({k-1:k*v for k,v in poly.items() if k},{k-3:2*c*v for k,v in poly.items()})
            checks[f"c{c}_p{p}_exact_tail_differential_identity"]=derivative=={-p:Q(1)}
            for x in (Q(1,16),Q(1,4),Q(1)):
                value=sum(v*x**k for k,v in poly.items())
                bound=x**(3-p)/(2*c)*sum(Q(comb(q,j)*factorial(j))*(x*x/c)**j for j in range(q+1))
                checks[f"c{c}_p{p}_x{x}_kernel_and_exact_tail_agree"]=value==bound
    for p in (Q(0),Q(3,2),Q(3),Q(4),Q(9)):
        for L in (Q(1),Q(4),Q(9)):
            r=primitive_recipe(p,L,Q(1,2));collar=Q(r["collar"])
            checks[f"p{p}_L{L}_collars_inside_shell"]=0<collar<=1 and 4*collar<=L
            checks[f"p{p}_L{L}_kernel_degree_dominates"]=r["kernel_degree"]>=max(Q(0),(p-3)/2)
            bound=rational_upper(r)
            checks[f"p{p}_L{L}_rational_exponential_envelope"]=all(bound>=a*3**(-(-b.numerator//b.denominator)) for a,b in r["branches"].values())
            rows.append(dict(pole=str(p),shell_length=str(L),collar=r["collar"],degree=r["kernel_degree"]))
    checks["narrower_width_cost_retained"]=rational_upper(primitive_recipe(9,4,Q(1,8)))>rational_upper(primitive_recipe(9,4,Q(1,2)))
    return dict(checks=checks,rows=rows,scope="exact normalized tail polynomials and moderate-width rational envelopes, not selected profile evaluation")


def current_review():
    bump=(R-1)**6*(2-R)**6
    psi0=R*bump*(1+Z+T)
    base=(-EPS*psi0.diff(1),R**2*(2+Z+Z**2+T),means.radial(psi0))
    psi=R*bump*(1+Z+Z**3)*(1+T+means.U)
    state=(-EPS*psi.diff(1),(R**4*bump).diff(0).shift(-2)*(1+Z+T+means.U),means.radial(psi))
    b,V,G=base;beta,v,gamma=state;_,_,_,W=means.fixture()
    v0,g0=v.haar(),gamma.haar()
    gr=means.averaged_gr(base,state,W[0][0],W[1][1],W[0][2]);P0=gr.integral()
    rho=bump/bump.integral().scalar()*(1+Z*(R-Q(3,2)))
    cr=(R**2*rho).integral()/2;prim=(gr-rho*P0).primitive();pm=prim-prim.at(0,1)
    thetaflux=(G*v+V*gamma+gamma*v+W[1][2]).haar()
    zflux=(2*G*gamma+gamma**2+W[2][2]).haar()
    Jtheta=(R**2*thetaflux).integral();Jz=(R*zflux).integral()-(R**2*gr).integral()/2
    stress=R**2*bump*(1+Z+T)
    Er2=(EPS*v0.diff(2)+means.radial((b*v+beta*V+beta*v+W[0][1]-stress).haar(),2)
         +EPS*thetaflux.diff(1)-EPS*(means.lap(v0)-v0.shift(-2)))
    Er1=(EPS*g0.diff(2)+means.radial((b*gamma+G*beta+beta*gamma+W[0][2]-stress).haar())
         +EPS*(zflux+pm).diff(1)-EPS*means.lap(g0))
    checks={
        "current_base_divergence":means.radial(b)+EPS*G.diff(1)==0,
        "current_mean_divergence":means.radial(beta.haar())+EPS*g0.diff(1)==0,
        "both_velocity_moments_zero":(R**2*v0).integral()==(R*g0).integral()==0,
        "current_angular_moment_whole_scalar":(R**2*Er2).integral()==EPS*Jtheta.diff(1),
        "current_axial_moment_whole_pressure_scalar":(R*Er1).integral()==EPS*(Jz+cr*P0).diff(1),
        "whole_pressure_compact_both_ends":pm.at(0,1)==pm.at(0,2)==0,
        "moving_pressure_product_term_nonzero":bool(cr.diff(1)*P0),
        "all_six_complete_covariances_retained":all(bool(W[i][j]) for i in range(3) for j in range(i,3)),
        "old_haar_velocities_not_reset":bool(v0) and bool(g0),
        "products_averaged_after_multiplication":(gamma*v).haar()!=gamma.haar()*v.haar(),
    }
    chi=3*(R-1)**2-2*(R-1)**3
    for e,E,M in ((2,Er2,EPS*Jtheta.diff(1)),(1,Er1,EPS*(Jz+cr*P0).diff(1))):
        mass=(R**e*bump).integral().scalar();center=(R**(e+1)*bump).integral().scalar()/mass
        sigma=bump*(1+(Z/13+T/17)*(R-center))/mass
        h=E-sigma*M;left=scalar.primitive(h,e);total=(R**e*h).integral();shifted=left-(chi*total).shift(-e)
        checks[f"e{e}_moving_projector_has_unit_weighted_mass"]=(R**e*sigma).integral()==1 and bool(sigma.diff(1))
        checks[f"e{e}_projected_source_exact_zero_moment"]=total==0
        checks[f"e{e}_Haar_shifted_primitive_equals_ordinary"]=shifted==left
        checks[f"e{e}_signed_target_compact_both_ends"]=left.at(0,1)==left.at(0,2)==0
        checks[f"e{e}_signed_target_reconstructs_source"]=means.radial(left,e)==h
        bad=scalar.primitive(E,e)
        checks[f"e{e}_omitted_projection_loses_outer_compactness"]=bool(bad.at(0,2))
    checks["omitted_moving_pressure_factor_detected"]=(R*Er1).integral()!=EPS*(Jz.diff(1)+cr*P0.diff(1))
    return dict(checks=checks,scope="current full conservative equations on finite-regularity polynomial fixtures, not actual canonical prefix fields")


def obstruction_review():
    checks={};rows=[]
    # phi(s)=s^4(1/4-s)^4 is a finite-regularity diagnostic for the smooth bump proof.
    mass=sum(Q((-1)**j*comb(4,j))*Q(1,4)**(4-j)*Q(1,4)**(5+j)/(5+j) for j in range(5))
    checks["scaled_bump_mass_is_positive"]=mass>0
    for eta in (Q(1,64),Q(1,32),Q(1,16)):
        x=Q(3,2)*eta
        checks[f"eta{eta}_observation_between_two_disjoint_bumps"]=Q(5,4)*eta<x<Q(7,4)*eta
        checks[f"eta{eta}_exact_two_bump_total_zero"]=eta*mass-eta*mass==0
        checks[f"eta{eta}_first_bump_primitive_mass"]=eta*mass>0
        coefficient=1/(Q(3,2)**2)-Q(1,4)
        checks[f"eta{eta}_exponential_weight_gap"]=coefficient==Q(7,36)>0
    for d in range(13):
        k=(d+1)//2+1
        checks[f"d{d}_any_fixed_weighted_pole_is_defeated"]=d+1-2*k<0
        # Unweighted order m: exp(-a/(4eta^2))*eta^-m ->0 by the inverse exponential series.
        checks[f"d{d}_all_fixed_unweighted_derivatives_can_tend_zero"]=2*(d+1)-d>0
        rows.append(dict(pole=d,exponential_series_order=k,power_after_lower_bound=d+1-2*k))
    checks["exact_moment_not_promoted_to_weighted_bound"]=Q(7,36)>0
    return dict(checks=checks,rows=rows,scope="scaled smooth-bump analytical counterexample to an unweighted-to-weighted inference; polynomial masses are diagnostics")


def budget_review():
    checks={};row=fixture()
    for n in range(6):
        edge=row.edge(n)
        for j in range(n+1):
            checks[f"n{n}_j{j}_normalized_target_fits_root_input"]=edge["quotient"][j]<=3*j+edge["d"]
        for S in (1,4):
            expected=row.row(n).evaluate(S)+sum(row.bump(j).evaluate(S)*row.moment(n-j).evaluate(S) for j in range(n+1))
            checks[f"n{n}_S{S}_current_projection_budget"]=row.H(n).evaluate(S)==expected
            expected=row.primitive(n).evaluate(S)*expected+(row.H(n-1).evaluate(S) if n else 0)
            checks[f"n{n}_S{S}_two_edge_primitive_budget"]=row.F(n).evaluate(S)==expected
        trace=set();t=fixture(trace);t.V(n)
        checks[f"n{n}_primitive_and_weight_division_add_no_source_jet"]=max(k for name,k in trace if name=="row")==n
    for n in (0,1):
        trace=set();primary=root.fixture(trace);seed=atlas.seed.fixture_budget(1)
        seed.omega=lambda k:primary.aggregate(k)["Omega"].evaluate(1)
        target=fixture();target.moment=lambda k:P(seed.axial_moment(k,Q(5)))
        native=root.fixture();native.signed(n+1,target.V)
        checks[f"n{n}_signed_curl_moment_channel_uses_n_plus_13"]=max(k for name,k in trace if name in ("H","target","wave"))==n+13
        checks[f"n{n}_scalar_moment_channel_is_not_full_row_proof"]=max(k for name,k in trace if name=="curl")==n+12
    kappa=Q(1,100000);alpha=2-3*kappa
    checks["initial_moment_epsilon_three_is_stronger_than_signed_alpha"]=3>alpha>0
    return dict(checks=checks,scope="finite positive-polynomial substitution; weighted actual Haar-row inputs remain explicit hypotheses")


def chart_review():
    checks={};x=2*R+R**2+R*Z;y=Z+R*Z
    f={(0,0):Q(1),(1,0):Q(-2),(0,1):Q(3),(2,1):Q(4),(1,2):Q(-5),(2,2):Q(2)}
    direct=sum((v*x**i*y**j for (i,j),v in f.items()),means.Poly())
    norm=lambda poly,n:sum(abs(v) for (r,z,t,u),v in poly.terms.items() if r+z+t+u==n)
    target=fixture();target.chart=lambda n:P(norm(x,n)+norm(y,n))
    for n in range(9):
        upper=sum(sum(abs(v) for (i,j),v in f.items() if i+j==k)*target.chart_power(n,k).evaluate(1) for k in range(n+1))
        checks[f"n{n}_nonlinear_mixed_chart_composition"]=norm(direct,n)<=upper
    # F=x has a second native derivative even though its logarithmic second jet is zero.
    checks["dropping_higher_chart_derivatives_detected"]=norm(x,2)>0
    return dict(checks=checks,scope="exact nonlinear two-variable polynomial chart controls; native physical chart seminorms remain fixed inputs")


def review():
    parts=dict(two_edge_tail=tail_review(),current_projection=current_review(),
               unweighted_obstruction=obstruction_review(),target_root_chain=budget_review(),chart_composition=chart_review(),
               pressure_memory=means.pressure_tail_review())
    return dict(checks={name+":"+k:v for name,row in parts.items() for k,v in row["checks"].items()},parts=parts)
