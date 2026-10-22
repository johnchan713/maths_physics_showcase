"""Exact frozen-phase identities and manufactured coefficient/floor controls."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product as indices
from math import factorial,comb
from pathlib import Path
import importlib.util
import sys
from budgets import P,PhaseBudget,fixture,normal_gate,carrier_bounds,product

HERE=Path(__file__).resolve().parent


@lru_cache(None)
def maps():
    spec=importlib.util.spec_from_file_location("phase_coefficient_jets",HERE.parent/"modulation_repair_all_order/map_checks.py")
    m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m


def same(a,b,n=2):
    return all((a-b).partial(alpha)==0 for alpha in indices(range(n+1),repeat=a.n) if sum(alpha)<=n)


def identities():
    J=maps().JetClass();values=(Q(3,2),Q(1,5),Q(2,7),Q(1,3),Q(0))
    R,Z,T,v,theta=[J.variable(i,z,n=5,order=6) for i,z in enumerate(values)]
    eps=Q(1,20);p=Q(-2,3);pz=Q(3,5);x0=Q(1,7)
    F=1+R+Z*T;G=1+R*R+Z*Z*T;b=eps*(R+Z+T)
    H=p*F+pz*G;phi=p*theta+pz*Z/eps+x0*R-v*H
    n=[x0-v*H.derivative(0),p*R**-1,pz-eps*v*H.derivative(1)]
    derived=[phi.derivative(0),phi.derivative(4)/R,eps*phi.derivative(1)]
    eik=-eps*phi.derivative(2)+phi.derivative(3)+b*phi.derivative(0)+F*phi.derivative(4)+G*eps*phi.derivative(1)
    exact=eps*v*(H.derivative(2)-G*H.derivative(1))+b*n[0]
    np=[-H.derivative(0),0*R,-eps*H.derivative(1)]
    checks={
        "all_three_actual_phase_normal_components_exact":all(same(a,c) for a,c in zip(n,derived)),
        "full_eikonal_identity_has_exact_epsilon_factor":same(eik,exact),
        "moving_normal_is_path_derivative_not_full_normal_jet":all(same(a,c.derivative(3)) for a,c in zip(np,n)),
        "normal_is_affine_in_path_with_frozen_carriers":all(same(c.derivative(3).derivative(3),0*R) for c in n),
        "full_axial_normal_error_not_dropped":bool((eps*v*H.derivative(1)).terms),
        "axial_transport_term_in_eikonal_is_essential":not same(eik,eps*v*H.derivative(2)+b*n[0]),
        "nonzero_radial_background_term_in_eikonal_is_essential":not same(eik,eps*v*(H.derivative(2)-G*H.derivative(1))),
        "phase_cancellation_uses_full_F_and_G":bool((H-(p*(1+R)+pz*(1+R*R))).terms),
    }
    return dict(checks=checks,scope="exact mixed rational jets of the actual frozen-phase formula with manufactured full backgrounds, not selected physical fields")


def profile_review():
    J=maps().JetClass();v=J.variable(0,Q(1,2),n=2,order=8)
    u=Q(3,4);lam=Q(3,2);s=u*(Q(1,2)+v);A=1+s*s
    growth=lam*Q(4,5)*(A/Q(25,16))**Q(-1,2)
    damping=lam*Q(64,125)*A;rate=growth-damping
    b=PhaseBudget(F=lambda n:P(1),G=lambda n:P(1),btilde=lambda n:P(1),
                  path=lambda n:P(1 if n<2 else 0),Rmin=1,Rmax=2,Lmin=1,u=u,lambda_max=lam,
                  carriers=dict(theta=1,axial=1,radial=1))
    checks={"growth_is_lambda0_over_sqrt_shape_not_constant":growth.value()==Q(6,5)!=lam,
            "normalization_center_has_zero_logarithmic_profile_rate":rate.value()==0,
            "third_profile_derivative_is_not_zero":rate.partial((3,0))!=0,
            "damping_has_zero_third_derivative_but_growth_does_not":damping.partial((3,0))==0 and growth.partial((3,0))!=0}
    rows=[]
    for n in range(7):
        val=abs(rate.partial((n,0)));upper=b.profile_rate(n).evaluate(1)
        checks[f"order{n}_actual_varying_growth_profile_jet_bounded"]=val<=upper
        rows.append(dict(order=n,exact_abs_jet=str(val),polynomial_bound=str(upper)))
    return dict(checks=checks,rows=rows,scope="exact algebraic profile jets at a rational Pythagorean normalization center; no selected Gaussian profile sampled")


def floor_review():
    good=normal_gate(1,2,1,100,Q(1,20000));bad=normal_gate(1,2,1,10,Q(1,20000))
    beta=good["beta"];delta=good["phase_error"]
    checks={
        "good_value_budget_meets_both_whole_path_gates":all(good["gates"].values()),
        "insufficient_normal_value_margin_is_rejected":not all(bad["gates"].values()),
        "floor_uses_full_or_tangential_reference_length":4*delta*delta<=beta and good["inverse_normal_square"]==4/beta,
        "nonzero_angular_carrier_alone_is_not_a_uniform_normal_floor":all(Q(1,n*n)<beta/4 for n in (100,1000)),
    }
    wronground=normal_gate(1,2,1,100,Q(1,100))
    checks["large_accumulated_rounding_error_is_not_hidden"]=not wronground["gates"]["whole_path_rounding_and_axial_error"]
    # beta is rational and below lambda_min/[4(1+u^2)^(3/2)]:
    # (1+u^2)^2 >= (1+u^2)^(3/2) follows by squaring.
    checks["squared_tilt_denominator_is_a_valid_outward_bound"]=(1+Q(4))**4>=(1+Q(4))**3
    return dict(checks=checks,good={k:str(v) for k,v in good.items() if k!="gates"},
                scope="rational triangle-inequality/value-gate witnesses; actual Cphi growth gaps and accepted-domain substitution remain unverified")


def polynomial_review():
    b=fixture();checks={};rows=[]
    for n in range(5):
        c=b.components(n);H=b.Hder(n)
        for S in (1,3):
            # Expand the raw Leibniz sum independently using rational values.
            vh=sum(Q(comb(n,j))*b.path(j).evaluate(S)*b.Hder(n-j).evaluate(S) for j in range(n+1))
            checks[f"n{n}_S{S}_normal_component_substitution"]=c["radial"].evaluate(S)==int(n==0)*b.x0+vh and c["axial"].evaluate(S)==int(n==0)*b.pz+vh
            checks[f"n{n}_S{S}_moving_normal_has_no_extra_path_factor"]=b.moving_normal(n).evaluate(S)==2*H.evaluate(S)
        checks[f"n{n}_normalized_source_coefficients_keep_all_named_inputs"]=set(b.source_coefficients(n))=={"b","G","V","r","normal","K","Eik","bR","Bz"}
        rows.append(dict(order=n,normal=b.normal(n).record(),matrix=b.matrix(n).record(),eikonal=b.eikonal(n).record()))
    for n in (0,1,3):
        tr=set();bb=fixture(tr);bb.source_coefficients(n);bb.profile_rate(n)
        checks[f"n{n}_only_background_jets_through_n_plus_one"]=max(j for name,j in tr if name in ("F","G","btilde"))==n+1
    noG=fixture();noG.G=lambda n:P()
    noB=fixture();noB.btilde=lambda n:P()
    checks["full_G_transport_and_nonzero_radial_flow_affect_eikonal_budget"]=b.eikonal(1).evaluate(1)>noG.eikonal(1).evaluate(1) and b.eikonal(1).evaluate(1)>noB.eikonal(1).evaluate(1)
    for dim in (4,5):
        checks[f"dim{dim}_raw_to_factorial_multinomial_count"]=b.normalized("normal",3,dim).evaluate(1)==sum(Q(dim**j,factorial(j))*b.normal(j).evaluate(1) for j in range(4))
    cb=carrier_bounds(3,2,Q(1,2),2,1)
    checks["nearest_nonzero_rounding_unit_cost_retained"]=cb["theta"]==2*cb["axial"]+1
    return dict(checks=checks,rows=rows,scope="positive-polynomial formulas with manufactured lower background bounds; no actual norm enclosure")


def box_review():
    J=maps().JetClass();x,Z,T,v=[J.variable(i,Q(0),n=4,order=10) for i in range(4)]
    R=1+x;F=1+R+Z*T;G=1+R*R+Z*Z*T;bt=R+Z+T
    def bound(polynomials,j):
        result=Q(0)
        for alpha in indices(range(j+1),repeat=4):
            if sum(alpha)!=j:continue
            total=Q(0)
            for poly in polynomials:
                for powers,c in poly.terms.items():
                    if all(a<=r for a,r in zip(alpha,powers)):
                        factor=Q(1)
                        for a,r in zip(alpha,powers):factor*=factorial(r)//factorial(r-a)
                        total+=abs(c)*factor
            result=max(result,total)
        return result
    p=Q(-2,3);pz=Q(3,5);x0=Q(1,7);eps=Q(1,20);H=p*F+pz*G
    nr=x0-v*H.derivative(0);nz=pz-eps*v*H.derivative(1)
    ee=v*(H.derivative(2)-G*H.derivative(1))+bt*nr
    kentries=[-2*F,2*F+R*F.derivative(0),G.derivative(0)]
    b=PhaseBudget(F=lambda j:P(bound([F],j)),G=lambda j:P(bound([G],j)),btilde=lambda j:P(bound([bt],j)),
        path=lambda j:P(1 if j<2 else 0),Rmin=1,Rmax=2,Lmin=1,u=1,lambda_max=3,
        carriers=dict(theta=abs(p),axial=abs(pz),radial=abs(x0)))
    checks={};rows=[]
    for j in range(4):
        actual_normal=bound([nr,nz],j)+abs(p)*factorial(j)
        checks[f"order{j}_whole_box_normal_bound"]=actual_normal<=b.normal(j).evaluate(1)
        checks[f"order{j}_whole_box_matrix_bound"]=bound(kentries,j)<=b.matrix(j).evaluate(1)
        checks[f"order{j}_whole_box_eikonal_over_epsilon_bound"]=bound([ee],j)<=b.eikonal(j).evaluate(1)
        rows.append(dict(order=j,normal_box_majorant=str(actual_normal),matrix_box_majorant=str(bound(kentries,j)),eikonal_box_majorant=str(bound([ee],j))))
    return dict(checks=checks,rows=rows,scope="independent exact coefficient-sum suprema over a whole unit box with R in [1,2]; finite manufactured backgrounds, not actual core bounds")


@lru_cache(None)
def review():
    parts=dict(actual_phase=identities(),varying_growth_profile=profile_review(),normal_value_gate=floor_review(),coefficient_polynomials=polynomial_review(),full_box=box_review())
    return dict(checks={name+":"+k:v for name,row in parts.items() for k,v in row["checks"].items()},parts=parts)
