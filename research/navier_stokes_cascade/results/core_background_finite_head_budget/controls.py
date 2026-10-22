"""Independent exact finite-head, cutoff-curl, chart and phase controls."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import factorial
from budgets import HERE,P,load,fixture,finite_head,increment_constant,chart_constants,length_compatibility


@lru_cache(None)
def helpers():
    return (load("finite_head_maps",HERE.parent/"modulation_repair_all_order/map_checks.py"),
            load("finite_head_geometry",HERE.parent/"core_reserved_geometry_recheck/proof_checks.py"),
            load("finite_head_selector",HERE.parent/"core_primary_tail_recheck/proof_checks.py"))


def tail_review():
    checks={};rows=[]
    for m in range(7):
        N=max(2,m);head={n:Q(n+m+1) for n in range(1,N)};out=finite_head(m,head)
        checks[f"m{m}_coefficient_one_always_in_finite_head"]=1 in head
        for t in (Q(1),Q(1,2),Q(1,17)):
            actual=Q(1,2**N)*t**(N-2)/(1-t/2)
            checks[f"m{m}_t{t}_entire_geometric_tail_bound"]=actual<=out["tail"]
            # Finite terms retain their original 2nh powers before normalization.
            total=sum((c*t**(2*n-2) for n,c in head.items()),Q(0))+actual
            checks[f"m{m}_t{t}_whole_interval_correction_constant"]=total<=out["correction"]
        rows.append(dict(order=m,start=N,head=str(out["head"]),tail=str(out["tail"])))
    checks["starting_diagonal_tail_at_one_loses_q_to_2h_rate"]=all(Q(1,2)/t/(1-t/2)>1 for t in (Q(1,4),Q(1,32)))
    return dict(checks=checks,rows=rows,scope="exact geometric series with t=q^h in (0,1]; analytic monotonic bounds cover the entire interval")


def cutoff_review():
    maps,_,_=helpers();J=maps.JetClass();checks={}
    q,X,eta=[J.variable(i,v,n=3,order=12) for i,v in enumerate((Q(1),Q(1,2),Q(1,3)))]
    chi=1-q*q/5+q**3/11;W=1+X+eta*eta;M=X*X+eta**3
    E=lambda z:q*z.derivative(0)
    h=Q(1,200);n=3;lam=2*n*h;rad=q**lam*(chi*W-E(chi)*M)
    bad=q**lam*chi*W
    checks["radial_cutoff_curl_commutator_nonzero"]=rad.value()!=bad.value()
    checks["commutator_sign_is_essential"]=rad.value()!=(q**lam*(chi*W+E(chi)*M)).value()
    for m in range(5):
        c=increment_constant(h,n,m,100,1+Q(2**(m+1),5)+Q(3**(m+1),11))
        z=rad
        rhs=chi*W-E(chi)*M
        for _ in range(m):z=E(z);rhs=lam*rhs+E(rhs)
        checks[f"m{m}_Euler_product_rule_retains_commutator"]=all((z-q**lam*rhs).partial(a)==0 for a in product(range(3),repeat=3) if sum(a)<=2)
        checks[f"m{m}_cutoff_extra_Euler_jet_bound"]=abs(z.value())<=c
    checks["Euler_second_jet_differs_from_frozen_second_jet"]=(E(E(q**lam))).value()==lam*lam and (q*q*(q**lam).derivative(0).derivative(0)).value()==lam*(lam-1)!=lam*lam
    return dict(checks=checks,scope="manufactured local polynomial cutoff jets; no claim that this polynomial is the actual compact cutoff or actual pressure coefficient")


def chart_review():
    maps,geometry,_=helpers();J=maps.JetClass();checks={};rows=[]
    h=Q(1,200);A=Q(1,2)+h;D=Q(1,2)-h
    for m in range(4):
        C=chart_constants(m,1,2)
        for eta0 in (Q(0),Q(3,5)):
            w,z,t=[J.variable(i,0,n=3,order=m) for i in range(3)]
            implicit=geometry.implicit_jet(h,eta0,m)
            s=J(1,n=3,order=m)+J({(0,i,j):v for (i,j),v in implicit.c.items()},n=3,order=m)
            R=1+w;X=R*R/(2*s);eta=(eta0+z)*s**(-D)
            normalized=s**Q(1,5)*(1+X+eta*eta)
            fields=dict(F=s**(-A)*normalized/R,G=s**(-A)*normalized,btilde=s**Q(-1,2)*normalized)
            # On X<=2, |eta|<=1, every Euler/X/eta jet through m
            # of q^(1/5)(1+X+eta^2), at q<=1, is at most four.
            for name,f in fields.items():
                actual=max(abs(f.partial(a)) for a in product(range(m+1),repeat=3) if sum(a)<=m)
                key=f"m{m}_eta{eta0}_{name}"
                checks[key+"_finite_chart_composition_bounded"]=actual<=4*C[name]
                rows.append(dict(order=m,eta=str(eta0),field=name,actual_max_raw_jet=str(actual)))
    checks["radial_reciprocal_is_not_dropped"]=chart_constants(3,Q(1,2),2)["F"]>chart_constants(3,1,2)["F"]
    return dict(checks=checks,rows=rows,scope="exact implicit coordinate jets and manufactured C^m profile jets through order three; analytic majorant, not sampled field norms, supplies general orders")


def comparison_review():
    b=fixture();out=b.comparison(Q(3,4),1,1,2,2,1);checks={};rows=[]
    C=out["Cphi"];pplus=out["carriers"]["theta"];qz=out["carriers"]["axial"]
    U=Q(3,4);Bs=Q(1,2);R0=Q(3,2);K=(Q(3,5),Q(4,5));N=(Q(4,5),Q(-3,5))
    for S in (Q(1),Q(3),Q(10)):
        L=2*S;eps=Q(1,16);k=S**3+3;R=R0+Q(1,2)/S**3
        for sign in (-1,1):
            pt=R0*Bs*(K[0]-sign*U*N[0]/L);pz=Bs*(K[1]-sign*U*N[1]/L);p=pt+1/k
            hr0=pt*N[0]/R0+pz*N[1]
            checks[f"S{S}_sign{sign}_frozen_shear_cancellation"]=hr0==-sign*Bs*U/L
            d=3*b.leading_bound(2)/S**3+eps*eps*b.bound(1)["error_FG"]
            fr=N[0]/R0+d/2;gr=N[1]+d/2
            for v in (Q(0),L/2,L):
                normal=(sign*Bs*U/2-v*(p*fr+pz*gr),p/R,pz-eps*v*(p+pz)*b.bound(1)["full_FG"]/2)
                ref=(Bs*sign*U*(Q(1,2)+v/L),Bs*K[0],Bs*K[1])
                discrepancy=sum(abs(x-y) for x,y in zip(normal,ref))
                bound=C*(1/S+1/S**2+S*(eps+eps*eps+1/k))
                key=f"S{S}_sign{sign}_v{v}"
                checks[key+"_full_normal_error_envelope"]=discrepancy<=bound
                rows.append(dict(S=str(S),sign=sign,v=str(v),L1_error=str(discrepancy)))
    checks["every_error_group_has_positive_cost"]=all(v>0 for v in out["groups"].values())
    checks["Cphi_sums_all_seven_groups"]=len(out["groups"])==7 and C==sum(out["groups"].values())
    trace=set();fixture(trace).comparison(1,1,1,1,2,1)
    checks["comparison_uses_only_coefficient_one_jets_through_one"]=max(t[2] for t in trace if t[0]=="profile")==1 and all(t[1]==1 for t in trace if t[0]=="profile")
    checks["comparison_uses_leading_two_and_cutoff_two"]=max(t[1] for t in trace if t[0]=="leading")==2 and max(t[1] for t in trace if t[0]=="cutoff")==2
    return dict(checks=checks,rows=rows,Cphi=str(C),groups={k:str(v) for k,v in out["groups"].items()},
                scope="manufactured admissible reference/carrier/error witnesses over the whole path endpoints and interior; continuum estimate follows seven analytic bounds, not these witnesses")


def selector_review():
    _,_,selector=helpers();p=selector.fixture();chosen=selector.select_tail(p)
    # Use the SAME selector inputs and its output. No new threshold is selected.
    C=Q(7);out=length_compatibility(1,2,C,p["radius"],p["length"],p["slow_cost"],p["slow_tolerance"])
    S=Q(chosen["ell_tail"]**2);beta=out["beta"]
    bad=length_compatibility(1,2,C,p["radius"],p["radius"],p["slow_cost"],p["slow_tolerance"])
    checks={"unchanged_selector_all_original_checks_pass":all(chosen["checks"].values()),
            "unchanged_length_budget_joins_parent_normal_margin":all(out["gates"].values()) and 36*C*C<=beta*S*S,
            "slow_gate_already_dominates_parent_rounding_gate":p["slow_tolerance"]/p["slow_cost"]<=1,
            "insufficient_original_length_is_reported_not_reselected":not bad["gates"]["original_length_budget_dominates_normal_margin"],
            "normal_floor_same_as_parent":out["inverse_normal_square"]==4/beta}
    return dict(checks=checks,ell_tail=str(chosen["ell_tail"]),scope="exact retained selector fixture, not substitution of actual finite background constants into the selected physical domain")


@lru_cache(None)
def review():
    parts=dict(tail=tail_review(),cutoff=cutoff_review(),chart=chart_review(),comparison=comparison_review(),selector=selector_review())
    return dict(checks={name+":"+k:v for name,p in parts.items() for k,v in p["checks"].items()},parts=parts)
