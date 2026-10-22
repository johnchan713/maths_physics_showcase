"""Independent exact complete-current equations and scoped norm controls."""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import importlib.util
from budgets import P,KEYS,row_budget,weighted_covariance,fast_weighted,weighted_operator_budget,fixture,root,atlas,edge

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("haar_row_scalar",HERE.parent/"core_scalar_target_pressure_join/proof_checks.py")
scalar=importlib.util.module_from_spec(spec);spec.loader.exec_module(scalar)
means=scalar.means;R,Z,T,Y=means.R,means.Z,means.T,means.U;EPS=means.EPS


def pressure(gr,rho):
    P0=gr.integral();v=(gr-rho*P0).primitive()
    return P0,v-v.at(0,1)


def full_rows(base,state,W,stress,rho):
    b,V,G=base;B,v,g=state
    gr=means.averaged_gr(base,state,W[0][0],W[1][1],W[0][2]);P0,pm=pressure(gr,rho)
    et=(EPS*v.haar().diff(2)+means.radial((b*v+V*B+B*v+W[0][1]-stress[0]).haar(),2)
        +EPS*(G*v+V*g+g*v+W[1][2]).haar().diff(1)-EPS*(means.lap(v.haar())-v.haar().shift(-2)))
    ez=(EPS*g.haar().diff(2)+means.radial((b*g+G*B+B*g+W[0][2]-stress[1]).haar())
        +EPS*((2*G*g+g*g+W[2][2]).haar()+pm).diff(1)-EPS*means.lap(g.haar()))
    return gr,P0,pm,et,ez


def decomposed(base,fast,slow,W,stress,rho,drop_fast=False):
    b,V,G=base;bf,vf,gf=fast;B,u,g=slow
    avg=lambda p:means.Poly() if drop_fast else p.haar()
    s=2*b*B+B*B+avg(bf*bf)+W[0][0].haar()
    h=b*g+G*B+B*g+avg(bf*gf)+W[0][2].haar()
    t=2*V*u+u*u+avg(vf*vf)+W[1][1].haar()
    gr=EPS*B.diff(2)-means.radial(s)-EPS*h.diff(1)+t.shift(-1)+EPS*(means.lap(B)-B.shift(-2))
    p,pm=pressure(gr,rho)
    et=(EPS*u.diff(2)+means.radial(b*u+V*B+B*u+avg(bf*vf)+W[0][1].haar()-stress[0],2)
        +EPS*(G*u+V*g+g*u+avg(gf*vf)+W[1][2].haar()).diff(1)-EPS*(means.lap(u)-u.shift(-2)))
    ez=(EPS*g.diff(2)+means.radial(b*g+G*B+B*g+avg(bf*gf)+W[0][2].haar()-stress[1])
        +EPS*(2*G*g+g*g+avg(gf*gf)+W[2][2].haar()+pm).diff(1)-EPS*means.lap(g))
    return gr,p,pm,et,ez


def current_review():
    bump=(R-1)**6*(2-R)**6
    psi0=R*bump*(1+Z+T)
    base=(-EPS*psi0.diff(1),R**2*(2+Z+T),means.radial(psi0))
    fp=R*bump*(1+Z+Z**2)*(1+T)*Y
    fast=(-EPS*fp.diff(1),R*bump*(1+Z+T)*Y,means.radial(fp))
    sp=R*bump*(1+2*Z)*(1+T)
    slow=(-EPS*sp.diff(1),(R**4*bump).diff(0).shift(-2)*(1+Z+T),means.radial(sp))
    state=tuple(a+b for a,b in zip(fast,slow));_,_,_,W=means.fixture()
    stress=(R**2*bump*(1+Z),R**2*bump*(1+T))
    rho=bump/bump.integral().scalar()*(1+Z*(R-Q(3,2)))
    actual=full_rows(base,state,W,stress,rho);split=decomposed(base,fast,slow,W,stress,rho)
    nofast=decomposed(base,fast,slow,W,stress,rho,True)
    gr,P0,pm,et,ez=actual;B,u,g=state
    Jt=(R**2*(base[2]*u+base[1]*g+g*u+W[1][2]).haar()).integral()
    Jz=(R*(2*base[2]*g+g*g+W[2][2]).haar()).integral()-(R**2*gr).integral()/2
    cr=(R**2*rho).integral()/2
    checks={
        "actual_fast_all_Haar_means_zero":all(v.haar()==0 for v in fast),
        "actual_slow_Haar_means_nonzero":all(bool(v.haar()) for v in slow),
        "fast_and_slow_divergence_exact":all(means.radial(v[0])+EPS*v[2].diff(1)==0 for v in (fast,slow,state)),
        "both_velocity_moments_zero":(R**2*u).haar().integral()==(R*g).haar().integral()==0,
        "full_radial_equation_matches_split":actual[0]==split[0],
        "full_pressure_defect_matches_split":actual[1]==split[1],
        "full_pressure_profile_matches_split":actual[2]==split[2],
        "full_current_angular_row_matches_split":actual[3]==split[3],
        "full_current_axial_row_matches_split":actual[4]==split[4],
        "fast_linear_time_and_viscosity_Haar_zero":all((EPS*v.diff(2)-EPS*means.lap(v)).haar()==0 for v in fast),
        "fast_slow_crosses_zero_only_after_averaging":all((f*s).haar()==0 for f in fast for s in slow) and bool(fast[1]*slow[1]),
        "all_six_complete_covariances_nonzero":all(bool(W[i][j]) for i in range(3) for j in range(i,3)),
        "omitted_fast_pressure_squares_detected":actual[2]!=nofast[2],
        "omitted_fast_angular_correlations_detected":actual[3]!=nofast[3],
        "omitted_fast_axial_correlations_detected":actual[4]!=nofast[4],
        "current_pressure_compact_both_ends":pm.at(0,1)==pm.at(0,2)==0,
        "angular_scalar_identity_retained":(R**2*et).integral()==EPS*Jt.diff(1),
        "whole_axial_pressure_scalar_identity_retained":(R*ez).integral()==EPS*(Jz+cr*P0).diff(1),
        "moving_pressure_factor_cannot_be_frozen":bool(cr.diff(1)*P0),
        "base_radial_velocity_not_zero":bool(base[0]),
    }
    badstate=tuple(v+s for v,s in zip(fast,(means.Poly(),R*bump,means.Poly())))
    checks["nonzero_old_Haar_breaks_zero_Haar_split"]=full_rows(base,tuple(a+b for a,b in zip(badstate,slow)),W,stress,rho)!=decomposed(base,badstate,slow,W,stress,rho)
    return dict(checks=checks,scope="exact finite-regularity complete-current polynomial fixtures; fast coefficients are not the actual inverse values")


def exponent_review():
    checks={};rows=[]
    for k in (Q(1,100000),Q(1,1000),Q(1,100)):
        H=1-k;beta=2-2*k;alpha=2-3*k
        groups={"selected_real_curl_square":beta,"slow_time":2,"base_radial_slow":2,
                "base_swirl_radial":2,"slow_radial_tangential":3,"fast_radial_tangential":2*H+1,
                "base_axial_slow":2,"slow_axial_square_flux":3,"fast_axial_square_flux":2*H+1,
                "principal_variance_axial_flux":2,"whole_pressure_axial_flux":2,"slow_viscosity":2}
        for name,a in groups.items():checks[f"k{k}_{name}_fits_sharper_Haar_rate"]=a>=beta
        checks[f"k{k}_unaveraged_fast_graph_viscosity_not_promoted"]=1+H-2*k==alpha<beta
        checks[f"k{k}_pressure_itself_has_order_one"]=1<beta
        checks[f"k{k}_full_mean_pressure_fast_square_is_higher"]=2*H>1
        rows.append(dict(kappa=str(k),Haar_order=str(beta),full_row_order=str(alpha),limiting=[name for name,a in groups.items() if a==beta]))
    return dict(checks=checks,rows=rows,scope="general algebraic inequalities for 0<kappa<1/2 proved in the note; rational fixtures are not an actual field rate lower bound")


def norm_review():
    checks={};data=fixture()
    for n in (0,1,3):
        out=row_budget(n,data)
        for S in (1,4,16):
            direct=row_budget(n,{k:v.evaluate(S) for k,v in data.items()})
            checks[f"n{n}_S{S}_full_positive_row_substitution"]=all(v.evaluate(S)==direct[k] for k,v in out.items() if isinstance(v,P))
            nofast=dict(data,Ft=P(),Fr=P());lower=row_budget(n,nofast)
            checks[f"n{n}_S{S}_fast_correlations_in_both_row_budgets"]=out["Etheta"].evaluate(S)>lower["Etheta"].evaluate(S) and out["Ez"].evaluate(S)>lower["Ez"].evaluate(S)
            no_pressure=dict(data,Kp=P());reduced=row_budget(n,no_pressure)
            checks[f"n{n}_S{S}_whole_pressure_not_removed_from_axial_budget"]=out["Ez"].evaluate(S)>reduced["Ez"].evaluate(S)
        checks[f"n{n}_pressure_jet_input_is_n_plus_one_not_N_plus_two"]=out["field_order"]==n+3 and out["pressure_input_order"]==n+1
    for n in range(5):
        e=lambda k:P(k+1);kt=lambda k:P(k+2);ka=lambda k:P(k+3)
        f=fast_weighted(n,e,e,kt,ka,P({1:1}))
        checks[f"n{n}_five_variable_raw_to_normalized_count"]=comb(n+5,5)>=atlas.seed.fast.slow_count(n)
        checks[f"n{n}_actual_axial_A1_budget_retained"]=f["Ft"].evaluate(1)>fast_weighted(n,e,e,kt,lambda k:P(),P({1:1}))["Ft"].evaluate(1)
        checks[f"n{n}_weighted_fast_radial_extra_input_jet"]=f["Fr"].evaluate(1)>0
    primary=root.fixture()
    for n in (0,1,2):
        w=weighted_covariance(n,primary)
        checks[f"n{n}_real_selected_curl_square_is_not_full_covariance"]=w["W"].evaluate(1)>w["Dcurl"].evaluate(1)>0
        # The weighted bound does not depend on an unweighting constant.
        old=primary.width;primary.width=Q(1,1000)
        other=weighted_covariance(n,primary);primary.width=old
        checks[f"n{n}_weighted_primary_norms_not_unweighted"]=w["W"].record()==other["W"].record()
    # Feed the derived row into the immediate parent's actual target routing.
    for n in (0,1):
        target=edge.fixture();target.row=lambda k:row_budget(k,fixture())["Etheta"]+row_budget(k,fixture())["Ez"]
        value=root.fixture().signed(n,target.V)
        checks[f"n{n}_current_rows_feed_two_edge_and_original_root_chain"]=value.evaluate(1)>0
    return dict(checks=checks,scope="positive-polynomial norm substitution with manufactured finite inputs; no selected parameter evaluation")


def operator_review():
    checks={};bump=(R-1)**6*(2-R)**6;chi=3*(R-1)**2-2*(R-1)**3
    for e in (0,1,2):
        mass=(R**e*bump).integral();left=scalar.primitive(bump,e)
        cut=left-(chi*mass).shift(-e)
        checks[f"e{e}_nonzero_moment_left_primitive_has_wrong_outer_end"]=bool(left.at(0,2)) and bool(mass)
        checks[f"e{e}_edge_directed_cutoff_primitive_is_compact"]=cut.at(0,1)==cut.at(0,2)==0
        checks[f"e{e}_interior_compact_error_retained"]=means.radial(cut,e)==bump-(chi.diff(0)*mass).shift(-e)
    for p in (0,3,8):
        w=weighted_operator_budget(p,4,Q(1,2),P({0:2,1:1}),P({0:3,1:2}),Q(1))
        narrow=weighted_operator_budget(p,4,Q(1,2),P({0:2,1:1}),P({0:3,1:2}),Q(1,2))
        checks[f"p{p}_weighted_compact_constants_retain_profile_and_gap"]=narrow["KT"].evaluate(1)>w["KT"].evaluate(1)>0 and narrow["KA"].evaluate(1)>w["KA"].evaluate(1)
        checks[f"p{p}_primitive_collar_must_fit_actual_cutoff_plateau"]=narrow["proof_collar"]==Q(1,2) and w["proof_collar"]==1
        checks[f"p{p}_weighted_A1_not_assigned_unweighted_constant"]=w["KA"].evaluate(1)>5
    return dict(checks=checks,scope="exact one-sided cutoff algebra and moderate fixed-profile constant recipes; not arbitrary-operator weight preservation")


def review():
    parts=dict(current_complete_equations=current_review(),sharper_Haar_rate=exponent_review(),
               weighted_norm_chain=norm_review(),weighted_fixed_maps=operator_review(),complete_pressure_tails=means.pressure_tail_review())
    return dict(checks={name+":"+k:v for name,row in parts.items() for k,v in row["checks"].items()},parts=parts)
