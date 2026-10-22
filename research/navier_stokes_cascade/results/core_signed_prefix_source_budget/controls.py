"""Full local mean equations, exponent witnesses and fresh-prefix controls."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import importlib.util
import sys
from budgets import P,FIXED,STATE,LINEAR,full_mean,source_recipe,fixture,haar,pole_ledger

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module;spec.loader.exec_module(module);return module


@lru_cache(None)
def modules():
    return (load("prefix_mean_exact",HERE.parent/"complete_covariance_initial_means/means.py"),
            load("prefix_actual_fast",HERE.parent/"core_fast_inverse_budget/proof_checks.py"),
            load("prefix_regular_exact",HERE.parent/"regular_source_foundation/source_checks.py"))


def full_equations():
    e,_,_=modules();R,Z,T,Y=e.R,e.Z,e.T,e.U;eps=e.EPS
    bump=(R-1)**6*(2-R)**6
    dr=lambda f:f.diff(0)+R**2*f.diff(3)
    dz=lambda f:eps*f.diff(1)
    dt=lambda f:-eps*f.diff(2)+Q(2,7)*f.diff(3)
    de=lambda f,j:dr(f)+j*f.shift(-1)
    lap=lambda f:dr(dr(f))+dr(f).shift(-1)+dz(dz(f))
    def from_potential(psi,v):return (-dz(psi),v,de(psi,1))
    base=from_potential(R**2*bump*(1+Z+T),R**2*(2+Z+T))
    mean=from_potential(R*bump*(1+Z+Y+Z*Y)*(1+T),R*bump*(1+Y+Y**2)*(1+Z))
    b,V,G=base;x,u,z=mean;total=tuple(a+c for a,c in zip(base,mean));bt,Vt,Gt=total
    _,_,_,W= e.fixture();W=[[v*(1+Y)**2 for v in row] for row in W]
    st=R**2*bump*(1+Y+Z);sz=R*bump*(1+Y+T);pm=R*bump*(1+Z+Z**2)*(1+Y)
    srr=2*b*x+x*x+W[0][0];srz=b*z+G*x+x*z+W[0][2];stt=2*V*u+u*u+W[1][1]
    gr=-dt(x)-de(srr,1)-dz(srz)+stt.shift(-1)+eps*(lap(x)-x.shift(-2))
    et=dt(u)+de(b*u+V*x+x*u+W[0][1]-st,2)+dz(G*u+V*z+z*u+W[1][2])-eps*(lap(u)-u.shift(-2))
    ez=dt(z)+de(b*z+G*x+x*z+W[0][2]-sz,1)+dz(2*G*z+z*z+W[2][2]+pm)-eps*lap(z)
    advr=lambda state:state[0]*dr(state[0])+state[2]*dz(state[0])-state[1]**2*R.shift(-2)
    # R.shift(-2)=R^-1, consistent with the exact polynomial convention.
    advt=lambda state:state[0]*dr(state[1])+state[2]*dz(state[1])+(state[0]*state[1]).shift(-1)
    advz=lambda state:state[0]*dr(state[2])+state[2]*dz(state[2])
    direct_gr=-dt(x)-(advr(total)-advr(base))-de(W[0][0],1)-dz(W[0][2])+W[1][1].shift(-1)+eps*(lap(x)-x.shift(-2))
    direct_et=dt(u)+advt(total)-advt(base)+de(W[0][1]-st,2)+dz(W[1][2])-eps*(lap(u)-u.shift(-2))
    direct_ez=dt(z)+advz(total)-advz(base)+de(W[0][2]-sz,1)+dz(W[2][2]+pm)-eps*lap(z)
    nofasttime=et-Q(2,7)*u.diff(3)
    flatlap=lambda f:f.diff(0).diff(0)+f.diff(0).shift(-1)+dz(dz(f))
    checks={
        "full_local_base_and_mean_divergence_exact":all(de(s[0],1)+dz(s[2])==0 for s in (base,mean,total)),
        "full_radial_conservative_and_advective_equations_match":gr==direct_gr,
        "full_angular_conservative_and_advective_equations_match":et==direct_et,
        "full_axial_conservative_and_advective_equations_match":ez==direct_ez,
        "pointwise_auxiliary_time_term_survives":et!=nofasttime,
        "full_graph_viscosity_cannot_be_replaced_by_Haar_laplacian":lap(u)!=flatlap(u),
        "graph_coefficient_derivative_survives_second_derivative":dr(dr(u))!=u.diff(0).diff(0)+2*R**2*u.diff(0).diff(3)+R**4*u.diff(3).diff(3),
        "radial_mean_and_all_six_covariances_nonzero":bool(x) and all(bool(W[i][j]) for i in range(3) for j in range(i,3)),
        "current_pressure_gradient_not_zero_or_deleted":bool(dz(pm)) and ez!=ez-dz(pm),
        "mean_cross_products_nonzero_before_any_Haar_average":bool(x*u) and bool(z*u) and bool(x*z),
    }
    return dict(checks=checks,scope="exact finite local polynomial differential algebra with a nonconstant commuting graph derivative; not a global torus or selected mean-map field")


def exponent_review():
    checks={};rows=[]
    for k in (Q(1,100000),Q(1,1000),Q(1,100)):
        mu=1-k;a=Q(1,2);alpha=2-3*k
        groups={"radial_mean_fast_time":mu+1,"radial_W_derivative":1-k,"swirl_mean_linear":mu,
            "radial_mean_viscosity":mu+2-2*k,"tangential_mean_fast_time":mu,
            "tangential_mean_viscosity":mu+1-2*k,"radial_W_minus_stress":1-k,
            "base_radial_mean_flux_derivative":mu+1-k,"axial_mean_flux":mu+1,
            "full_current_pressure_axial_flux":mu+1,"full_primary_linear":a+Q(1,2)-k,
            "full_signed_linear":alpha-k,"complete_wave_wave":2*a-k,"mean_wave_phase":a+mu-Q(1,2),
            "mean_wave_amplitude_derivative":a+mu-k}
        checks.update({f"k{k}_{name}_meets_mu":value>=mu for name,value in groups.items()})
        checks[f"k{k}_mean_wave_phase_loss_is_essential"]=a+mu-Q(1,2)<a+mu-k
        checks[f"k{k}_no_improved_full_row_assumption_needed"]=mu<2-3*k<2-2*k
        rows.append(dict(kappa=str(k),source_order=str(mu),limiting=[n for n,v in groups.items() if v==mu]))
    return dict(checks=checks,rows=rows,scope="rational exponent witnesses; all 0<kappa<1/4 inequalities derived in note, not actual rate lower bounds")


def budgets_review():
    checks={};b=fixture();trace=set();traced=fixture(trace)
    for n in (0,1,2):
        N=n+3;w=b.wave(N)["U"]
        data=dict(b.fixed(N),**b.initial(N),w=w*w);row=full_mean(n,data)
        for S in (1,3):
            rational=full_mean(n,{k:v.evaluate(S) for k,v in data.items()})
            checks[f"n{n}_S{S}_full_mean_positive_substitution"]=all(v.evaluate(S)==rational[k] for k,v in row.items() if isinstance(v,P))
            for key in ("Kp","c","graph","w"):
                dropped=dict(data);dropped[key]=P();bad=full_mean(n,dropped)
                checks[f"n{n}_S{S}_current_{key}_omission_detected"]=row["Ez"].evaluate(S)>bad["Ez"].evaluate(S)
        fixed=b.fixed(N);wave=b.wave(N);state=b.final(N)
        data={k:fixed[k] for k in ("b","G","graph","c","r")}
        source_data=dict(data,**b.linear(N),**wave,M=sum(state.values(),P()))
        result=source_recipe(n,source_data)
        for S in (1,2):
            direct=source_recipe(n,{k:v.evaluate(S) for k,v in source_data.items()})
            checks[f"n{n}_S{S}_source_recipe_exact_substitution"]=all(v.evaluate(S)==direct[k] for k,v in result.items() if isinstance(v,P))
        for key in ("Pi","R","M"):
            omitted=dict(source_data);omitted[key]=P()
            checks[f"n{n}_source_{key}_omission_detected"]=result["F"].evaluate(1)>source_recipe(n,omitted)["F"].evaluate(1)
        before=b.rows(n);after=b.rows(n,True)
        checks[f"n{n}_full_pressure_recomputed_after_fast"]=after["pm"].evaluate(1)>before["pm"].evaluate(1)
        checks[f"n{n}_fresh_slow_uses_postfast_scalars"]=all(after[k].evaluate(1)>before[k].evaluate(1) for k in ("P","Jtheta","Jz"))
        checks[f"n{n}_total_means_accumulate_both_maps"]=all(b.final(n)[k].evaluate(1)>b.postfast(n)[k].evaluate(1)>b.initial(n)[k].evaluate(1) for k in STATE)
    traced.source(0)
    highest=max(n for name,n in trace if name=="initial")
    checks["source_zero_uses_finitely_many_initial_field_jets"]=highest==15
    checks["source_does_not_request_later_complex_or_forced_states"]=set(name for name,n in trace)<=set(("fixed","initial","wave","linear","KT1","KA1","time_inverse","slow"))
    missingA=fixture();missingA.KA1=lambda n:P()
    checks["actual_fast_A1_changes_final_source_budget"]=b.source(0)["F"].evaluate(1)>missingA.source(0)["F"].evaluate(1)
    for p in (Q(0),Q(3,2),Q(4)):
        poles=pole_ledger(p)
        checks[f"p{p}_nonlinear_poles_enlarged_before_weighted_operator_calls"]=poles["pre_fast_rows"]==2*p and poles["fresh_post_fast_rows_and_moments"]==4*p
        checks[f"p{p}_closed_hull_slow_response_and_final_source_poles"]=poles["slow_response"]==0 and poles["source"]==3*p and poles["final_means"]==2*p
    _,fast,regular=modules()
    checks.update({"actual_fast:"+k:v for k,v in fast.mode_review()["checks"].items()})
    transport=regular.transport_checks()
    checks.update({"complete_transport:"+k:transport[k] for k in (
        "complete_harmonic_transport_matches","exact_complete_curl_longitudinal_identity","complete_amplitude_not_transverse","zero_and_negative_outputs_present")})
    linear=regular.linear_checks()
    checks["complete_linear_operator_exact_viscosity_pressure_connections"]=linear["principal_plus_remainder_matches"] and linear["all_omission_controls_nonzero"]
    return dict(checks=checks,highest_initial_field_jet_for_source_zero=highest,
                source_order_zero=b.source(0)["F"].record(),scope="explicit finite-prefix polynomial substitutions and inherited exact complete-curl/A1 controls; not actual parameter sampling")


@lru_cache(None)
def review():
    e,_,_=modules()
    parts=dict(full_current_equations=full_equations(),coarse_rates=exponent_review(),fresh_prefix_and_source=budgets_review(),
               complete_pressure_tails=e.pressure_tail_review())
    return dict(checks={name+":"+k:v for name,row in parts.items() for k,v in row["checks"].items()},parts=parts)
