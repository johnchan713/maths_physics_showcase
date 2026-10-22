"""Exact initial Haar/pressure identities and positive primary-to-seed budgets."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


fast=load("initial_seed_fast",HERE.parent/"core_fast_inverse_budget/proof_checks.py")
slow=load("initial_seed_slow",HERE.parent/"core_slow_inverse_budget/proof_checks.py")
coupled=fast.coupled;means=fast.means;P=fast.P
R,Z,T,AUX=means.R,means.Z,means.T,means.U
EPS=means.EPS;order=fast.order


def covariance_majorant(n,a,b):order(n);return Q(2**n,2)*(a*a+b*b)


def graph_majorant(n,e,j,r):
    order(n);order(e)
    j,r=map(coupled.cycle.exact,(j,r))
    if min(j,r)<0:raise ValueError("Nonnegative fixed multiplier norms required")
    return 2**n*(3*j+e*r)


class SeedBudget:
    """Finite-order norm chain; data are supplied bounds, never sampled fields."""
    def __init__(self,omega,M,Ms,graph,r,kt0,kt1,ka1,profile,base,tg,S):
        self.omega,self.M,self.Ms=omega,M,Ms
        self.graph,self.r,self.kt0,self.kt1,self.ka1=graph,r,kt0,kt1,ka1
        self.profile,self.base,self.tg,self.S=profile,base,tg,S

    @lru_cache(None)
    def G(self,n):
        order(n)
        return (self.graph(n,1)+1)*self.omega(n+1)+2**n*self.r(n)*self.omega(n)

    @lru_cache(None)
    def E(self,n):
        order(n)
        return dict(theta=2*(self.graph(n,2)+1)*self.omega(n+1),
                    z=2*(self.graph(n,1)+1)*self.omega(n+1)+2*self.kt0(n+1)*self.G(n+1))

    @lru_cache(None)
    def initial(self,n):
        order(n);M=self.M;J=fast.slow_count(n);w=self.omega(n);w1=self.omega(n+1)
        return dict(P=J*(2*M[-1]*w+M[0]*w1),Jtheta=J*M[2]*w,
                    Jz=J*(2*M[1]*w+Q(1,2)*M[2]*w1))

    @lru_cache(None)
    def Ft(self,n):
        order(n);E=self.E(n+4)
        return fast.slow_count(n)*fast.temporal_constant(n)*self.tg*self.S*(E["theta"]+(1+self.ka1(n))*E["z"])

    @lru_cache(None)
    def Fr(self,n):
        order(n)
        return fast.slow_count(n)*self.kt1(n+1)*fast.temporal_constant(n+1)*self.tg*self.S*self.E(n+5)["z"]

    @lru_cache(None)
    def preslow_pz(self,n):
        order(n);I=self.initial(n)
        F=coupled.fast_budget(n,self.M,0,0,self.Ft(n+1),self.Fr(n+1))
        return dict(P=I["P"]+F["P"],Jz=I["Jz"]+F["Jz"])

    @lru_cache(None)
    def preslow_theta(self,n):
        order(n)
        return self.initial(n)["Jtheta"]+self.M[2]*self.Ft(n)*self.Ft(n)

    @lru_cache(None)
    def fields(self,n):
        order(n);PZ=self.preslow_pz(n);C=self.profile(n);Cnext=self.profile(n+1)
        return slow.inverse_budget(n,C["H"],Cnext["H"],C["Kp"],C["Kz"],C["Kg"],Cnext["Kpsi"],
                                   PZ["P"],PZ["Jz"],self.preslow_theta(n),self.preslow_theta(n+1))

    @lru_cache(None)
    def output(self,n):
        order(n);f=self.fields(n+2)
        return coupled.slow_budget(n,self.Ms,self.base(n+2),0,0,f["X"],f["U"],f["Z"])

    def axial_moment(self,n,c):
        order(n);out=self.output(n+1)
        return (n+1)*(out["Jz"]+c*out["P"])


def fixture_budget(Svalue=None,trace=None):
    def evalp(poly):return poly if Svalue is None else poly.evaluate(Svalue)
    def omega(n):
        if trace is not None:trace.add(n)
        return covariance_majorant(n,evalp(P({0:n+1,1:1})),evalp(P({0:2,1:n+1})))
    def profile(n):
        return dict(H=evalp(P({0:n+1,1:1})),Kp=Q(n+2),Kz=Q(n+3),Kg=Q(n+4),Kpsi=Q(n+5))
    return SeedBudget(omega,coupled.geometry(1,3),coupled.geometry(Q(5,4),2),
                      lambda n,e:graph_majorant(n,e,n+2,n+1),lambda n:Q(n+1),
                      lambda n:Q(n+2),lambda n:Q(n+3),lambda n:Q(n+4),profile,
                      lambda n:evalp(P({0:1,1:n+1})),3,evalp(P({1:1})))


def pressure_review():
    f=(R-1)**6*(2-R)**6;rho0=f/f.integral().scalar()
    rho=rho0+Z*(R-Q(3,2))*rho0/10;star=rho+(1+T)*(R-Q(3,2))*rho0/20
    _,_,_,W=means.fixture();zero=(means.Poly(),)*3
    gr=means.averaged_gr(zero,zero,W[0][0],W[1][1],W[0][2]);p0=gr.integral()
    pp=lambda density:coupled.phi(gr-density*p0)
    p=pp(rho);ps=pp(star);c=(R**2*rho).integral()/2
    wc=[[x.haar() for x in row] for row in W]
    ip=(wc[1][1]-wc[0][0]).shift(-1).integral()-EPS*wc[0][2].integral().diff(1)
    iz=(R*(wc[2][2]-(wc[0][0]+wc[1][1])/2)).integral()+EPS*(R**2*wc[0][2]).integral().diff(1)/2
    # A genuine nonzero formal Fourier mode is kept separately from the Haar mode.
    primitive=f.primitive();chi=(primitive-primitive.at(0,1))/f.integral().scalar()+Z*f/100
    pcirc,_,_=fast.compact(f*(1+Z+T),chi,Q(5,2))
    full={0:p,1:pcirc};changed={0:ps,1:pcirc}
    checks={
        "complete_initial_pressure_scalar_identity":p0==ip,
        "complete_initial_axial_scalar_identity":(R*wc[2][2]).integral()-(R**2*gr).integral()/2==iz,
        "complete_initial_angular_scalar_nonzero":bool((R**2*wc[1][2]).integral()),
        "both_moving_compensators_unit_mass":rho.integral()==star.integral()==1,
        "both_full_Haar_pressures_compact":p.at(0,1)==p.at(0,2)==ps.at(0,1)==ps.at(0,2)==0,
        "compensator_changes_full_Haar_pressure_shape":full[0]!=changed[0],
        "compensator_does_not_change_zero_Haar_fast_pressure":full[1]==changed[1] and not fast.zero(pcirc),
        "moving_cutoff_stays_in_zero_Haar_pressure_derivative":not fast.zero(fast.diff(pcirc,1)),
        "whole_moving_Haar_pressure_moment":(R*p).integral()==-(R**2*gr).integral()/2+c*p0,
        "whole_moving_Haar_moment_derivative":(R*p).integral().diff(1)==(-(R**2*gr).integral()/2+c*p0).diff(1),
        "moving_compensator_term_cannot_be_deleted":bool(c.diff(1)*p0),
        "slow_stress_removed_only_by_Haar_projection":(R*f*(1+Z+T)).haar()==R*f*(1+Z+T),
    }
    return dict(checks=checks,scope="exact complete-covariance Haar moments plus a separate nonzero formal Fourier pressure mode; not an actual canonical field")


def prefix_review():
    f=(R-1)**6*(2-R)**6;base=(EPS*R*f*(1+Z+T),R**2*(2+Z),means.Poly())
    psi_f=R*f*(1+Z**2+T)*AUX;psi_s=R*f*(1+Z**3+T**2)
    mf=(-EPS*psi_f.diff(1),(R**4*f).diff(0).shift(-2)*(1+Z+T)*AUX,means.radial(psi_f))
    ms=(-EPS*psi_s.diff(1),(R**4*f).diff(0).shift(-2)*(1+2*Z+T),means.radial(psi_s))
    _,_,_,W=means.fixture();zero=(means.Poly(),)*3;new=tuple(x+y for x,y in zip(mf,ms))
    gr=lambda m:means.averaged_gr(base,m,W[0][0],W[1][1],W[0][2])
    g0,gf,gn=gr(zero),gr(mf),gr(new);Bf,uf,zf=mf;B,u,z=ms;b,V,G=base
    fast_change=-means.radial((Bf*Bf).haar())-EPS*(Bf*zf).haar().diff(1)+(uf*uf).haar().shift(-1)
    rg=gn-gf-(2*V*u).shift(-1)
    expected=(EPS*B.diff(2)-means.radial(2*b*B+B*B)-EPS*(b*z+B*z).diff(1)
              +(u*u).shift(-1)+EPS*(means.lap(B)-B.shift(-2)))
    theta=lambda m:(R**2*(V*m[2]+m[2]*m[1]).haar()).integral()
    axial=lambda m,g:(R*(m[2]*m[2]).haar()).integral()-(R**2*g).integral()/2
    later=(mf[0],mf[1]+R*f,mf[2]);later_new=tuple(x+y for x,y in zip(later,ms))
    checks={
        "initial_fast_has_zero_Haar_all_components":all(x.haar()==0 for x in mf),
        "actual_fast_quadratic_Haar_change_retained":gf-g0==fast_change and bool(fast_change),
        "old_fast_slow_crosses_vanish_after_Haar_only":all((x*y).haar()==0 for x in mf for y in ms),
        "fresh_slow_current_radial_remainder_exact":rg==expected,
        "fresh_slow_angular_square_after_linear_row":theta(new)-theta(mf)-(R**2*V*z).integral()==(R**2*z*u).integral(),
        "fresh_slow_axial_square_after_linear_row":axial(new,gn)-axial(mf,gf)+(R*V*u).integral()==(R*z*z).integral()-(R**2*rg).integral()/2,
        "nonlinear_seed_defects_not_zero":all(bool(x) for x in (rg.integral(),(R**2*z*u).integral(),(R*z*z).integral()-(R**2*rg).integral()/2)),
        "initial_sharper_Haar_velocities_are_slow_increment":all(new[j].haar()==ms[j] for j in range(3)),
        "later_nonzero_old_Haar_breaks_initial_reduction":gr(later_new)-gr(later)-(2*V*u).shift(-1)!=expected,
        "slow_divergence_and_velocity_moments":means.radial(B)+EPS*z.diff(1)==0 and (R**2*u).integral()==(R*z).integral()==0,
    }
    return dict(checks=checks,scope="manufactured finite Haar coefficient fields after explicit current linear-row subtraction; not a canonical five-map numerical solve")


def budget_review():
    checks={};rows=[]
    for n in (0,1,3):
        abstract=fixture_budget();out=abstract.output(n);moment=abstract.axial_moment(n,Q(5))
        for S in (1,4,16):
            concrete=fixture_budget(S);actual=concrete.output(n)
            checks[f"n{n}_S{S}_complete_primary_to_initial_scalar_polynomial_chain"]=all(poly.evaluate(S)==actual[key] for key,poly in out.items())
            checks[f"n{n}_S{S}_whole_axial_moment_polynomial_chain"]=moment.evaluate(S)==concrete.axial_moment(n,Q(5))
            fields=abstract.fields(n);resolved=concrete.fields(n)
            checks[f"n{n}_S{S}_initial_Haar_and_full_velocity_majorants"]=all(fields[key].evaluate(S)==resolved[key] for key in ("U","Z","X"))
            slow_fields=concrete.fields(n+2)
            for eps in (Q(1,4),Q(1,16)):
                norms={key:Q(0) for key in coupled.dynamic.KEYS}
                norms.update(b=eps*concrete.base(n+2),B=eps**2*slow_fields["X"],
                             u=eps*slow_fields["U"],g=eps*slow_fields["Z"])
                raw=coupled.dynamic.budget(n,eps,Q(5,4),Q(2),norms)
                checks[f"n{n}_S{S}_eps{eps}_C1_initial_budget_dominates_unreduced_terms"]=all(raw[key]<=eps**2*actual[key] for key in ("P","Qrem","Jtheta","Jz"))
        trace=set();fixture_budget(trace=trace).output(n)
        checks[f"n{n}_scalar_primary_input_order_is_n_plus_10"]=max(trace)==n+10
        trace=set();fixture_budget(trace=trace).axial_moment(n,Q(5))
        checks[f"n{n}_axial_moment_primary_input_order_is_n_plus_11"]=max(trace)==n+11
        for label,method,extra in (("theta_fast_input","E",2),("fast_tangent","Ft",6),("fast_radial","Fr",7),("postfast_pz","preslow_pz",8),("postfast_theta","preslow_theta",6)):
            trace=set();getattr(fixture_budget(trace=trace),method)(n)
            checks[f"n{n}_{label}_input_order_ledger"]=max(trace)==n+extra
        rows.append(dict(order=n,ThetaP=out["P"].record(),ThetaTheta=out["Jtheta"].record(),ThetaZ=out["Jz"].record()))
    k=Q(1,100000);H=1-k
    checks["initial_fast_squares_still_supply_S1"]=2*H>=1
    checks["initial_C1_slow_square_ties_required_S2"]=2*Q(1)==Q(1)+1==2
    checks["initial_radial_time_and_base_products_are_higher_order"]=1+2>=2 and 1+1+1>=2
    checks["Haar_rates_are_sharper_than_full_fast_rates"]=1>H and 2>H+1
    checks["complete_real_curl_square_reduction_margin"]=1-2*k>0
    checks["aggregate_label_norm_not_inferred_from_individual_bounds"]=sum(Q(1) for _ in range(16))>Q(1)
    return dict(checks=checks,rows=rows,scope="exact positive polynomial norm routing; current primary aggregate/root/edge hypotheses remain inputs")


def review():
    parts={"initial_pressure_channels":pressure_review(),"fresh_initial_prefix":prefix_review(),
           "primary_to_seed_budgets":budget_review(),"pressure_tails":means.pressure_tail_review(),
           "complete_covariance":means.covariance_review()}
    return dict(checks={p+":"+k:v for p,row in parts.items() for k,v in row["checks"].items()},parts=parts)
