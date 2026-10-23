"""Current angular-row normal form; actual field/row/support gates stay open."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
mixed=load('angular_pressure_mixed',HERE.parent/'core_mixed_compact_time_budget/bounds.py')
parent=mixed.parent
P,order,exact=mixed.P,mixed.order,mixed.exact

def temporal_recipe(n):
    """Exact deferred 2 Q * 6^n Gamma_n (n+4)! raw/through conversion."""
    order(n)
    return dict(output_order=n,raw_source_order=n+4,scalar=48,two_power=n+4,
        six_power=n,source_factorial=n+4,
        factors=['original c^-1 floor Tg S','original weighted/raw temporal domain and profile conversions'],
        materialized=False)

def rates(kappa=Q(1,100000)):
    r=parent.rates(kappa);k=exact(kappa);B,E,C0=(r[x] for x in ('B','E','C0'))
    J=dict(bv=B+1,Vbeta=B+1,beta_v=2*B+1,Wrt=Q(1))
    F=dict(current_angular_row=C0,slow_time=B+1,axial_flux=B+1,
           full_vector_viscosity=B+1-2*k,axial_viscosity=B+3)
    return dict(B=B,E=E,C0=C0,J_groups=J,F_groups=F,J=min(J.values()),F=min(F.values()),
        quadratic=2*B,rest=Q(1),H=Q(1),K=E,local_pressure=Q(1),reduced_pressure_source=Q(1),
        current_mean_pressure=Q(1),whole_regular=E,parent_A1=B+1,
        actual_inputs_absolute_decay_or_blowup_certified=False)

def request_counts(m,*,physical_power=None,h=Q(1,200),kappa=Q(1,100000)):
    order(m);k=exact(kappa);r=rates(k)
    t=r['C0'] if physical_power is None else parent.compact_physical_orders(0,m,physical_power,h,k)['normalized_target']
    out={}
    for name,a,extra in (('J',r['J'],2),('F',r['F'],2),('quadratic',r['quadratic'],3),('rest',r['rest'],3)):
        p=parent.mean.strict_integrations(a,t,k)
        out[name]=dict(source_order=a,integrations=p,raw_order=m+2*p+extra,strict_margin=a+p*k-t)
    out.update(target=t,output_order=m,physical_power=physical_power,
        J_multiplier_V_order=out['J']['raw_order']+1,
        F_current_field_order=out['F']['raw_order']+2,
        rest_beta_order=out['rest']['raw_order']+2,
        older_entire_DAG_retained=True,actual_input_or_giant_constant_evaluated=False)
    return out

def pressure_envelope(n,data):
    """Positive factorial-through-order bound; current fields through n+6.

    temporal bounds c^-1 Nt^-1 Q INCLUDING four-jet/raw conversions and
    original weights/profile domains. All scalar inputs are nonnegative
    scaled norm bounds, not newly assumed actual field enclosures.
    """
    order(n)
    keys=('x','u','z','uHaar','b0','V','G','r','graph','c','Etheta','W','temporal','KT0','Kmom','rho')
    a=parent.mean.prefix.require(data,keys)
    x,u,z,u0,b,V,G,r,g,c,e,W,I,KT,KM,rho=(a[k] for k in keys)
    d=n+7;D=3*d*(1+g)
    av=2*V*r;aq=6*V*r*r+2*d*g*V*r
    J=b*u+V*x+x*u+W
    axial=G*u+V*z+z*u+W
    vis=D*D+r*D+r*r+d*d
    F=e+d*u+d*axial+vis*u
    H=I*J;K=I*F
    S=2*b*x+x*x+W;rz=b*z+G*x+x*z+W
    rest=d*(1+2*c)*x+r*S+d*rz+r*W+d*d*x
    local=S+(D+r)*x+av*H
    source=av*u0+r*u*u+rest+aq*H+av*K
    moment=KM*source;pressure=local+KT*(source+rho*moment)
    return dict(J=J,F=F,H=H,K=K,rest=rest,local=local,source=source,
        moment=moment,pressure=pressure,current_field_input_order=n+6,
        current_angular_row_input_order=n+4,temporal_input_cost=4,temporal_constant=temporal_recipe(n),
        angular_row_is_existing_derived_final_theta_not_a_new_fine_callback=True)

def certificate(m=0):
    order(m);r=rates()
    return dict(parent_commit='29195306c4bccab9dc19dd5528ae03c2a6b84d51',
        current_stage='FINAL current means, SAME complete current wave covariance and selected stress',
        angular_row='Etheta=tstar v+D2 J+Dz L-epsilon(D_r D1 v+Dz^2 v)-D2 Sigma_theta',
        J='b v+V beta+beta v+W_rtheta',L='G v+V gamma+gamma v+W_theta,z',
        F='-Etheta-epsilon partial_T v+Dz L-epsilon(D_r D1 v+Dz^2 v)',
        H='c^-1 Nt^-1 Qcirc J',K='c^-1 Nt^-1 Qcirc F',
        exact_swirl='Qcirc v=-D2 H-K; Sigma_theta is auxiliary-independent',
        a='2V/R',a1='2a/R-partial_R a=6V/R^2-2V_R/R',
        local='K1=K0-a H, K0=-S+epsilon D1 beta',
        source='r1=a Haar(v)+v^2/R+Rrest-a1 H-a K',
        radial='gr=r1+D_r K1; A0 Qcirc gr=-A0(a1 H)-A0(a K)+A0 Qcirc(v^2/R+Rrest)',
        pressure='P=int Haar r1; pm=K1+T0(r1-rho P), SAME moving rho, base Pi1/Pi2 and all old full pressure memory retained',
        Haar='Haar r1=a Haar(v)+Haar(v^2)/R+Haar Rrest; variance and all old/new correlations retained',
        noncircular='use the already derived final ANGULAR row C0; its equation and exact fast cancellation have no current pm or current A0 tail. Older A0 dependencies remain in the archived source graph',
        current_mean_pressure_order=str(r['current_mean_pressure']),whole_regular_order=str(r['whole_regular']),parent_direct_A1_order=str(r['parent_A1']),
        pressure_inputs='original weighted current fields/covariance, sharper Haar v order1, derived final Etheta C0, all raw temporal/profile/compact operator and WHOLE support contracts; not actually enclosed',
        domains='retain entire previous FINAL12q common-pole graph, all archived shifts and source/output weight conversions; add fixed multiplier poles before selecting q',
        normalized_requests=request_counts(m),physical_requests_illustration=request_counts(m,physical_power=0),
        total_original_minimal_jet_actual_source_validity_absolute_decay_or_blowup_proved=False)
