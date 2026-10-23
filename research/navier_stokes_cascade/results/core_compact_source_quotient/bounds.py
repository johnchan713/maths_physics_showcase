"""Exact compact-source quotient and the actual axial-potential A1 gain.

All norm assertions retain the original coarse actual input/support gates.
No current canonical norm callback or actual input enclosure is supplied.
"""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

sharp=load('source_quotient_sharp',HERE.parent/'core_compact_input_sharpness/bounds.py')
parent=sharp.parent
P,order,exact=parent.P,parent.order,parent.exact

def rates(kappa=Q(1,100000)):
    r=parent.rates(kappa);k=exact(kappa);B,E,C0=(r[x] for x in ('B','E','C0'))
    radial=dict(time_beta=B+1,S_over_R=Q(1),axial_H=min(B+2,2*B+2,Q(2)),
        centripetal=B,axial_viscosity_beta=B+4)
    return dict(B=B,E=E,C0=C0,reduced_radial_groups=radial,reduced_radial_source=min(radial.values()),
        complete_reduced_pressure_flux=B,actual_A1_error=B+1,
        gain_over_parent_A1_source=B+1-E,margin_over_C0=B+1-C0,
        normalized_oscillatory_A1_integrations=0,
        full_pressure_or_regular_source_order_promoted=False,
        actual_input_validity_absolute_decay_or_blowup_certified=False)

def reduced_envelope(n,data):
    """A1 output through n from coarse inputs THROUGH n+3.

    Products are submultiplicative in the inherited five-variable norm.
    The entire finite coefficient list dominates every fixed/moment/profile
    call through n+3 on one common weighted chart. W bounds ALL six entries
    of the complete current covariance. No radial graph factor is needed
    in r0. This does not delete its pressure primitive K0 from the field.
    """
    order(n)
    keys=('x','u','z','b','V','G','c','r','rho','Kp','Kmom','KA1','W')
    data=parent.mean.prefix.require(data,keys)
    x,u,z,b,V,G,c,r,rho,Kp,Kmom,KA1,W=(data[k] for k in keys)
    d=n+4
    S=2*b*x+x*x+W;H=b*z+G*x+x*z+W;C=2*V*u+u*u+W
    r0=d*(1+2*c)*x+r*S+d*H+r*C+d*d*x
    pressure=Kp*(r0+rho*Kmom*r0)
    h1=2*G*z+z*z+W+S+pressure
    return dict(S=S,H=H,C=C,r0=r0,P=Kmom*r0,tilde_pressure=pressure,
        H1=h1,A1_coefficient=2*KA1*d*h1,field_input_order=n+3,
        internal_r0_pressure_H1_order=n+1,explicit_epsilon_power=1,
        source_is_full_pressure_flux_not_reduced_fine_E_row=True)

def A1_request(m,*,physical_power=None,h=Q(1,200),kappa=Q(1,100000)):
    order(m);r=rates(kappa);a=r['actual_A1_error']
    target=r['C0'] if physical_power is None else parent.compact_physical_orders(1,m,physical_power,h,kappa)['normalized_target']
    p=max(0,(target-a)//exact(kappa)+1)
    # Direct fixed-profile operator for p=0; positive-p raw kernel (8)
    # on Dz H1 needs ONE extra H1 derivative. No factorial is evaluated.
    J=m+1 if p==0 else m+2*p+4
    old=parent.rates(kappa)['A1_integrations'] if physical_power is None else parent.compact_physical_orders(1,m,physical_power,h,kappa)['integration_order']
    return dict(output_order=m,physical_power=physical_power,effective_source_order=a,
        normalized_target=target,integration_order=p,raw_H1_input_order=J,
        strict_margin=a+p*exact(kappa)-target,
        previous_full_E_row_integration_order=old,previous_full_E_row_raw_order=m+2*old+3,
        original_profile_and_weighted_KA1_at_6q_retained=True,
        older_A0_source_root_normal_and_physical_request_DAGs_erased=False)

def certificate(n=0):
    order(n)
    return dict(
        parent_commit='e441859258f915144b59ca4f967b10e28d49c8f0',
        exact_operator_identities=['A_e D_e=0','T_e D_e=I on compact potentials','D_e T_e=I-A_e','A_e^2=A_e','T_e A_e=0'],
        boundary_rule='A_e D_b h=R^-e chi_R [boundary(r^e h along the full drift)+ (b-e) J_e(h/r)]',
        radial=dict(stage='CURRENT FINAL means after latest eight-mode wave and fresh fast/slow',
            identity='gr=r0+D_r K0; K0=-S+epsilon D1 beta',
            r0='-tstar beta-S/R-Dz H+C/R+epsilon Dz^2 beta',
            pressure='P=int Haar r0; pm=K0+T0(r0-rho P), SAME moving unit-mass rho',
            A0='A0 Qcirc gr=A0 Qcirc r0, provided the WHOLE K0 support/extension contract holds',
            source_order=str(rates()['B'])),
        axial=dict(stage='CURRENT POSTWAVE Ez BEFORE the latest fresh fast/slow; do not mix stages',
            exact_mean_route='gamma=D1 Psi, beta=-Dz Psi, Psi the finite sum of SAME actual compact fast/slow potentials',
            original_row='Ez=tstar gamma+D1 H+Dz(2G gamma+gamma^2+Wzz+pm)-epsilon(D1 D_r gamma+Dz^2 gamma)-D1 Sigma_z',
            identity='A1 Qcirc Ez=A1 Dz Qcirc H1',
            H1='2G gamma+gamma^2+Wzz-S+T0(r0-rho P)',
            radial_pressure_memory='-S remains; epsilon D1 beta is retained in pm but its Dz image is annihilated ONLY under A1',
            full_forced_wave_covariance_and_old_pressure_memory_retained=True,
            error_order=str(rates()['actual_A1_error'])),
        moving_cutoff='[partial_s,A_e]f=R^-e chi_Rs J_e f for s=Z,T; never divide by chi_R or commute it without this term',
        domains='postwave old means2q, complete current wave3q, covariance6q, reduced r0/pressure/H1 at6q; SAME Kp/Kmom/KA1 at6q and all original shifts. Final r0 remains12q with original A0 raw profiles. Retain unweighted fixed multiplier contracts, or add every multiplier pole at products before choosing q.',
        normalized_A1=A1_request(n),physical_A1_illustration=A1_request(n,physical_power=0),
        radial_obstruction='pure swirl comparison is divergence free with both exact Haar velocity moments, yet its fundamental centripetal A0 source remains Pell resonant; no canonical-source counterexample',
        complete_regular_source_order=str(rates()['E']),
        actual_high_jet_support_inputs_or_absolute_residual_improvement_validated=False)
