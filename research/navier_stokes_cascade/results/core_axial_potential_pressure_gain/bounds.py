"""Axial potential normal form and the SAME postwave A1 pressure flux."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
angular=load('axial_potential_angular',HERE.parent/'core_angular_pressure_normal_form/bounds.py')
parent=angular.parent;quotient=angular.mixed.quotient
P,order,exact=angular.P,angular.order,angular.exact
ANGULAR_KEYS=('x','u','z','uHaar','b0','V','G','r','graph','c','Etheta','W','temporal','KT0','Kmom','rho')
POTENTIAL_KEYS=('x','z','zHaar','b0','G','r','graph','Ez','pmCoarse','W','KT1','temporal')
JOIN_KEYS=tuple(dict.fromkeys(ANGULAR_KEYS+POTENTIAL_KEYS+('KA1',)))
COMPLEMENT_KEYS=tuple(k for k in JOIN_KEYS if k not in ('uHaar','zHaar','Kmom','rho'))

def stage_rates(stage='postwave',kappa=Q(1,100000)):
    r=angular.rates(kappa);B,E,C0=(r[x] for x in ('B','E','C0'))
    if stage not in ('postwave','final'):raise ValueError('Use the SAME postwave or final current state')
    a=E if stage=='postwave' else C0
    groups=dict(current_Ez=a,T1_Dz_X=B+1,slow_potential_time=B+1,
        radial_flux=Q(1),radial_viscosity=B+1-exact(kappa),axial_viscosity=B+3)
    return dict(stage=stage,row_order=a,pole_multiplier=6 if stage=='postwave' else 12,
        potential_source_groups=groups,potential_source=min(groups.values()),
        Q_potential=Q(1),Haar_potential=Q(1),full_potential=Q(1),full_beta=Q(2),
        gamma=B,Q_mean_pressure=Q(1),projected_axial_flux=Q(1),A1_error=Q(2),
        previous_A1_error=B+1,whole_regular=E,
        sharper_Haar_history_required_only_for_full_potential_and_beta=True,
        actual_all_band_input_support_or_blowup_certified=False)

def A1_request(m,*,physical_power=None,h=Q(1,200),kappa=Q(1,100000)):
    order(m);r=stage_rates('postwave',kappa);k=exact(kappa)
    target=angular.rates(k)['C0'] if physical_power is None else parent.compact_physical_orders(1,m,physical_power,h,k)['normalized_target']
    p=max(0,(target-2)//k+1);J=m+1 if p==0 else m+2*p+4
    old=quotient.A1_request(m,physical_power=physical_power,h=h,kappa=k)
    return dict(stage='postwave',output_order=m,physical_power=physical_power,source_order=Q(2),
        normalized_target=target,integration_order=p,raw_Phi_order=J,strict_margin=2+p*k-target,
        previous_raw_H1_order=old['raw_H1_input_order'],previous_source_order=old['effective_source_order'],
        current_fields_for_potential=J+6,current_Ez_for_potential=J+4,
        coarse_pressure_for_potential=J+5,coarse_pressure_conservative_field_recipe=J+8,
        current_Etheta_for_fine_pressure=J+4,fixed_G_for_GR=J+1,
        entire_older_A0_root_normal_and_source_graph_retained=True,
        giant_factorials_or_actual_source_norms_evaluated=False)

def direct_A1_recipe(m):
    order(m)
    return dict(output_order=m,raw_Phi_order=m+1,scalar=2,six_power=m,source_factorial=m+1,
        factors=['epsilon','ORIGINAL raw KA1_m at postwave6q','factorial-through norm Phi_(m+1)'],
        cutoff_not_commuted_past_Dz=True,materialized=False)

def potential_envelope(n,data,stage='postwave'):
    """Potential at n: fields/coarse T1 through n+6, Ez n+4, old pm n+5.

    Uses the old M_B pressure bound; no fine pm/desired-potential callback.
    KT1 dominates the FINITE list of original required orders at ONE pole.
    temporal is the parent's explicit four-jet raw/through recipe.
    """
    order(n);sr=stage_rates(stage)
    a=parent.mean.prefix.require(data,POTENTIAL_KEYS)
    x,z,zh,b,G,r,g,ez,pm,W,KT,I=(a[k] for k in POTENTIAL_KEYS)
    d=n+7;D=3*d*(1+g)
    psi0=KT*z;X=2*G*z+z*z+W+pm
    Hrz=b*z+G*x+x*z+W
    A=KT*(ez+d*X)+d*psi0+Hrz+D*z+d*d*psi0
    circ=I*A;haar=KT*zh
    return dict(coarse_potential=psi0,coarse_X=X,Hrz=Hrz,A=A,Q_potential=circ,
        Haar_potential=haar,potential=circ+haar,stage=stage,pole_multiplier=sr['pole_multiplier'],
        current_fields=n+6,current_Ez=n+4,old_coarse_pressure=n+5,
        no_new_fine_output_callback=True)

def joined_envelope(n,data,stage='postwave'):
    """Phi=-2G_R Psi+gamma^2+Wzz+pm, with BOTH fine outputs derived."""
    order(n);stage_rates(stage);a=parent.mean.prefix.require(data,JOIN_KEYS)
    pot=potential_envelope(n,{k:a[k] for k in POTENTIAL_KEYS},stage)
    # The angular lemma needs Etheta>=E>1, so applies at BOTH stages.
    pressure=angular.pressure_envelope(n,{k:a[k] for k in ANGULAR_KEYS})
    GR=(n+7)*a['graph']*a['G']
    Phi=2*GR*pot['potential']+a['z']*a['z']+a['W']+pressure['pressure']
    noHaar={k:a[k] for k in ANGULAR_KEYS};noHaar['uHaar']=a['uHaar']*0
    qpressure=angular.pressure_envelope(n,noHaar)
    # Q rho P=0, Q(a Haar v)=0. Q costs two on the other local/source
    # terms. The zero Haar placeholder is NOT a full-pressure estimate.
    pmQ=2*(qpressure['local']+a['KT0']*qpressure['source'])
    PhiQ=2*GR*pot['Q_potential']+2*a['z']*a['z']+2*a['W']+pmQ
    return dict(potential=pot,pressure=pressure,Phi=Phi,Q_pressure=pmQ,Q_Phi=PhiQ,stage=stage,
        same_current_pm_and_full_covariance_required=True)

def complement_envelope(n,data,stage='postwave'):
    """Q-only route: no finer Haar input, moment bound or rho norm required."""
    order(n);a=parent.mean.prefix.require(data,COMPLEMENT_KEYS)
    full={k:(a[k] if k in a else a['x']*0) for k in JOIN_KEYS}
    out=joined_envelope(n,full,stage)
    return dict(Q_potential=out['potential']['Q_potential'],Q_pressure=out['Q_pressure'],
        Q_Phi=out['Q_Phi'],stage=stage,current_fields=n+6,current_Ez=n+4,
        coarse_pressure=n+5,sharper_Haar_inputs_used=False)

def tail_envelope(m,data,stage='postwave'):
    order(m)
    if stage!='postwave':raise ValueError('The actual latest fast A1 tail uses POSTWAVE fields, not final means')
    a=complement_envelope(m+1,data,stage)
    # KA1 is 6^m m! times its ORIGINAL raw profile constant, with weights.
    return dict(A1=2*(m+1)*data['KA1']*a['Q_Phi'],joined=a,
        explicit_A1_epsilon_power=2,
        through_order_operator_recipe=direct_A1_recipe(m))

def beta_envelope(m,data,stage='final'):
    order(m);a=parent.mean.prefix.require(data,POTENTIAL_KEYS)
    out=potential_envelope(m+1,a,stage)
    return dict(beta=(m+1)*out['potential'],explicit_epsilon_power=2,
        potential=out,sharper_Haar_gamma_history_required=True)

def certificate(m=0):
    order(m)
    return dict(parent_commit='1999c84d2476b7b332a08bc108d158078f3a7f6c',
        actual_potential='gamma=D1 Psi, beta=-Dz Psi; SAME finite sum of compact actual fast/slow potentials; T1 gamma=Psi',
        X='2G gamma+gamma^2+Wzz+pm, COMPLETE current pm including epsilon D1 beta',
        exact_axial_row='Ez=D1[tstar Psi+Hrz-epsilon D_r gamma-epsilon Dz^2 Psi-Sigma_z]+Dz X',
        normal_form='Q Psi=c^-1 Nt^-1 {T1 Q(Ez-Dz X)+epsilon partial_T Q Psi-Q Hrz+epsilon D_r Q gamma+epsilon Dz^2 Q Psi}',
        noncircular='every Psi and pm on the RHS uses the OLD coarse M_B recipe, and Ez is the already derived SAME-stage row E or C0; mean-pressure improvement is independently from the angular row',
        pressure_viscosity='if Xflat=X-epsilon D1 beta, the explicit axial potential term must be -2epsilon Dz^2 Psi instead of -epsilon Dz^2 Psi',
        Haar='Haar Psi=T1 Haar gamma at order1 from preserved finite slow sum, not from relabelled full gamma B',
        full_potential_order='1',full_beta_order='2',full_gamma_order=str(stage_rates()['gamma']),
        projected_identity='A1 Q Ez=A1 Dz Q Phi; Phi=-2 G_R Psi+gamma^2+Wzz+pm; Q Phi order1 needs only Q Psi and Q pm, not finer Haar inputs; never commute moving A1 through Dz',
        omitted_radial_potential='2G gamma=D1(2G Psi)-2G_R Psi; G_R retained on WHOLE compact potential hull',
        postwave=dict(row_order=str(stage_rates()['row_order']),operator_pole='6q',actual_A1_error_order='2',
            source='SAME postwave source before latest fast/slow; no final pressure substitution'),
        final=dict(row_order=str(stage_rates('final')['row_order']),operator_pole='12q',potential_order='1',beta_order='2'),
        pressure='Q pm=Q K1+T0 Q r1, Q rho P=0 and Q(a Haar v)=0; only this channel removes rho P. FULL pm still SAME stage K1+T0(r1-rho P), complete moving rho P, full variance, axial pressure memory and retained base Pi1/Pi2; no total-pressure estimate',
        normalized=A1_request(m),physical_illustration=A1_request(m,physical_power=0),
        whole_regular_order=str(stage_rates()['whole_regular']),
        actual_all_scale_high_jet_support_norms_absolute_decay_infinite_field_or_blowup_proved=False)
