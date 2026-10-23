"""Joint radial/temporal divisor estimate and the current centripetal split.

All actual coarse field, weighted/profile and whole support gates remain.
No giant norm program or actual canonical source is materialized here.
"""
from fractions import Fraction as Q
from math import ceil
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

quotient=load('mixed_time_quotient',HERE.parent/'core_compact_source_quotient/bounds.py')
sharp,parent=quotient.sharp,quotient.parent
P,order,exact=quotient.P,quotient.order,quotient.exact

def divisor_data(k):
    if len(k)!=2 or any(type(v) is not int for v in k) or k==(0,0) or list(k)==[0,0]:
        raise ValueError('Nonzero integer torus mode required')
    u,v=k;N=(u+v)**2-2*v*v
    return dict(k=list(k),l1=abs(u)+abs(v),quadratic_norm=N,
        mixed_product='(sqrt(2)-1)*quadratic_norm',
        lower_product='sqrt(2)-1',same_mode_in_both_divisors_required=True)

def kernel_recipe(m,p,*,e=0):
    order(m);order(p)
    if p<1 or e not in (0,1):raise ValueError('Original e=0/1 and positive integer integration count required')
    J=m+2*p+2
    return dict(operator='A_e(a*c^-1*Nt^-1 F)',e=e,output_order=m,integration_order=p,
        raw_F_input_order=J,source_factorial=J,
        scalar=8,two_power=2*m+J+p+2,three_power=p,
        factors=['c^-1','original U-hull length L',f'N_{m}(k_e chi_U)',f'N_{J}(f_e*a)',f'N_{J}(F)'],
        multiplier_contract='a is auxiliary-independent; retain all its raw derivatives, not an arbitrary multiplier norm',
        Fourier_tail='r^-3 after p radial integrations and a_y+p+2 torus derivatives',
        sequential_F_input_order=m+2*p+7,saved_F_derivatives_vs_sequential=5,
        original_width_gap_weights_c_floor_M_floor_and_all_archived_shifts_retained=True,
        actual_input_or_factorial_evaluated=False)

def sharp_boundary(m,J,p):
    order(m);order(J);p=exact(p)
    if p<0:raise ValueError('Nonnegative M gain required')
    gap=m+2*p-J-1
    return dict(output_order=m,F_input_order=J,M_gain=p,
        ratio_power_in_P=gap,strict_finite_jet_obstruction=gap>0,
        necessary_F_order=max(0,ceil(m+2*p-1)),
        equality_or_optimal_upper_rule_resolved=False,
        fast_map_compliant_comparison_not_actual_canonical_source=True)

def radial_profile_recipe(J,d,u0):
    """Exact all-jet upper recipe for f0/R=1/(d U), without allocating J!."""
    order(J);d,u0=map(exact,(d,u0))
    if min(d,u0)<=0:raise ValueError('Original positive d and U-hull floor required')
    return dict(profile='f0/R=1/(d*U)',raw_order=J,factorial_index=J,
        prefactor=1/(d*u0),base=max(Q(1),1/u0),power=J,
        exact_derivative='(-1)^j j!/(d U^(j+1))',
        bound='prefactor * factorial(factorial_index) * base^power',
        numerical_factorial_or_actual_source_materialized=False)

def radial_counts(m,*,physical_power=None,h=Q(1,200),kappa=Q(1,100000)):
    order(m);r=parent.rates(kappa);k=exact(kappa)
    target=r['C0'] if physical_power is None else parent.compact_physical_orders(0,m,physical_power,h,k)['normalized_target']
    pc=parent.mean.strict_integrations(r['B'],target,k)
    pr=parent.mean.strict_integrations(Q(1),target,k)
    JG=m+2*pc+2;Jr=m+2*pr+3;old=m+2*pc+3
    return dict(output_order=m,physical_power=physical_power,normalized_target=target,
        centripetal_source_order=r['B'],rest_source_order=Q(1),
        mixed_integrations=pc,rest_integrations=pr,G_raw_order=JG,rest_raw_order=Jr,
        generic_Fsum_from_v_order=JG+1,rest_beta_order=Jr+2,
        local_current_field_order=max(JG+1,Jr+2),
        previous_full_gr_raw_order=old,previous_full_gr_field_program_order=old+3,
        mixed_strict_margin=r['B']+pc*k-target,rest_strict_margin=1+pr*k-target,
        mixed_epsilon_target_necessary_F_order=max(0,ceil(m-1+2*(target-r['B'])/k)),
        older_root_normal_source_DAG_expanded_or_erased=False,
        sufficient_physical_transfer_loss_not_proved_necessary=True)

def split_envelope(n,data):
    """Positive through-n G/rest envelopes from fields THROUGH n+2.

    b=O(epsilon), V,G fixed multiplier O(1), beta M_(B+1),
    v,gamma M_B, complete covariance M1. Fsum=-c Nt v is derived,
    not an independently supplied desired old-row callback.
    """
    order(n);keys=('x','u','z','b','V','G','c','r','W')
    data=parent.mean.prefix.require(data,keys)
    x,u,z,b,V,G,c,r,W=(data[k] for k in keys);d=n+3
    Fsum=2*c*d*u
    Gtheta=2*(V+u)*Fsum
    S=2*b*x+x*x+W;H=b*z+G*x+x*z+W
    rest=d*(1+2*c)*x+r*S+d*H+r*W+d*d*x
    return dict(Fsum=Fsum,Gtheta=Gtheta,rest=rest,
        Gtheta_input_order=n+1,rest_input_order=n+2,
        Gtheta_order='B',rest_order='1',
        Haar_pressure='2V Haar(v)+Haar(v^2), including every fast variance and old/fast cross; not reconstructed from Nt inverse')

def certificate(m=0):
    order(m)
    return dict(parent_commit='57aabbb712e34ce4e4a81363281a924f09f63308',
        exact_divisor_identity='(vr.k)(vt.k)=(sqrt(2)-1)*[(k1+k2)^2-2k2^2], nonzero integer bracket',
        mixed_rule='|vr.k|^-p |vt.k|^-1 <= 3^p |k|_1^(p-1), p>=1',
        current_stage='FINAL current means after latest eight-mode wave and fresh fast/slow',
        Fsum='-c Nt v; actual fast residual sum only with original cover/label/c normalization intertwiners retained; generic bound uses one v derivative',
        Gtheta='2(V+v)Fsum, V auxiliary-independent, v the FULL accumulated mean including all slow/fast contributions',
        radial_identity='A0 Qcirc gr=-c^-1 A0(R^-1 Nt^-1 Qcirc Gtheta)+A0 Qcirc Rrest',
        Rrest='-tstar beta-S/R-Dz H+Wtt/R+epsilon Dz^2 beta',
        pressure='r0=Haar(2Vv+v^2)/R-c^-1 R^-1 Nt^-1 Qcirc Gtheta+Rrest; P=int Haar r0; pm=K0+T0(r0-rho P), same moving rho and K0=-S+epsilon D1 beta',
        no_free_commutation='Do not pull an auxiliary-dependent multiplier through Nt^-1; Nt(a h)=a Nt h+(Nt a)h. No new slow commutation is used.',
        domains='Gtheta and Rrest at final12q under retained fixed multiplier contracts; original raw profiles, f0/R reciprocal jets, output edge floor and ALL old common shifts/request DAGs retained',
        actual_new_profile='f0/R=1/(d U); N_J <= J! max(1,u0^-1)^J/(d u0) on the SAME original U>=u0>0; no giant factorial evaluated',
        normalized_requests=radial_counts(m),physical_requests_illustration=radial_counts(m,physical_power=0),
        complete_regular_order=str(parent.rates()['E']),parent_A1_direct_order=str(quotient.rates()['actual_A1_error']),
        canonical_row_high_jet_or_support_norms_absolute_decay_infinite_field_or_blowup_validated=False)
