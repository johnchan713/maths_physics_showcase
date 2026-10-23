"""Reference-scaled normal floor from the SAME retained length proof."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,HERE.parent/path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
tangent=load('reference_floor_tangent','core_tangent_propagator_closure/budgets.py')
normal=tangent.normal
transfer=tangent.transfer
P,exact=transfer.P,transfer.exact

def relative_error_floor(bmin,Bs,error):
    """Rational algebra control; actual Bs may be real in the theorem."""
    b,B,e=map(exact,(bmin,Bs,error))
    if not 0<b<=B or not 0<=e<b:raise ValueError('Valid reference floor and error below bmin required')
    theta=e/b
    return dict(triangle_lower=B-e,relative_lower=(1-theta)*B,
        factor=1-theta,comparison_margin=B-e-(1-theta)*B,
        coarse_uniform_lower=b-e)

def certificate(p,selector_inputs,ell):
    """Recompute ORIGINAL gates; do not select a band or change old inputs.

    This proves an implication for valid original primitive bounds and the
    original Bs calibration. It does not numerically certify those bounds.
    """
    p=normal.frame_inputs(p)
    old=normal.verify_existing(p,selector_inputs,ell)
    if not old['conditional_normal_floor_implication_passes']:
        raise ValueError('Supplied original band and original length semantics must pass')
    b=old['original_budget']['original_budgets'];L=exact(selector_inputs['length'])
    theta=b['Dt']/(L*p['bmin']);cap=1/(32*(1+b['Sref']))
    beta=p['lam_min']/(4*(1+p['u']**2)**2)
    qref=1/((1-theta)**2*beta);qold=old['inverse_normal_square']
    quse=min(qold,qref)
    gates=dict(original_normal_and_selector=old['conditional_normal_floor_implication_passes'],
        original_Riccati_slack_gives_relative_error=0<theta<=cap<=Q(1,32),
        positive_phase_beta=beta>0,
        calibrated_reference_inverse_improves_requested_phase_bound=qref<=Q(1024,961)/beta<4/beta,
        original_inverse_bound_remains_available=0<quse<=qold and quse<=qref)
    if not all(gates.values()):raise ValueError('Reference-scaled original-length implication failed')
    return dict(gates=gates,original_normal_certificate=old,existing_ell=ell,
        existing_Lreq=L,relative_error_fraction=theta,relative_error_cap=cap,
        reference_normal_factor=1-theta,phase_beta=beta,
        reference_scaled_inverse=qref,uniform_reference_scaled_inverse=Q(1024,961)/beta,
        old_inverse_normal_square=qold,inverse_normal_square=quse,
        phase_beta_bound_derived_from_full_error_and_reference_calibration=True,
        phase_beta_bound_inferred_from_coarse_floor_alone=False,
        extra_native_margin_gate_required_for_this_original_frame_route=False,
        original_Mz_Db_frame_and_calibration_validity_retained=True,
        actual_primitive_enclosures_or_unknown_historical_band_certified=False,
        bands_lengths_norm_bounds_roots_or_fields_selected_or_changed=False)

def forced_from_reference(p,selector_inputs,ell,*,F,G,btilde,path,source,harmonic=1):
    """Rebind only the positive inverse majorant; the actual inverse is unchanged.

    The existing tangent theorem supplies its original T and length. The
    materializer is for supplied manageable inputs; giant constants may
    remain the exact symbolic recipes in the analytical statement.
    """
    c=certificate(p,selector_inputs,ell)
    phase,old,tc=tangent.forced_from_frame(p,selector_inputs,ell,
        F=F,G=G,btilde=btilde,path=path,source=source,harmonic=harmonic)
    new=transfer.inverse.ForcedBudget(normal=old.normal,moving_normal=old.moving_normal,
        matrix=old.matrix,source=old.source,profile_rate=old.profile_rate,
        tangent=old.tangent,length=old.length,
        inverse_normal_square=c['inverse_normal_square'],harmonic=old.harmonic)
    return phase,new,old,c,tc

def coefficientwise_leq(a,b):
    if not isinstance(a,P) or not isinstance(b,P):raise TypeError('Positive scale polynomials required')
    return all(a.c.get(k,Q(0))<=b.c.get(k,Q(0)) for k in set(a.c)|set(b.c))

def certificate_scope():
    return dict(parent_commit='0ddb0ef928e3e51b76173f847a18b150e3687cba',
        correction='optional low-order native margin is not a prerequisite for the original full-frame route',
        core='|n_tan-Bs K|<=Dt/L<=theta*bmin<=theta*Bs; original Lreq>=LR gives theta<=1/[32(1+3u/2)]',
        calibration='Bs^2=lambda0/[epsilon*k^2*(1+u_actual^2)^(3/2)]>=phase_beta=lambda_min/[4(1+u_upper^2)^2]',
        conclusion='|n|^-2<=1/[(1-theta)^2 phase_beta]<=1024/(961 phase_beta); min with original 4/bmin^2',
        recurrence='all same ForcedBudget recipes are coefficientwise nondecreasing in the reciprocal constant with other inputs fixed',
        next='complete generated path/source and weighted target/freezing/root/cone estimates, full covariance/mean and residual',
        actual_primitive_enclosures_complete_residual_or_blowup_certified=False)
