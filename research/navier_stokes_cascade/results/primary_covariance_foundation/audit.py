#!/usr/bin/env python3
"""Primary covariance positivity, normalized inverse, weighted roots and limits."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from covariance_bounds import corner_review,scale_review
from covariance_checks import jet_review,gaussian_review,edge_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='explicit-primary-covariance-root-gap-inverse-and-weighted-jets-from-reviewed-leading-frame-inputs'
DEPENDENCIES=('sourced_inverse_foundation','global_pressure_moment_foundation','weighted_stress_foundation')
FLAGS=dict(explicit_positive_primary_covariance_roots_from_reviewed_columns=True,
    sharp_entrywise_error_threshold_gamma_over_two_minus_gamma=True,
    safe_normalized_inverse_below_three_with_Haar_radius_costs=True,
    finite_homogeneous_ODE_existence_from_parent_Picard_argument=True,
    freezing_and_length_budgets_independent_of_derivative_order=True,
    original_leading_target_root_branches_frequency_and_width_retained=True,
    all_fixed_root_and_reciprocal_jets_keep_exact_half_edge_weight=True,
    actual_outer_normalized_target_polynomial_pole_retained=True,
    selected_principal_signed_cross_covariance_not_full_tensor_update=True,
    full_radial_axial_and_Pi1_Pi2_mean_pressure_obligations_retained=True,
    historical_evidence_unchanged_and_controls_not_actual_huge_fields=True,
    normalized_target_or_root_quotient_bounded_on_closed_shell=False,
    strict_positive_root_gap_at_zero_stress_endpoints=False,
    all_retained_core_background_and_frame_floors_independently_proved=False,
    all_generated_complete_covariance_mean_and_residual_inputs_closed=False,
    unconditional_actual_pulse_family_or_blowup_verified=False,
    formal_certificate=False)
REMAINING=[
    'independently_close_retained_exact_core_actual_background_and_frame_floor_inputs',
    'join_actual_generated_complete_covariance_mean_targets_including_radial_axial_velocity_and_Pi1_Pi2',
    'verify_all_viscous_cutoff_interaction_residuals_and_single_field_realization_before_unconditional_claims']


def protocol():
    return dict(parent_commit='b17351e44a98b0ca8ef02a936ea7a0d0e94a4fc4',
        fixed_choices='same Md64,h,N,positive matching lambda,unit viscosity,original primary roots,r0,frequency,width C^-20,supports and gauge',
        normalized_matrix='C=[[1+a,1+b],[-1+c,1+d]]; entry errors <=eta; target (V,W), V>0, |W|<=(1-gamma)V',
        sharp_positive_threshold='eta<gamma/(2-gamma); determinant 2 witness loses positivity at equality',
        safe_budget='eta<=gamma/4: detC>=9/8, inverse2<=20/9<3, gamma*V/8<=h_sigma*y_sigma<=3V',
        freezing='tau<=gamma0/8 retains gamma=gamma0/2; exact half-margin threshold gamma0/(4-gamma0)',
        length='L>=[4*Cdir*max(Amin^-1,u^-1)/gamma]^2; existing padded length dominates',
        Haar_scale='s0=r0^2/sqrt(L); H_inverse2<=3/(s0*jmin*min(A,u))',
        weighted_roots='y=(zeta/s0)*b; b>=gamma*vstar/(8*jmax); sqrt(y)=sqrt(zeta)*L^(1/4)/r0*sqrt(b)',
        edge='outer T0/zeta has delta^-3 pole; root/sqrt(zeta) may grow delta^-3/2; half constants t1^2/2 and 2 unchanged',
        signed='d=H^-1(Sigma/epsilon); delta_a=d/(2*a_original); no updated square root; retain signed square',
        scope='finite primary root foundation relative to reviewed leading and frame inputs; not unconditional pulses full covariance mean residual or blowup')


@lru_cache(None)
def parent_module():
    p=HERE.parent/DEPENDENCIES[0];before=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('primary_covariance_current_parent',p/'audit.py')
        module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
        return module
    finally:sys.path[:]=before


@lru_cache(None)
def parent_review():
    before=sys.path[:]
    try:
        module=parent_module();r=module.generate();module.validate(r)
        frozen=json.loads((HERE.parent/DEPENDENCIES[0]/'evidence.json').read_text())
        module.validate(frozen);module.compare_records(r,frozen)
    finally:sys.path[:]=before
    history=dict(r['inherited_provenance'])
    for name in DEPENDENCIES:
        record=json.loads((HERE.parent/name/'evidence.json').read_text())
        for path,digest in record['source_hashes'].items():
            history[name+':direct:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==digest
    return r,history


def sources():
    paths=[HERE/name for name in ('README.md','audit.py','covariance_bounds.py','covariance_checks.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/name/'evidence.json' for name in DEPENDENCIES]
    paths += [HERE.parent/name for name in ('geometric_pulse_replay/estimates.py','auxiliary_geometry_audit/geometry.py','physical_pulse_foundation/pulse_bounds.py')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def reviews():
    return deepcopy(dict(exact_perturbed_cone_and_sharp_obstruction=corner_review(),
        target_freezing_length_and_Haar_radius_scales=scale_review(),
        independent_positive_root_and_signed_Taylor_jets=jet_review(),
        manufactured_homogeneous_Gaussian_Haar_columns=gaussian_review(),
        actual_weighted_outer_pole_and_half_edge_limits=edge_review()))


def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    parent,history=parent_review();r=reviews();g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    g.update(current_inverse_and_retained_input_parent_reproduced=all(parent['gates'].values()),historical_bytes_preserved=bool(history) and all(history.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=dict(history),source_hashes=sources(),remaining_obligations=REMAINING[:])


def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:raise ValueError('Primary-root scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Remaining core/frame or nonlinear obligations erased')
    expected=reviews();keys=set(expected)|{'current_inverse_and_retained_input_parent_reproduced','historical_bytes_preserved'}
    if set(r.get('gates',{}))!=keys or not all(r['gates'].values()):raise ValueError('Primary-root gate failed')
    compare_records(expected,r.get('reviews',{}),'reviews')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if r.get('inherited_provenance')!=parent_review()[1] or not all(r['inherited_provenance'].values()):raise ValueError('History changed')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--verify-record',type=Path);args=parser.parse_args(argv)
    if args.output.exists():raise SystemExit('Use a fresh output; frozen evidence is not overwritten')
    r=generate();validate(r)
    if args.verify_record:
        frozen=json.loads(args.verify_record.read_text());validate(frozen);compare_records(r,frozen)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} bindings, {len(r['inherited_provenance'])} inherited checks pass")
    print('Positive primary roots and weighted jets from reviewed inputs; exact core/frame and full nonlinear residual obligations remain.')
if __name__=='__main__':main()
