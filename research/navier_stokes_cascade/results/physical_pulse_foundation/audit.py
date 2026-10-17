#!/usr/bin/env python3
"""Actual primary physical-pulse foundation, retaining later source obligations."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from pulse_bounds import input_review,pulse_review,physical_review,support_review
from pulse_checks import diagnostic_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='actual-primary-physical-pulses-all-fixed-jets-and-common-domain-tails-relative-to-exact-core-and-finite-lemmas'
FLAGS=dict(actual_background_inputs_joined_to_primary_pulse_frame_and_covariance=True,
    one_starting_band_for_all_later_scales_and_every_fixed_physical_derivative_order=True,
    uniform_primary_complete_curl_pressure_and_fixed_potential_bounds_derived=True,
    uncut_radial_ODE_pulse_has_positive_lower_bound_relative_to_P=True,
    exact_physical_graph_frame_and_moving_pressure_terms_retained=True,
    same_domain_Gaussian_tail_flatness_constants_derived_for_every_fixed_order=True,
    exact_closed_supports_and_one_sided_terminal_jets_preserved=True,
    complete_background_radial_flow_and_I4_pressure_tails_retained=True,
    exact_core_and_finite_existence_inverse_gap_lemmas_remain_dependencies=True,
    historical_sources_evidence_parameters_width_frequency_and_supports_unchanged=True,
    finite_manufactured_controls_not_continuum_suprema=True,
    positive_pointwise_lower_bound_for_oscillating_cutoff_velocity=False,
    generated_nonlinear_source_bounds_and_supports_all_verified=False,
    unconditional_all_actual_pulse_families_verified=False,
    unconditional_complete_covariance_mean_corrections_verified=False,
    unconditional_complete_residual_improvement_verified=False,
    actual_giant_parameters_fields_or_starting_band_numerically_evaluated=False,
    blowup_verified=False,formal_certificate=False)
REMAINING=['independently_audit_complete_wave_covariance_means_with_axial_radial_velocity_and_Pi1_Pi2',
    'verify_bounds_and_full_harmonic_endpoint_supports_for_every_actual_generated_nonlinear_source',
    'verify_complete_residual_viscosity_all_cutoffs_and_interactions_before_unconditional_claims']
DEPENDENCIES=('weighted_stress_foundation','physical_pulse_transfer','actual_band_pulse_bounds',
    'geometric_pulse_replay','auxiliary_geometry_audit','physical_wave_interface','background_jet_foundation','coefficient_support_audit')


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m


def protocol():
    return dict(parent_commit='585342bc4225cd69cf07a0a9a24d65ef01729804',
        unchanged_parameters='Md64; unit viscosity; fixed h,N,t1=C^-20; r0=1/6057072',
        enlarged_chart='1/4<=q/Q<=4; T>=0; 4Q<qstar; representatives in closed active intersection',
        physical_losses='sx=1/2+h/2; st=1+3h/2',
        primary_powers='velocity Q^(-(a+1)*sx-b*st); pressure Q^(-(a+2)*sx-b*st); fixed potential Q^(-a*sx-b*st)',
        tail_same_start='D=B+max(0,Lambda+K); exp(D*ell-a0*ell^2)<=exp(D^2/(4*a0))',
        source_mapping_not_generated_source_verification=True,
        discrete_data_frozen_before_all_differentiation=True,
        scope='uniform actual PRIMARY pulse foundation relative to exact core and finite lemmas; complete means generated sources full residual and blow-up not promoted')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'weighted_stress_foundation';sys.path.insert(0,str(p))
    try:return load('physical_pulse_foundation_parent',p/'audit.py')
    finally:sys.path.pop(0)


@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    f=json.loads((HERE.parent/'weighted_stress_foundation/evidence.json').read_text());m.validate(f);m.compare_records(r,f)
    inherited=dict(r['inherited_provenance'])
    for p,digest in f['source_hashes'].items():
        inherited['weighted_stress_foundation:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==digest
    return r,inherited


def sources():
    paths=[HERE/n for n in ('README.md','pulse_bounds.py','pulse_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/'core_all_order_foundation/coordinate_checks.py',HERE.parent/'physical_pulse_transfer/transfer.py',HERE.parent/'geometric_pulse_replay/estimates.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def reviews():return deepcopy(dict(actual_background_and_common_band_inputs=input_review(),
    moving_frame_positive_pulse_and_Haar_weights=pulse_review(),
    all_fixed_physical_curl_pressure_potential_jets=physical_review(),
    exact_supports_and_same_domain_Gaussian_tails=support_review(),
    independent_exact_physical_frame_and_tail_controls=diagnostic_review()))


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    p,inherited=parent_review();r=reviews()
    gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    gates.update(actual_background_parent_and_retained_chain_reproduced=all(p['gates'].values()),
        historical_bytes_preserved=bool(inherited) and all(inherited.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(inherited),source_hashes=sources(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('flags')!=FLAGS or r.get('protocol')!=protocol():raise ValueError('Primary pulse scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Mean/source/full residual gates erased')
    expected=reviews();keys=set(expected)|{'actual_background_parent_and_retained_chain_reproduced','historical_bytes_preserved'}
    if set(r.get('gates',{}))!=keys or not all(r['gates'].values()):raise ValueError('A primary pulse gate failed')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=855 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
    compare_records(expected,r.get('reviews',{}),'root.reviews')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    args=p.parse_args(argv)
    if args.output.exists():raise SystemExit('Use a fresh output; historical evidence is never overwritten')
    r=generate();validate(r)
    if args.verify_record:
        f=json.loads(args.verify_record.read_text());validate(f);compare_records(r,f)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} bindings, {len(r['inherited_provenance'])} inherited checks pass")
    print('Primary bounds relative to retained lemmas; generated sources, complete means/residual and blow-up remain open.')


if __name__=='__main__':main()
