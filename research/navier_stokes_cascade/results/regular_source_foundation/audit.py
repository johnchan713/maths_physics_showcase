#!/usr/bin/env python3
"""Regular source foundation; no full nonlinear residual promotion."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from source_bounds import review
from source_checks import reviews as exact_reviews

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='regular-source-foundation-sharper-linear-gain-and-actual-initial-supported-input-relative-to-core-and-pulse-lemmas'
FLAGS=dict(all_fixed_coefficient_jets_used_without_repeated_graph_loss=True,
    regular_complete_linear_gain_one_half_minus_kappa_derived=True,
    complete_curl_longitudinal_identity_not_principal_transversality_used=True,
    genuine_mean_wave_half_power_loss_retained=True,
    actual_primary_plus_initial_mean_regular_source_one_minus_kappa=True,
    full_source_endpoint_extension_and_closed_label_supports_preserved=True,
    full_viscosity_pressure_graph_coefficients_and_cylindrical_connections_retained=True,
    complete_Pi1_Pi2_background_pressure_and_mean_moments_retained=True,
    physically_flat_cutoff_tails_kept_additive_not_forward_inverted=True,
    exact_core_and_finite_background_pulse_lemmas_remain_dependencies=True,
    historical_sources_evidence_parameters_supports_and_carriers_unchanged=True,
    finite_manufactured_jet_checks_not_continuum_suprema=True,
    all_actual_later_source_families_and_complex_targets_verified=False,
    sharper_full_cycle_rate_established=False,
    unconditional_all_pulse_family_bounds_verified=False,
    unconditional_complete_residual_improvement_verified=False,
    actual_giant_fields_or_first_band_numerically_evaluated=False,
    blowup_verified=False,formal_certificate=False)
REMAINING=['join_actual_first_sourced_correction_full_complex_covariance_and_recomputed_mean_pressure_defects',
    'independently_verify_later_complex_target_bounds_and_every_actual_source_family',
    'close_complete_residual_viscosity_cutoffs_all_interactions_before_unconditional_claims']
DEPENDENCIES=('compact_mean_foundation','physical_pulse_foundation','actual_band_pulse_bounds',
    'physical_wave_interface','residual_cycle_closure','first_sourced_wave_update',
    'finite_stage_induction','geometric_pulse_replay')


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m


def protocol():
    return dict(parent_commit='10b5bdff358cd6c54474be533dec1618f829d940',
        unchanged_parameters='Md64; unit viscosity; fixed h,N,positive lambda,t1=C^-20; original radius supports and carriers',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        regular_linear_gain='1/2-kappa; old 1/2-3*kappa remains conservative',
        curl_principal_order='alpha+1/2-kappa; c_i*N_t is a coefficient derivative, not another Dr',
        nonlinear_orders='same-label alpha+alpha_prime-kappa; mean-wave alpha+mu-1/2',
        complete_divergence='div_star(a)+i*k*m*n.a=0; complete amplitude is generally not transverse',
        initial_regular_source='primary W_(1/2); means M_(1-kappa), radial M_(2-kappa); source W_(1-kappa)',
        exact_tail_split='psi*regular_linear_remainder + (1-psi)*f + T_psi; last two additive, not inverse sources',
        scope='initial supported source and regular operator grammar relative to retained analytical inputs; full later complex iteration and residual remain open')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'compact_mean_foundation';sys.path.insert(0,str(p))
    try:return load('regular_source_foundation_parent',p/'audit.py')
    finally:sys.path.pop(0)


@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    f=json.loads((HERE.parent/'compact_mean_foundation/evidence.json').read_text());m.validate(f);m.compare_records(r,f)
    history=dict(r['inherited_provenance'])
    for p,digest in f['source_hashes'].items():
        history['compact_mean_foundation:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==digest
    return r,history


def sources():
    paths=[HERE/n for n in ('README.md','source_bounds.py','source_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/'actual_band_pulse_bounds/diagnostics.py',HERE.parent/'modulation_repair_all_order/map_checks.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def reviews():return deepcopy(dict(exponent_accounting=review(),**exact_reviews()))


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    p,history=parent_review();r=reviews()
    gates=dict(regular_linear_exponent_accounting=all(r['exponent_accounting']['checks'].values()),
        exact_complete_linear_identity=r['exact_full_viscous_linear_expansion']['principal_plus_remainder_matches'] and r['exact_full_viscous_linear_expansion']['all_omission_controls_nonzero'],
        exact_complete_curl_transport=all(r['exact_complete_curl_nonlinear_transport'][k] for k in ('complete_harmonic_transport_matches','exact_complete_curl_longitudinal_identity','complete_amplitude_not_transverse','zero_and_negative_outputs_present')),
        endpoint_and_flat_tail_scope=r['endpoint_and_tail_split']['cut_potential_and_pressure_endpoint_jets'] and r['endpoint_and_tail_split']['flat_tails_must_remain_additive_not_forward_sources'],
        actual_primary_and_mean_foundations_reproduced=all(p['gates'].values()),historical_bytes_preserved=bool(history) and all(history.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(history),source_hashes=sources(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:raise ValueError('Source foundation scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Later source/full residual obligations erased')
    keys={'regular_linear_exponent_accounting','exact_complete_linear_identity','exact_complete_curl_transport','endpoint_and_flat_tail_scope','actual_primary_and_mean_foundations_reproduced','historical_bytes_preserved'}
    if set(r.get('gates',{}))!=keys or not all(r['gates'].values()):raise ValueError('A source foundation gate failed')
    compare_records(reviews(),r.get('reviews',{}),'root.reviews')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=889 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
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
    print('Initial regular sources only; all later complex states, complete residual and blow-up remain open.')


if __name__=='__main__':main()
