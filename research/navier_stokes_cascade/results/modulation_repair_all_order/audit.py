#!/usr/bin/env python3
"""All fixed orders of the SAME finite-N modulation, phase and exact repair maps."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from all_order_maps import map_review,transfer_review
from map_checks import diagnostic_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='same-finite-N-modulation-moving-phase-and-exact-repair-all-order-map-estimate'
FLAGS=dict(all_fixed_orders_of_signed_variance_and_moving_phase_inverse_derived=True,
    same_accepted_C1_repair_has_all_fixed_derivatives_by_pointwise_Jacobian=True,
    exact_five_row_full_pressure_memory_and_normalization_losses_retained=True,
    field_and_moment_radial_frequency_powers_distinguished=True,
    finite_N_and_all_supports_and_historical_evidence_unchanged=True,
    low_order_cone_acceptance_identifies_same_map_conditionally=True,
    smooth_input_on_whole_compact_region_remains_explicit=True,
    independent_foundation_parent_covers_inner_shape_not_whole_exterior=True,
    diagnostics_are_manufactured_exact_jets_not_actual_continuum_certificates=True,
    simultaneous_Ck_contractions_or_frequency_reselection_required=False,
    modulated_compact_profile_globally_complex_eta_analytic=False,
    all_radial_derivative_errors_small_or_N_uniform=False,
    whole_annular_exterior_and_global_profile_foundation_independently_verified=False,
    unconditional_actual_pulse_bounds_verified=False,
    unconditional_complete_residual_improvement_verified=False,
    blowup_verified=False,formal_certificate=False)
REMAINING=['independently_supply_remaining_annular_exterior_pre_modulation_all_order_inputs',
    'complete_exact_global_pressure_moment_identities_weighted_edges_and_reserved_support_join',
    'only_then_promote_uniform_actual_pulse_mean_and_complete_residual_claims_beyond_conditionals']
DEPENDENCIES=('prescribed_continuation_foundation','loop_modulation_bounds','repair_state_bounds',
    'higher_jet_budget','stress_realization_audit','compact_jet_envelope','regional_cone_join')


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m


def protocol():
    return dict(parent_commit='9f4ae8973fdd1384530085505010597033e206cf',
        unchanged_frequency='N=1+floor(H^32)',unchanged_repair='same C1 small root of exact Bc+Q(c,c)=z_N',
        input='C-infinity pre-modulation profile on WHOLE compact modulation-through-repair rectangle with positive floors',
        root='signed F=d0*mu*sqrt(R2(mu*p)), not unsigned sqrt(V)',
        field_frequency_power='N^(r-1)',moment_pressure_frequency_power='N^max(r-2,-1)',
        repair_angular_frequency_power='N^-1 for every fixed k via the pointwise Jacobian',
        no_parameter_support_or_frozen_evidence_change=True,
        scope='all-order MAP theorem relative to smooth input; annular/exterior independent input audit remains open')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'prescribed_continuation_foundation';sys.path.insert(0,str(p))
    try:return load('modulation_all_order_continuation_parent',p/'audit.py')
    finally:sys.path.pop(0)


@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    f=json.loads((HERE.parent/'prescribed_continuation_foundation/evidence.json').read_text());m.validate(f);m.compare_records(r,f)
    inherited=dict(r['inherited_provenance'])
    for p,digest in f['source_hashes'].items():
        inherited['prescribed_continuation_foundation:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==digest
    return r,inherited


def sources():
    paths=[HERE/n for n in ('README.md','all_order_maps.py','map_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/'core_all_order_foundation/coordinate_checks.py',HERE.parent/'stress_realization_audit/loop.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def reviews():return deepcopy(dict(implicit_root_phase_and_same_repair_maps=map_review(),
    frequency_powers_full_memory_and_input_limits=transfer_review(),independent_exact_higher_jet_controls=diagnostic_review()))


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    parent,inherited=parent_review();r=reviews()
    gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    gates.update(parent_foundation_and_retained_conditional_chain_reproduced=all(parent['gates'].values()),
        historical_source_bytes_preserved=bool(inherited) and all(inherited.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(inherited),source_hashes=sources(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('flags')!=FLAGS or r.get('protocol')!=protocol():raise ValueError('All-order map scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Global input gate erased')
    expected=reviews();keys=set(expected)|{'parent_foundation_and_retained_conditional_chain_reproduced','historical_source_bytes_preserved'}
    if set(r.get('gates',{}))!=keys or not all(r['gates'].values()):raise ValueError('An all-order map gate failed')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=801 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
    compare_records(expected,r.get('reviews',{}),'root.reviews')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    args=p.parse_args(argv)
    if args.output.exists():raise SystemExit('Use a fresh output; historical frozen evidence is never overwritten')
    r=generate();validate(r)
    if args.verify_record:
        f=json.loads(args.verify_record.read_text());validate(f);compare_records(r,f)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} source/record bindings and {len(r['inherited_provenance'])} inherited checks pass")
    print('Relative all-order MAP theorem only; remaining annular/exterior input, global profile, pulses, residual and blow-up stay open.')


if __name__=='__main__':main()
