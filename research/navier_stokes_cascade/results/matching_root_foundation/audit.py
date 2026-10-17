#!/usr/bin/env python3
"""Same-root annular and compact five-moment matching foundation."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from matching_bounds import HERE,load,annular_review,target_review,repair_review,repair_budget_review,pressure_review
from matching_checks import reviews as exact_reviews
PROJECT=HERE.parent.parent
STATUS='actual-annular-and-compact-matching-roots-audited-relative-to-retained-core-and-compact-input-estimates'
FLAGS=dict(actual_five_row_maps_traced=True,annular_inverse_below_610=True,repair_inverse_below_800=True,
    annular_quartic_eta_coefficient_derivatives_retained=True,repair_lambda_inverse_target_loss_retained=True,
    full_sum_norm_distinguished_from_old_max_norm=True,same_roots_parameters_frequency_supports_gauge_and_width=True,
    all_five_moments_and_complete_pressure_memory_retained=True,all_fixed_one_sided_implicit_jets_supplied=True,
    historical_frozen_evidence_preserved_and_replayed=True,
    quadrature_diagnostics_promoted_to_continuum_proof=False,selected_huge_coordinates_evaluated=False,
    complete_core_global_backward_pressure_identification_discharged=False,
    all_finite_background_and_pulse_lemmas_discharged=False,unconditional_complete_residual_or_blowup_verified=False,formal_certificate=False)
REMAINING=['complete_current_core_global_backward_pressure_identification_with_same_roots_and_gauge',
    'discharge_remaining_finite_background_and_actual_pulse_existence_inverse_gap_size_derivative_support_lemmas',
    'complete_wave_covariance_mean_and_full_viscous_cutoff_interaction_residual_before_unconditional_claims']
DEPENDENCIES=('axial_root_foundation','axis_matching_audit','shape_join_pressure','annular_exterior_foundation',
    'modulation_repair_all_order','stress_realization_audit','loop_modulation_bounds','repair_state_bounds')
def protocol():return dict(parent_commit='655f9510308d39d001415d176bf259e085baf638',
    parameters='same Md64,T=exp(64)+10,lambda=exp(-4T),h=lambda^2,Tw=240T,Tf=1000,co=.001,width C^-20',
    annular='P>=16; actual five normalized rows; beta610,q50; full factorial C2 target<.286epsilon, epsilon=1e-16; same root radius348.92epsilon',
    repair='0<lambda<=.01; same original e in deltaU and deltaE; beta800,q25000; full sum C1 target<=D/N, root<=1600D/N; unchanged N',
    jets='annular triple convolution retains quartic Q_eta jets; repair Q constant only after complete original-e normalization',
    pressure='same Pi_axis+Cp; all M,J,I,S,Cp rows and fixed eta jets restored; no pressure reset',
    diagnostics='exact Catalan versus implicit polynomial series through degree8; actual existing quadrature maps on manufactured finite inputs',
    scope='matching-root lemmas relative to retained core/compact-input estimates; global backward pressure identification and finite-pulse lemmas remain')
@lru_cache(None)
def parent_module():
    p=HERE.parent/'axial_root_foundation';sys.path.insert(0,str(p))
    try:return load('matching_root_axial_parent',p/'audit.py')
    finally:sys.path.pop(0)
@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    f=json.loads((HERE.parent/'axial_root_foundation/evidence.json').read_text());m.validate(f);m.compare_records(r,f)
    history=dict(r['inherited_provenance'])
    for p,d in f['source_hashes'].items():history['axial_root_foundation:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==d
    return r,history
def sources():
    paths=[HERE/n for n in ('README.md','matching_bounds.py','matching_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/'axis_matching_audit/moments.py',HERE.parent/'stress_realization_audit/repair.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
def reviews():
    r=dict(actual_annular_inverse_and_quartic_bound=annular_review(),actual_annular_target=target_review(),
        actual_compact_repair_inverse=repair_review(),actual_compact_target_and_frequency=repair_budget_review(),
        complete_pressure_and_jet_scope=pressure_review());r.update(exact_reviews());return deepcopy(r)
def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    p,h=parent_review();r=reviews();g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    g.update(current_axial_parent_reproduced=all(p['gates'].values()),historical_bytes_preserved=bool(h) and all(h.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,inherited_provenance=dict(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)
def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:raise ValueError('Matching scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Remaining obligations erased')
    if set(r.get('gates',{}))!=set(reviews())|{'current_axial_parent_reproduced','historical_bytes_preserved'} or not all(r['gates'].values()):raise ValueError('Matching gate failed')
    compare_records(reviews(),r.get('reviews',{}),'root.reviews')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=1069 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; historical evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:
        f=json.loads(a.verify_record.read_text());validate(f);compare_records(r,f)
    a.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} bindings, {len(r['inherited_provenance'])} inherited checks pass")
    print('Same matching roots supplied relative to retained inputs; global core pressure and finite-pulse obligations remain.')
if __name__=='__main__':main()
