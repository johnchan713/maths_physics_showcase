#!/usr/bin/env python3
"""Conditional finite-energy competitor comparison for the current localized field."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from comparison import load,exponent_review,limit_review,interface_review
from polynomial_checks import forced_difference_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='conditional-current-field-finite-energy-comparison-and-global-competitor-exclusion'
FLAGS=dict(same_current_localized_field_force_zero_datum_and_full_pressure_are_input=True,
    general_comparison_hypotheses_checked_on_every_fixed_preterminal_interval=True,
    pressure_gradient_recovered_without_individual_unforced_pressure_replacement=True,
    all_infinity_fluxes_controlled_without_global_competitor_derivative_bounds=True,
    conditional_current_global_smooth_bounded_energy_competitor_exclusion_supplied=True,
    inherited_foundation_is_explicitly_still_required=True,
    diagnostics_not_actual_fields_or_a_continuum_proof=True,
    historical_frozen_evidence_and_tolerances_unchanged=True,
    independent_analytic_foundation_verified=False,
    uniform_comparison_constants_through_terminal_time_assumed=False,
    future_velocity_or_pressure_constructed=False,
    unconditional_admissible_construction_or_exclusion_verified=False,
    unforced_singularity_verified=False,blowup_verified=False,formal_certificate=False)
REMAINING=['independently_trace_the_current_leading_profile_background_inverse_and_all_jet_physical_bounds_to_defining_maps',
    'identify_and_prove_the_first_unclosed_continuum_estimate_or_document_its_precise_obstruction',
    'unconditional_theorem_assessment_only_after_every_foundation_premise_is_verified']


def protocol():
    return dict(parent_commit='236056afbfd2d20e1839733708fdc53371835bc4',
        current_input='current_localized_force; never substitute the historical localized field',
        general_lemma='global_solution_comparison sections 1-6, reproduced with its own comparator',
        viscosity='1',cutoff_power=8,negative_Sobolev_orders=[2,3],
        commutator_Lp='4/3',commutator_radius_power='-3/4',energy_radius_decay='4/3',
        order_of_limits='fix T<1, then R->infinity; only afterward approach t=1 along protected core',
        external_input='Riesz Lp boundedness, Knapp Advanced Real Analysis III.5-6, Theorem 3.26 and Riesz example printed pp83,93 (PDF pp113,123), checked 2026-10-16',
        source='https://www.math.stonybrook.edu/~aknapp/download/a2-1-realanal-clickable.pdf',
        scope='conditional current competitor exclusion, not an independent construction or blow-up theorem')


@lru_cache(None)
def parent_module():
    p=HERE.parent/'current_localized_force';sys.path.insert(0,str(p))
    try:return load('current_comparison_force_parent',p/'audit.py')
    finally:sys.path.pop(0)


@lru_cache(None)
def parent_review():
    p=parent_module();r=p.generate();p.validate(r)
    frozen=json.loads((HERE.parent/'current_localized_force/evidence.json').read_text())
    p.validate(frozen);p.compare_records(r,frozen)
    inherited=dict(r['inherited_provenance'])
    for path,digest in frozen['source_hashes'].items():
        inherited['current_localized_force:'+path]=hashlib.sha256((PROJECT/path).read_bytes()).hexdigest()==digest
    return r,inherited


@lru_cache(None)
def general_review():
    p=HERE.parent/'global_solution_comparison';sys.path.insert(0,str(p))
    try:m=load('current_comparison_general_lemma',p/'audit.py')
    finally:sys.path.pop(0)
    r=m.generate();m.validate(r)
    frozen=json.loads((p/'evidence.json').read_text());m.validate(frozen);m.compare_records(r,frozen)
    return dict(historical_field_not_substituted_for_current=True,
        controls={k:{n:bool(v) for n,v in x['checks'].items()} for k,x in r['reviews'].items()},
        checks=dict(all_general_comparison_controls_reproduced=all(r['gates'].values()),
            every_historical_inherited_hash_agrees=all(r['inherited_provenance'].values()),
            all_general_lemma_source_bytes_unchanged=r['source_hashes']==frozen['source_hashes']))


def sources():
    names=('README.md','comparison.py','polynomial_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')
    paths=[HERE/n for n in names]+[HERE.parent/n/'evidence.json' for n in ('current_localized_force','global_solution_comparison')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def reviews():
    return dict(current_hypothesis_interface=interface_review(parent_review()[0]),
        general_pressure_and_infinity_controls=deepcopy(general_review()),
        exact_forced_difference_pressure=forced_difference_review(),
        independent_weight_kernel_and_absorption_exponents=exponent_review(),
        fixed_interval_and_protected_core_limits=limit_review())


def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    _,inherited=parent_review();r=reviews()
    gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    gates['current_parent_reproduced']=all(parent_review()[0]['gates'].values())
    gates['historical_provenance_preserved']=bool(inherited) and all(inherited.values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(inherited),source_hashes=sources(),remaining_obligations=REMAINING[:])


def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)


def validate(r):
    if r.get('status')!=STATUS or r.get('flags')!=FLAGS or r.get('protocol')!=protocol():raise ValueError('Conditional scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Foundation obligation erased')
    expected=reviews();keys=set(expected)|{'current_parent_reproduced','historical_provenance_preserved'}
    if set(r.get('gates',{}))!=keys or not all(r['gates'].values()):raise ValueError('Comparison gate missing or failed')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=748 or not all(r['inherited_provenance'].values()):raise ValueError('Provenance changed')
    compare_records(expected,r.get('reviews',{}),'root.reviews')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    args=p.parse_args(argv)
    if args.output.exists():raise SystemExit('Choose a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if args.verify_record:
        frozen=json.loads(args.verify_record.read_text());validate(frozen);compare_records(r,frozen)
    args.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['inherited_provenance'])} inherited hashes pass")
    print('Conditional current competitor exclusion only; independent foundation and blow-up remain unverified.')


if __name__=='__main__':main()
