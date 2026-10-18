#!/usr/bin/env python3
"""Leading pressure identification and quantitative regular angular recovery."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from pressure_recovery import load,inverse_review,obstruction_review,scale_and_memory_review,input_review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='leading-core-pressure-identified-with-regular-angular-inverse-and-retained-higher-memory'
DEPENDENCIES=('canonical_comparison_join','core_growth_activation','core_contraction_audit',
    'core_all_order_foundation','global_pressure_moment_foundation','annular_exterior_foundation',
    'coefficient_support_foundation')
FLAGS=dict(reference_datum_independent_of_radial_attachment_scale=True,
    exact_positive_measure_not_diagnostic_quadrature_used_for_core=True,
    leading_pressure_identification_already_supplied_by_same_root_ledger=True,
    core_datum_joined_to_same_completed_leading_backward_pressure_relative_to_accepted_maps=True,
    bounded_angular_inverse_and_weighted_endpoint_compatibility_derived=True,
    every_fixed_angular_defect_jet_has_explicit_majorant=True,
    singular_homogeneous_mode_detects_essential_endpoint_condition=True,
    higher_background_Pi1_Pi2_and_complete_canonical_pressure_memory_retained=True,
    next_actual_canonical_value_acceptance_inequality_named=True,
    historical_sources_fields_parameters_widths_and_evidence_unchanged=True,
    finite_rational_fixtures_distinguished_from_continuum_argument=True,
    full_corrected_NS_pressure_identified_with_leading_centrifugal_pressure=False,
    diagnostic_pressure_measure_identified_with_actual_reference_measure=False,
    leading_pressure_identification_still_an_independent_missing_equation=False,
    bounded_endpoint_condition_omitted=False,
    heat_endpoint_complex_analyticity_assumed=False,
    all_jet_constants_uniform_in_order_or_h=False,
    actual_core_or_giant_schedule_numerically_evaluated=False,
    actual_canonical_first_band_value_acceptance_verified=False,
    every_retained_core_background_and_finite_analytical_premise_discharged=False,
    unconditional_pulse_complete_residual_or_blowup_proved=False,formal_certificate=False)
REMAINING=['audit_common_canonical_root_value_acceptance_theta_le_one_over_two_kappa_plus_one_and_normal_Riccati_cone_thresholds_on_actual_unchanged_accepted_domain_not_manufactured_inputs',
    'discharge_remaining_actual_defining_existence_inverse_gap_and_full_construction_premises_before_unconditional_pulse_residual_or_blowup_claims']

def protocol():
    return dict(parent_commit='0b16f1b18aa6277dc507e39fd7cecbd9a6e4016a',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        fixed='same Md64,P,h,Lambda,C,N,positive lambda,width C_setup^-20,leading roots,pulse labels,supports,qstar and gauges',
        reference='Pi_axis=-1/2 int E_id,sched(y,eta)^2 dy; y=log(X/XR); dX/X=dy cancels radial attachment scale',
        family='exact positive pushforward measure on [0,1], mass<5, inner atom>=5/2 at theta=1; diagnostic quadrature not substituted',
        leading='same-root leading ledger Pi_axis+Cp_infinity=0 and all fixed eta jets; leading physical field not rebuilt by canonical wave-root change',
        defect='delta=Pi_forward-Pi_backward=Pi_axis+Cp_infinity; Delta(r*Tz)=X*K_A(delta)/L',
        operator='K_A delta=4A eta delta-(1-eta^2)delta_prime; A=1/2+h,0<h<=1/100',
        inverse='p=2A; (d^p delta)_prime=-d^(p-1)g; weighted integral g=0; right and left endpoint integrals recover unique smooth bounded delta',
        value_bound='||delta||_infinity<=||K_A delta||_infinity<=||Gamma||_infinity, Gamma=(K_A delta)/L, L=1-2h eta^2',
        jets='K_(A+k/2) delta^(k)=g^(k)-k(4A+k-1)delta^(k-1); B0=G0,Bk=2^k*(Gk+k*(4A+k-1)*B_(k-1))',
        source_jets='Gk=||Gamma^k||+4hk||Gamma^(k-1)||+2hk(k-1)||Gamma^(k-2)||; negative indices zero',
        endpoint='nonzero delta=d^(-2A) has zero operator on open interval but nonvanishing weighted endpoint trace; g=1 fails compatibility',
        pressure_memory='on I4, q^(2A)p_B=Pi0+chi(c1q)q^(2h)(K1-mu^2/(4X))+chi(c2q)q^(4h)K2; no full-pressure collapse',
        next_gate='theta=12CH/(ja_minus*min(Amin,u)*L)<=1/(2*kappa+1),kappa=24ja_plus/(gamma*ja_minus),plus common normal/Riccati/cone thresholds; actual domain acceptance not proved here',
        scope='leading datum compatibility relative to same accepted leading maps; new regular angular inverse is analytical; no unconditional full construction or blowup')

@lru_cache(None)
def parent_module():
    p=HERE.parent/'canonical_comparison_join';old=sys.path[:];sys.path.insert(0,str(p))
    try:return load('core_pressure_canonical_comparison_parent',p/'audit.py')
    finally:sys.path[:]=old

@lru_cache(None)
def records():
    r={n:json.loads((HERE.parent/n/'evidence.json').read_text()) for n in DEPENDENCIES}
    parent_module().validate(r['canonical_comparison_join'])
    return r

def history():return parent_module().provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])

def sources():
    paths=[HERE/n for n in ('README.md','audit.py','pressure_recovery.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    paths += [HERE.parent/'global_pressure_moment_foundation/identities.py',HERE.parent/'outer_pressure_pilot/schedule.py']
    paths += [HERE.parent/'canonical_primary_foundation'/f for f in ('README.md','protocol.json')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def reviews():
    return deepcopy(dict(same_reference_core_and_leading_pressure_trace=input_review(records()),
        exact_endpoint_inverse_and_fixed_jet_controls=inverse_review(),
        endpoint_obstructions_and_scope=obstruction_review(),
        radial_scale_stress_defect_and_higher_memory=scale_and_memory_review()))

def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)

def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=reviews();h=history();g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    g['historical_recursive_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=h,source_hashes=sources(),remaining_obligations=REMAINING[:])

def validate(r):
    if r.get('status')!=STATUS:raise ValueError('Pressure-identification scope changed')
    for a,b,path in ((protocol(),r.get('protocol'),'protocol'),(FLAGS,r.get('flags'),'flags'),
        (REMAINING,r.get('remaining_obligations'),'remaining'),(reviews(),r.get('reviews'),'reviews'),
        (history(),r.get('inherited_provenance'),'history'),(sources(),r.get('source_hashes'),'sources')):
        compare_records(a,b,path)
    g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in reviews().items()}
    g['historical_recursive_source_bytes']=bool(history()['bindings']) and all(history()['bindings'].values())
    compare_records(g,r.get('gates'),'gates')
    if not all(g.values()):raise ValueError('Pressure-identification gate failed')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:
        f=json.loads(a.verify_record.read_text());validate(f);compare_records(r,f)
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} direct bindings, {len(r['inherited_provenance']['bindings'])} historical bindings pass")
    print('Leading compatibility and bounded angular recovery only; actual value acceptance and unconditional proof remain open.')
if __name__=='__main__':main()
