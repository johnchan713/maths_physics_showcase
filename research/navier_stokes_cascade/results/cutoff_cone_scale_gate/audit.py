#!/usr/bin/env python3
"""Reproduce corrected-width affine acceptance bounds, not an accepted actual scale."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from gate_bounds import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='corrected-collar-width-and-max-affine-cutoff-cone-scale-gate'
DEPENDENCIES=('global_swirl_envelope','global_axial_source_envelope','cutoff_repair_derivative_budget',
    'first_seminorm_bound','leading_cone_margin','background_gap_certificate','regional_weighted_gate','core_growth_activation')
FLAGS=dict(historical_tc_equals_native_t1_identification_corrected_in_new_checkpoint=True,
    actual_tc_is_t1_over_ten_sqrt_logC_and_native_width_unchanged=True,
    complete_first_stress_bound_is_exact_max_of_four_nonnegative_affine_shape_rows=True,
    full_future_signed_pressure_moments_and_separate_pointwise_source_costs_retained=True,
    positive_viscosity_L2_coefficient_quantifies_fixed_budget_certificate_obstruction=True,
    finite_later_band_witness_conditional_on_supplied_shape_bounds_and_same_tilt_margin=True,
    original_profiles_roots_supports_cutoffs_frequency_gauges_and_mean_memories_unchanged=True,
    historical_frozen_evidence_and_original_comparators_preserved=True,
    historical_cutoff_width_misidentification_endorsed_as_true=False,
    actual_cutoff_shape_existing_tilt_or_retained_scale_acceptance_verified=False,
    failure_of_sufficient_upper_bound_is_actual_cone_failure=False,
    compressed_auxiliary_shape_replaces_original_cutoff=False,
    actual_huge_profile_or_physical_pulses_numerically_resolved=False,
    retained_core_and_regional_premises_independently_proved_unconditional=False,
    unconditional_pulse_complete_residual_improvement_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['supply_valid_enclosures_for_unchanged_coefficient_cutoff_L1_L2_with_correct_tc_not_native_t1',
    'verify_same_existing_tilt_margin_or_sharpen_it_then_check_four_affine_rows_at_the_retained_scale_and_all_other_acceptance_gates',
    'prove_uniform_actual_physical_pulse_size_derivative_inverse_gap_and_support_on_every_scale',
    'retain_complete_wave_covariance_axial_velocity_pressure_and_means_in_full_viscous_cutoff_interaction_residual',
    'independently_verify_retained_exact_core_and_regional_foundations_before_unconditional_claims']
def protocol():
    return dict(parent_commit='4de29921a9af3fc05e5ade7badb9a0c2361919ba',
        correction='historical cutoff README/protocol tc=C^-20 is wrong: native t1=C^-20, tc=t1/(10 sqrt(log C)), Delta=tc/16',
        unchanged_geometry='Xinner=Xa exp(tc/16),Xkeep=Xa exp(tc/8),Xcut=Xa exp(3tc/16),analytic edge Xa exp(tc/4),leading modulation Xa exp(tc/2)',
        shape_jets='S=s0;Sx=x0+x1 L1;Sxx=z0+z1 L1+z2 L2; coefficients independent of shape derivatives; positive-order repair additive',
        combined_gate='Cstar=max(Ru,Rv,Rd,Rw), four nonnegative affine rows; pressure and shifted viscosities enter constant rows once',
        viscous_obstruction='dRw/dL2=nu*c^3*(c+2)*J/(tau*Delta^2*Xkeep^2)>0 for positive nu,J; actual nu=1',
        scale_region='every row<=r_Kc*2^(2h*(ell-2)); L1>=1,L2>=4 necessary not sufficient for an actual shape',
        later_witness='m=max(0,ceil(log2(Cstar/r)));ell_suf=2+ceil(m/(2h)); conditional on unchanged Borel tail and all other gates',
        compressed_shapes='Kn=1-sigma(n(s-1/2)+1/2);L1>=8n,L2>=16n^2; n=1+ceil(B/(16A2)) defeats this finite fixed certificate, not actual stress',
        diagnostics='twelve auxiliary rational complete-interface comparisons, symbolic identities and exact boundary/dyadic controls; no actual-profile samples',
        scope='corrected conditional gate and precise uniform-certificate obstruction; unchanged cutoff shape, tilt, existing scale, all physical pulses and complete residual not verified')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('combined_gate_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','gate_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','global_axial_source_envelope/envelopes.py',
        'global_swirl_envelope/swirl_bounds.py','cutoff_repair_derivative_budget/derivatives.py','first_seminorm_bound/seminorms.py',
        'leading_cone_margin/cone_bounds.py','background_gap_certificate/gap_bounds.py','signed_pressure_source_moments/moments.py',
        'annular_first_correction/README.md','axis_core_attachment/README.md')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
@lru_cache(None)
def reviews():return review()
def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=deepcopy(reviews());h=history();g=dict(r['checks'])
    g['historical_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Combined gate record changed')
    if not all(expected['gates'].values()):raise ValueError('Combined gate check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; unchanged cutoff shape, tilt and actual scale remain unverified.')
if __name__=='__main__':main()
