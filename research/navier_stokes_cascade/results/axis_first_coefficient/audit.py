#!/usr/bin/env python3
"""Reproduce the analytic-axis and whole-inner first-coefficient jet bound."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from axis_bounds import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='complex-core-first-coefficient-axis-jets-and-conditional-whole-inner-collar-join'
DEPENDENCIES=('core_all_order_foundation','collar_operator_envelope','collar_third_coefficient',
    'higher_pressure_source_moments','cutoff_higher_jet_interface','first_correction_construction')
FLAGS=dict(explicit_complex_seed_core_operator_envelope_without_Phi_or_F_inverse_derived=True,
    regular_complex_ray_sparse_Picard_and_reflection_descend_first_four_components_to_X=True,
    actual_unchanged_core_first_coefficient_axis_jets_follow_by_uniqueness_and_Cauchy=True,
    explicit_pre_cutoff_whole_inner_X3_Ceta2_bound_joins_core_and_activated_collar=True,
    pressure_source_two_shifted_viscosities_and_distinct_core_collar_constants_retained=True,
    unchanged_cutoff_product_retains_third_cutoff_jet_and_global_repair_pressure_memory=True,
    historical_evidence_gauges_widths_roots_supports_frequency_and_schedule_preserved=True,
    radial_analyticity_across_flat_activation_or_cutoff_asserted=False,
    full_finalized_global_first_coefficient_angular_norms_quantified=False,
    actual_second_inner_J3_H23_or_stronger_same_schedule_budget_verified=False,
    E1_or_its_inverse_is_an_axis_regular_X_profile=False,
    actual_cutoff_scale_acceptance_and_all_physical_pulse_derivative_support_bounds_closed=False,
    complete_covariance_axial_pressure_means_and_residual_improvement_proved=False,
    independent_unconditional_core_Picard_regional_foundations_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['supply_unchanged_cutoff_jets_and_actual_global_repair_angular_inputs_with_complete_future_pressure_sources',
    'derive_second_inner_system_and_third_jet_with_all_lower_viscosity_sources_then_join_to_H23',
    'verify_retained_value_scale_acceptance_all_physical_potential_curl_pulse_jets_and_supports',
    'complete_full_wave_covariance_axial_velocity_pressure_means_and_viscous_cutoff_interaction_residual',
    'independently_check_inherited_infinite_core_and_sparse_Picard_regional_hypotheses']
def protocol():
    return dict(parent_commit='d93abfc5528b2b6b4826144973ba5aeed574eed9',
        core='analytic SEED |Lambda X|<=5;actual unchanged restriction only real X<=Xa=4/Lambda;delta=rho/(800S);no Phi or F reciprocal needed in literal first system',
        envelope='B=2(256S+9+S/Lambda)(1+Lambda)^2,H=B(1+delta^-1+2delta^-2);same complete complex absolute row/source bounds with x=5/Lambda,b=1+x',
        inverse='Gg=xi int_0^1 t^D0 g(t xi)dt;integer D0=(0,0,2,0,3,1);kernel norm<=1;no inverse xi estimate at axis',
        convergence='same sparse variable-kernel word argument on complex rays;z=2Cb,beta*=6z^2/delta;Mcore=(z+beta*)exp(beta*);normal convergence and uniqueness on smaller eta domains',
        parity='P0=diag(1,1,1,1,-1,-1);A_j(-xi)=-P0 A_j(xi)P0,f(-xi)=-P0 f(xi);first four even and last two odd;only even components descend to holomorphic X',
        Cauchy='core X disks radius1/(2Lambda),4.5<5;radial factors r!(2Lambda)^r;factorial angular C2 factor1+2/delta+4/delta^2',
        join='a=min(1,Xa),a^2<=Xa;positive scalar p_r Mcollar+q_r plus core r!(2Lambda)^r Mcore;separate Mcore,Mcollar,1 triples;pre-cutoff entire inner interval only',
        cutoff='same eta-independent kappa;binomial product keeps actual kappa third jet;no global exact repair angular bound supplied or future pressure reset',
        diagnostics='all structural reflection/kernel identities;12 manufactured analytic six-system orders and72 scalar recurrences;three original floating structure samples;not selected PDE coefficients or finite truncation proof',
        scope='conditional on retained infinite-core/Picard inputs;no radial analyticity across flat joins;actual cutoff,global repair,second H23,all-scale pulses,covariance means,residual and independent foundations remain open')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';prior=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('axis_first_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=prior
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','axis_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','core_all_order_foundation/foundation.py',
        'collar_third_coefficient/third_bounds.py','first_correction_construction/system.py',
        'first_correction_construction/majorant.py','physical_residual_budget/residual.py','cutoff_higher_jet_interface/jet_bounds.py')]
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
@lru_cache(None)
def reviews():return review()
def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=deepcopy(reviews());h=history();g=dict(r['checks']);g['historical_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Axis first-coefficient record changed')
    if not all(expected['gates'].values()):raise ValueError('Axis first-coefficient check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual cutoff, global repair and second H23 remain unverified.')
if __name__=='__main__':main()
