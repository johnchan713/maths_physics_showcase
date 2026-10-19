#!/usr/bin/env python3
"""All-finite repair-bump jets and higher-cutoff input boundary; no actual acceptance."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from jet_bounds import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='all-finite-repair-bump-jets-and-higher-cutoff-input-boundary'
DEPENDENCIES=('sharp_current_scale_floor','cutoff_repair_derivative_budget',
    'first_seminorm_bound','background_jet_foundation','regional_weighted_gate')
FLAGS=dict(all_finite_jets_of_same_actual_repair_bumps_bounded_analytically=True,
    conditional_cutoff_collar_repair_and_physical_radial_jet_interface_derived=True,
    low_two_shape_norms_do_not_bound_third_jets_even_for_monotone_cutoffs=True,
    exact_repair_disjoint_support_and_full_pressure_memory_retained=True,
    viscous_second_stress_derivative_principal_cost_is_third_velocity_jet=True,
    original_all_jet_Borel_requirements_still_required_not_reported_as_missing=True,
    historical_evidence_and_construction_choices_unchanged=True,
    unchanged_cutoff_higher_jets_or_actual_nondegenerate_collar_value_verified=False,
    actual_first_band_or_all_physical_pulse_derivative_gates_closed=False,
    all_jet_constants_can_be_replaced_by_low_order_Cstar=False,
    auxiliary_cutoffs_are_the_actual_field_or_actual_residual_counterexamples=False,
    complete_covariance_axial_pressure_means_and_residual_improvement_proved=False,
    independent_unconditional_core_regional_foundations_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['supply_unchanged_cutoff_Lj_inner_collar_Jj_and_full_moment_inverse_norms_for_each_required_derivative_order',
    'verify_actual_low_order_cutoff_raw_deficit_endpoint_constants_and_full_mixed_retained_scale_gate',
    'retain_order_dependent_Borel_curl_cutoff_constants_and_verify_all_physical_scale_pulse_geometry_supports_and_jets',
    'complete_full_wave_covariance_and_axial_velocity_pressure_mean_corrections',
    'prove_complete_viscous_cutoff_interaction_residual_and_independently_check_retained_foundations']
def protocol():
    return dict(parent_commit='f84076cef9550e7bd146493ce9a5a4f7afdc5b3d',
        fixed_width='tc=C^-20/(10sqrt(log C));Delta=tc/16;no cutoff or field replaced',
        log_jets='dX^m[X^-alpha b(log X)]=X^(-m-alpha) product_(i<m)(Dy-alpha-i)b;all Jacobian terms retained',
        repaired_jets='km=Xkeep^-m sum_j |s(m,j)|Delta^-j Lj;Pm=sum binom(m,j)kj Jm-j;Sm=Pm+Cp max(Um,Fm);disjoint exact repair supports',
        actual_bumps='same beta=10sigma_prime(10(y-c)+1/2),squared-flat sigma;Bell/logistic recurrence gives explicit finite Dm;Bm=10^(m+1)D_(m+1);old first three bounds retained',
        physical_radial='fixed physical z,t;X=r^2/(2q);d_r^m(Hf)=H q^-m/2 sum_k m!R^(2k-m)/(2^(m-k)(m-k)!(2k-m)!) dX^k f;not all Cartesian/axial/time/potential jets',
        auxiliary_control='Kn=1-sigma+a/n^2 psi sin(n(s-1/2));psi supported[1/4,3/4],one[3/8,5/8];fixed a>0,all Kn monotone,values[0,1],M<33,N<1217,central third=3712-an',
        repair_pressure='order-one leading w0 fixed;full source moments through Xv in every absolute ledger;cancel only same source in differences;repair coefficients O(n^-2),bumps vanish on cutoff collar',
        viscous_principal='T1theta_vis=2nu X F1X,T1z_vis=2nu X U1X/sqrt(2X);second X derivative has principal coefficients2nu X and nu sqrt(2X);other full residual terms remain',
        nondegeneracy='unbounded third velocity/second stress control requires nonzero fixed inner F1 or U1 at chosen plateau point;not verified for actual profile;no complete-sum cancellation claim',
        schedule='original all-jet Borel prescription already retains order-dependent constants and extra radial curl cutoff derivative;low Cstar alone cannot substitute',
        diagnostics='orders0..8 independent logarithmic/physical identities,orders0..5 actual bump budgets,old low-order calculator comparison,six exact rational auxiliary controls and full-pressure differences',
        scope='conditional finite-order interface plus continuum norm-class obstruction;actual cutoff/scale/pulses/full covariance/residual and unconditional foundations remain open')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('higher_cutoff_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','jet_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py',
        'sharp_current_scale_floor/cost_floor.py','first_seminorm_bound/seminorms.py',
        'cutoff_repair_derivative_budget/derivatives.py','background_jet_foundation/audit.py')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Higher cutoff jet record changed')
    if not all(expected['gates'].values()):raise ValueError('Higher cutoff jet check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual higher cutoff and physical pulse jets remain unverified.')
if __name__=='__main__':main()
