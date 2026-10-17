#!/usr/bin/env python3
"""Join current physical foundations to compatible punctured terminal jets."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
from terminal_bounds import load,atlas_review,seed_review,mixed_review,scope_review
from terminal_checks import integral_review,exact_review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
STATUS='current-source-punctured-terminal-jets-relative-to-retained-core-and-finite-analytic-lemmas'
FLAGS=dict(same_current_single_realization_is_the_field_under_audit=True,
    full_enlarged_closed_band_window_one_quarter_to_four_retained=True,
    unchanged_observation_margin_covers_all_enlarged_fourQ_charts=True,
    actual_seed_axis_potential_and_direct_swirl_factors_joined=True,
    complete_uniform_next_time_derivative_bounds_joined_from_current_foundations=True,
    compatible_potential_direct_swirl_velocity_and_pressure_terminal_jets_supplied_relative_to_retained_lemmas=True,
    all_fixed_mixed_exterior_heat_and_full_pressure_majorants_derived=True,
    both_background_pressure_tails_and_fixed_gauge_memory_retained=True,
    same_labels_cutoff_sequences_supports_width_and_parameters_preserved=True,
    exact_Gaussian_ledger_never_forward_inverted_or_newly_summed=True,
    numerical_Gamma_and_polynomial_controls_are_not_continuum_suprema=True,
    historical_frozen_evidence_preserved_and_replayed=True,
    retained_exact_core_and_finite_existence_inverse_gap_lemmas_still_dependencies=True,
    actual_infinite_coefficients_scales_or_terminal_field_jets_evaluated=False,
    fractional_auxiliary_map_claimed_smooth_at_axis=False,
    analytic_future_heat_or_zero_extension_of_nonzero_pressure_justified=False,
    regular_velocity_or_pressure_at_singular_origin_inferred=False,
    localization_initialization_admissible_global_force_energy_and_comparison_joined=False,
    all_retained_analytic_dependencies_discharged=False,
    unconditional_complete_residual_or_blowup_verified=False,
    formal_certificate=False)
REMAINING=['join_axisymmetric_spatial_time_localization_and_zero_initialization_to_this_same_current_complete_field',
    'retain_every_curl_pressure_viscous_and_nonlinear_term_then_glue_the_force_by_smooth_terminal_jets_not_analytic_or_zero_heat_continuation',
    'join_whole_space_energy_dissipation_and_global_solution_comparison_for_the_same_field_data_and_force',
    'discharge_retained_exact_core_and_finite_analytic_existence_inverse_gap_dependencies_before_unconditional_claims']
DEPENDENCIES=('realization_foundation','physical_pulse_foundation','compact_mean_foundation','annular_exterior_foundation',
    'core_all_order_foundation','coefficient_support_audit','auxiliary_geometry_audit','current_endpoint_jets')

def protocol():
    return dict(parent_commit='b87f53726641b240050a71789c025401e359abaa',
        manuscript_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
        unchanged_parameters='Md64; unit viscosity; same h,N,positive lambda,width C_setup^-20,roots,geometry,labels,cutoff sequences and pressure gauge',
        domain='same current realization restricted to a late observation cylinder; t<=1 away from x=0; no future velocity assertion',
        enlarged_window='1/4<=q/Q<=4; ceil(log2(1/(4*qmax)))<=ell<=floor(log2(4/c)); ell>=ell0',
        observation='q_e=min(qstar/16,2^(-ell0-2)); tau<=q_e/16; |z|<=(q_e/16)^D; 4Q<=16q<=6q_e<=3qstar/8',
        seed='A_theta e_theta=(1/2)*q^(-A)*(M/X)*(-y,x,0); B=q^(-A-1/2)*F*(-y,x,0); p=q^(-2A)*Pi near axis',
        endpoint_input='all fixed complete physical representative jets bounded uniformly on enlarged punctured compacts, including NEXT time derivative',
        exterior='K=C*r^(-1-2h)*H(4*tau/r^2); p=-integral_r^infinity K(rho,tau)^2/rho d rho; same fixed gauge',
        radial_recursion='c_(n+1,k)=-(a+n+2*k)*c_(n,k)+2*c_(n,k-1)',
        pressure_convolution='F_j=(-1)^j*H^(j); G_j=sum binomial(j,k)*F_k*F_(j-k); G_j_prime=-G_(j+1)',
        finite_Taylor='G(1-tau)=sum_(j=0)^k (-tau)^j*F_j/j! plus remainder <=M_(k+1)*tau^(k+1)/(k+1)!',
        diagnostics='36 exact terminal pressure identities; 24 nonterminal Gamma heat and 18 pressure product comparisons; polynomial seed and finite Taylor controls',
        numerical_relative_error_threshold='2e-40',
        scope='compatible punctured one-sided terminal jets relative to retained exact-core and finite analytic lemmas; no admissible global force or blow-up')

@lru_cache(None)
def parent_module():
    p=HERE.parent/'realization_foundation';sys.path.insert(0,str(p))
    try:return load('endpoint_foundation_current_realization_parent',p/'audit.py')
    finally:sys.path.pop(0)

@lru_cache(None)
def parent_review():
    m=parent_module();r=m.generate();m.validate(r)
    f=json.loads((HERE.parent/'realization_foundation/evidence.json').read_text());m.validate(f);m.compare_records(r,f)
    history=dict(r['inherited_provenance'])
    for p,d in f['source_hashes'].items():history['realization_foundation:'+p]=hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==d
    return r,history

def sources():
    paths=[HERE/n for n in ('README.md','terminal_bounds.py','terminal_checks.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'evidence.json' for n in DEPENDENCIES]
    paths += [HERE.parent/'annular_exterior_foundation/exterior_bounds.py',HERE.parent/'realization_foundation/realized_checks.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def reviews():return deepcopy(dict(closed_enlarged_terminal_atlas=atlas_review(),actual_seed_and_uniform_jet_join=seed_review(),
    mixed_exterior_heat_and_complete_pressure_bounds=mixed_review(),independent_nonterminal_Gamma_controls=integral_review(),
    exact_axis_and_finite_Taylor_controls=exact_review(),one_sided_scope_and_quantifiers=scope_review()))

def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    p,history=parent_review();r=reviews();gates={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    gates.update(current_realization_parent_reproduced=all(p['gates'].values()),historical_bytes_preserved=bool(history) and all(history.values()))
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=gates,
        inherited_provenance=dict(history),source_hashes=sources(),remaining_obligations=REMAINING[:])

def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)

def validate(r):
    if r.get('status')!=STATUS or r.get('protocol')!=protocol() or r.get('flags')!=FLAGS:raise ValueError('Current punctured endpoint scope changed')
    if r.get('remaining_obligations')!=REMAINING:raise ValueError('Global or analytic dependencies erased')
    if set(r.get('gates',{}))!=set(reviews())|{'current_realization_parent_reproduced','historical_bytes_preserved'} or not all(r['gates'].values()):raise ValueError('Endpoint foundation gate failed')
    compare_records(reviews(),r.get('reviews',{}),'root.reviews')
    if r.get('source_hashes')!=sources():raise ValueError('Source bytes changed')
    if len(r.get('inherited_provenance',{}))!=979 or not all(r['inherited_provenance'].values()):raise ValueError('Historical provenance changed')
    compare_records(parent_review()[1],r['inherited_provenance'],'inherited_provenance')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; historical evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:
        f=json.loads(a.verify_record.read_text());validate(f);compare_records(r,f)
    a.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} bindings, {len(r['inherited_provenance'])} inherited checks pass")
    print('One-sided punctured terminal jets relative to retained lemmas; global force, energy, comparison and blow-up remain open.')

if __name__=='__main__':main()
