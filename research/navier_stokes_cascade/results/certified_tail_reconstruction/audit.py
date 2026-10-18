#!/usr/bin/env python3
"""Certified-tail reconstruction with fresh, explicitly nonidentified diagonals and force."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,json
from pathlib import Path
import sys
from tail_bounds import load,domain_review,scale_review,core_and_force_review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='conditional-certified-tail-reconstruction-with-distinct-field-and-force'
DEPENDENCIES=('canonical_value_domain','canonical_finite_cycle_join','canonical_local_realization',
    'canonical_endpoint_join','canonical_force_join','canonical_comparison_join')
FLAGS=dict(actual_finite_real_value_budgets_define_accepted_tail_analytically=True,
    numerical_actual_first_band_not_needed_for_distinct_existential_tail_theorem=True,
    common_tail_domain_is_independent_of_fixed_jet_order_and_finite_stage=True,
    fresh_realization_diagonal_has_whole_interval_support_and_tail_control=True,
    fresh_Borel_diagonal_controls_changed_localization_terminal_norms=True,
    optional_history_dominating_scales_do_not_identify_fields_or_forces=True,
    same_protected_core_and_exact_exterior_retained=True,
    complete_pressure_axial_means_cutoff_curls_viscosity_and_interactions_retained=True,
    fresh_compact_force_energy_and_same_force_comparison_joined_conditionally=True,
    precise_old_Borel_scale_reuse_inequality_and_countercontrol_given=True,
    historical_source_bytes_fields_forces_and_schedules_preserved=True,
    analytical_derivations_distinguished_from_finite_parameter_controls=True,
    original_first_band_or_old_localized_force_acceptance_verified=False,
    new_annular_field_identified_with_historical_realization=False,
    new_localized_force_identified_with_historical_force=False,
    actual_infinite_fields_norms_or_diagonal_scales_numerically_evaluated=False,
    finite_cycle_generated_source_and_edge_hypotheses_all_independently_discharged=False,
    future_velocity_or_pressure_constructed=False,
    proof_assistant_absence_treated_as_mathematical_obstruction=False,
    unconditional_pulse_complete_residual_or_blowup_proved=False,formal_certificate=False)
REMAINING=['audit_actual_canonical_finite_cycle_generated_sources_complete_pressure_normal_inverse_class_weighted_edges_and_support_through_signed_and_mean_updates',
    'identify_and_discharge_every_still_unverified_defining_map_hypothesis_before_unconditional_claims_not_merely_demand_a_formal_certificate',
    'old_first_band_and_old_force_remain_uncertified_but_are_not_inputs_to_distinct_tail_witness']

def protocol():
    return dict(parent_commit='e34f6a5810cd2e59564a3540470d20f67f0cb5b7',
        domain='E=ceil(max(16,ell_old,512/h^2,sqrt(Lreq/(2r0)))); qv=min(qstar/64,2^(-E-2)); actual finite real budgets, not manufactured values',
        realization='fresh a_j_sharp; fixed-gauge complete potential, separate direct swirl and pressure differences; whole-interval constraints on 0<q<=1/a_j_sharp',
        scale='s_j>max(previous+1,ln(1/qv)+1,1,2(c+j)/gamma+P^2/gamma^2); optional old log scale only a lower bound',
        observation='qe=qv/16,tau0=qe/16,z0=(qe/16)^D; q<=3qv/128,4Q<=3qv/8<qv on support closure',
        force='new c_sharp and f_sharp; all twelve residual terms, full pressure gradient and unit viscosity; old force values/norms not reused',
        Borel='new n_j>max(previous+1,2,(c_jam_sharp+j+2)/(j-m)); optional old exponent only lower bound; one sequence for all fixed jets',
        obstruction='old exponents require c_jam_sharp+(m-j)*n_j_old<=-j-2; old norm inequalities do not imply new ones',
        comparison='new datum zero,new compact force,new energy and R^-2 comparison constants; fix T then radius limit then protected core',
        fixed='old qstar,h,N,lambda,pulse radius,width C_setup^-20,phase labels/carriers,operator gauges and historical schedules unchanged as historical objects',
        scope='distinct conditional existential tail witness; old first band not verified; actual finite-cycle/source/edge hypotheses remain; no unconditional theorem')

@lru_cache(None)
def parent_module():
    p=HERE.parent/'canonical_value_domain';old=sys.path[:];sys.path.insert(0,str(p))
    try:return load('certified_tail_value_parent',p/'audit.py')
    finally:sys.path[:]=old

@lru_cache(None)
def records():
    r={n:json.loads((HERE.parent/n/'evidence.json').read_text()) for n in DEPENDENCIES}
    parent_module().validate(r['canonical_value_domain']);return r

@lru_cache(None)
def history():return parent_module().parent_module().parent_module().provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])

def sources():
    paths=[HERE/n for n in ('README.md','audit.py','tail_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    paths += [HERE.parent/'canonical_local_realization/realization.py',HERE.parent/'canonical_force_join/force_bounds.py',
        HERE.parent/'canonical_value_domain/acceptance.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def routing_review():
    r=records()
    return dict(checks=dict(
        parent_does_not_certify_old_band=not r['canonical_value_domain']['flags']['original_first_band_or_old_localized_force_value_acceptance_verified'],
        finite_cycle_is_uniform_in_value_domain=r['canonical_finite_cycle_join']['flags']['all_fixed_jets_same_value_domain_and_original_supports_width_and_carriers'],
        finite_cycle_keeps_complete_pressure_and_squares=r['canonical_finite_cycle_join']['flags']['complete_pressure_viscosity_all_six_covariances_and_nonlinear_squares_retained'],
        endpoint_requires_next_jet_not_residual_alone=not r['canonical_endpoint_join']['flags']['flat_residual_or_bounded_velocity_alone_used_to_infer_terminal_trace'],
        force_extension_has_no_interior_Gevrey_assumption=not r['canonical_force_join']['flags']['complete_force_Gevrey_bound_or_analytic_Taylor_convergence_assumed'],
        all_records_retain_nontrivial_remaining_obligations=all(r[n]['remaining_obligations'] for n in DEPENDENCIES)))

@lru_cache(None)
def reviews():return dict(actual_value_tail_and_support_margins=domain_review(),fresh_two_diagonal_controls=scale_review(),
    core_plateau_and_force_nonidentity=core_and_force_review(),operator_routing_not_field_identification=routing_review())

def compare_records(a,b,path='root'):return parent_module().compare_records(a,b,path)

def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=deepcopy(reviews());h=history();g={k:bool(v['checks']) and all(v['checks'].values()) for k,v in r.items()}
    g['historical_recursive_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate();compare_records(expected,r)
    if not all(expected['gates'].values()):raise ValueError('Certified-tail reconstruction gate failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; historical evidence is not overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates, {len(r['source_hashes'])} direct bindings, {len(r['inherited_provenance']['bindings'])} historical bindings pass")
    print('Distinct conditional tail witness only; old force and first band not recertified, actual source/edge audit and unconditional proof remain.')
if __name__=='__main__':main()
