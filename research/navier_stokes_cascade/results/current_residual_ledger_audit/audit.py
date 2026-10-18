#!/usr/bin/env python3
"""Exact current residual ledger and conditional all-fixed-jet tail estimates."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from residual_ledger import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='conditional-current-nonzero-residual-ledger-and-fixed-jet-tail-bound'
DEPENDENCIES=('current_target_cycle_bridge','regular_source_foundation','sourced_inverse_foundation',
 'finite_induction_foundation','compact_mean_foundation','source_envelope_edge_audit','physical_pulse_foundation')
FLAGS=dict(exact_current_regular_plus_additive_tail_identity_keeps_full_cut_pressure_and_potential=True,
    all_old_new_cross_orientations_and_complete_new_velocity_squares_retained=True,
    complete_cut_velocity_divergence_requires_cutoff_curl_in_nonlinear_products=True,
    recomputed_angular_zero_mean_pressure_does_not_reclassify_nonzero_residual=True,
    angular_zero_pressure_and_Pi2_axial_gradient_still_change_full_mean_equation=True,
    actual_mean_increment_nonzero_response_is_both_current_mean_wave_products=True,
    seven_regular_groups_tie_conditional_half_minus_kappa_gain_on_whole_half_line=True,
    Gaussian_tails_have_explicit_all_fixed_physical_jet_constants_on_same_domain=True,
    physically_flat_does_not_supply_improving_normalized_source_bound=True,
    actual_temporal_tail_cutoff_distinguished_from_spatial_slow_localizers=True,
    spatial_localizer_derivatives_retained_in_regular_remainder_not_declared_Gaussian=True,
    conditional_complete_finite_residual_estimate_relative_to_listed_inputs=True,
    same_full_Pi1_Pi2_base_actual_A0_A1_fixed_width_and_corrected_freezing_domain_retained=True,
    historical_provenance_is_byte_binding_not_hidden_reexecution=True,
    manufactured_polynomial_operator_controls_not_actual_inverse_or_continuum_suprema=True,
    retained_actual_initial_core_background_inverse_and_Gaussian_hypotheses_independently_reproved_here=False,
    unconditional_all_family_pulse_or_complete_residual_improvement_proved=False,
    actual_giant_fields_first_band_or_old_force_acceptance_evaluated=False,
    flat_tail_erased_or_forward_inverted=False,
    all_fixed_jet_constants_uniform_in_stage_width_or_lambda=False,
    arbitrary_stage_diagonal_or_single_field_realized_here=False,
    all_spatial_diagonal_or_localization_cutoff_errors_Gaussian_flat=False,
    unconditional_blowup_proved=False,formal_certificate=False)
REMAINING=['identify_and_independently_discharge_retained_actual_canonical_initial_seed_core_background_common_domain_inverse_and_Gaussian_hypotheses_not_protocol_labels',
 'join_this_same_complete_finite_residual_to_the_distinct_accepted_tail_single_field_and_force_with_its_own_controlled_schedules',
 'retain_historical_first_band_and_old_force_acceptance_limits_and_do_not_promote_finite_stage_asymptotics_to_an_uncontrolled_diagonal_or_unconditional_proof']

def protocol():
    return dict(parent_commit='62cce31db067535582f9b22d4ea0727a3c158f04',
      projection='nonzero angular projection of cylindrical components, not auxiliary Haar; full mean pressures angular zero only in this projection',
      exact_wave='Gnew_reg=psi*r+Pperp[mean-new, new-mean, oldwave-new, new-oldwave, new-new]; Tnew=Told+(1-psi)*CURRENT_f+Tpsi',
      cutoff='vpsi=curl(psi A)=psi v+grad psi cross A; pressure psi pi; Tpsi retains time transport viscosity pressure and L_B cutoff curl',
      classification='Gaussian-ledger psi is temporal with frozen Ls and grad_x v_pulse=0; spatial/transverse/slow localizers stay inside uncut potential/pressure and regular r; diagonal/localization cutoffs have separate estimates',
      mean='actual angular-zero mean update changes nonzero residual by N(d,w)+N(w,d); exact nonzero tail ledger unchanged',
      pressure='full n.Kt-nprime.t+n.f normal pressure; pure pressure pi grad psi nonzero; same full Pi1/Pi2 base and current mean moments retained',
      orders='B>=3/2-2kappa; complete linear/cross/fast rows B+1/2-kappa, squares 2B-kappa, slow mean B+1/2',
      Gaussian='q^-D ell^b exp(-c ell^2)<=n! (2/c)^n exp((D+N)^2/(2c)) q^N, ell>=1,n=ceil(b/2)',
      obstruction='physical flatness does not imply W_(B+delta): exp(h delta ell)<=C ell^a fails; old tails never inverse sources',
      physical='conditional fixed j,m residual q^(h Bnext-Km),Km=2A+3/2+2m; full base/exterior bounds still explicit inputs',
      domain='same corrected initial tail with ell^6>=8*Cfreeze/gamma0 and frozen-margin Lreq; actual C_setup^-20 width and radius retained',
      validation='exact independent cylindrical Fourier residual checks; relevant historical numerical evidence replayed separately, recursive hashes alone do not reexecute it',
      scope='conditional complete finite ledger and fixed-jet bound only; retained initial/core/background/inverse hypotheses and single-field forcing remain unpromoted')

@lru_cache(None)
def parent_module():
    p=HERE.parent/'current_target_cycle_bridge';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('residual_current_target_parent',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('current_residual_byte_provenance',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():
    # Relevant parent test/evidence replays are separate CI/local commands.
    # This function verifies every bound historical source byte, not its numerics.
    return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','residual_ledger.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/'finite_induction_foundation/finite_checks.py',HERE.parent/'first_sourced_foundation/sourced_checks.py',
              HERE.parent/'complete_covariance_initial_means/means.py',HERE.parent/'regular_source_foundation/source_checks.py',
              HERE.parent/'harmonic_closure_foundation/audit.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
@lru_cache(None)
def reviews():return review()
def generate():
    if json.loads((HERE/'protocol.json').read_text())!=protocol():raise ValueError('Protocol changed')
    r=deepcopy(reviews());h=history();g=dict(r['checks'])
    g['historical_recursive_source_bytes']=bool(h['bindings']) and all(h['bindings'].values())
    return dict(status=STATUS,protocol=protocol(),flags=dict(FLAGS),reviews=r,gates=g,
        inherited_provenance=deepcopy(h),source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Current-residual record changed')
    if not all(expected['gates'].values()):raise ValueError('Current-residual gate failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; historical frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f"{STATUS}: {len(r['gates'])} gates; {len(r['reviews']['exact_complete_residual_ledger']['identities'])} exact residual identities")
    print('Conditional finite residual only; actual initial/core/background/inverse premises and unconditional proof remain.')
if __name__=='__main__':main()
