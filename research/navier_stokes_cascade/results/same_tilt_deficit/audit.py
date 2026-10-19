#!/usr/bin/env python3
"""Same original tilt, explicit deficit budget and independent freezing; no actual acceptance."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from deficit_bounds import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='same-tilt-raw-deficit-radius-and-joint-freezing-scale-gate'
DEPENDENCIES=('cutoff_shape_feasibility','cutoff_cone_scale_gate','background_gap_certificate',
    'leading_cone_margin','source_envelope_edge_audit','primary_covariance_foundation')
FLAGS=dict(original_strict_tilt_choice_implies_qualitative_raw_deficit_relative_to_retained_leading_lemmas=True,
    new_normalized_cone_radius_bypasses_only_stronger_coarse_sufficient_tilt_test=True,
    explicit_deficit_requirement_contains_both_full_background_error_and_independent_freezing=True,
    complete_corrected_width_pressure_and_viscosity_rows_and_sharp_shape_program_retained=True,
    full_four_Q_chart_two_band_cost_and_dyadic_rounding_are_retained=True,
    same_original_tilt_carriers_roots_cutoff_coefficients_frequency_supports_gauges_and_mean_memories_unchanged=True,
    fixed_order_value_gates_not_reselected_for_derivative_order=True,
    historical_frozen_evidence_and_original_comparators_preserved=True,
    original_admissibility_implies_the_stronger_coarse_numeric_tilt_test=False,
    actual_raw_supremum_deficit_cutoff_or_Cfreeze_numerically_verified=False,
    retained_scale_or_all_other_geometry_length_inverse_and_pulse_gates_verified=False,
    auxiliary_sharp_shape_program_identifies_unchanged_cutoff=False,
    sufficient_gate_failure_establishes_actual_cone_or_PDE_failure=False,
    actual_huge_profile_or_physical_pulses_numerically_resolved=False,
    retained_core_and_regional_premises_independently_proved_unconditional=False,
    unconditional_pulse_complete_residual_improvement_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['supply_actual_unchanged_cutoff_derivative_enclosures_and_quantitative_lower_bound_d0_for_original_raw_deficit',
    'bound_relevant_Cfreeze_and_verify_retained_scale_joint_deficit_gate_plus_all_geometry_length_inverse_and_carrier_conditions',
    'prove_uniform_actual_physical_pulse_value_derivative_support_and_inverse_gap_on_every_scale',
    'retain_complete_wave_covariance_axial_velocity_pressure_and_means_in_full_viscous_cutoff_interaction_residual',
    'independently_verify_retained_exact_core_and_regional_foundations_before_unconditional_claims']
def protocol():
    return dict(parent_commit='abfaf6c133a1b24a88c3eb55f54ac0546178a770',
        original_tilt='same u>0 fixed before manuscript (7.2),printed page74; kappa=u^2/(1+u^2)>Ract=sup c0^2*(K0.n0)^2/H0^2,including flat-edge limits; actual deficit positive but unquantified',
        normalized_polynomial='Psi=kappa*H^2-c^2*J^2=Phi/(1+u^2);Lpsi=3kappa LH+4Lc+3cbar2 LH; geometry constants exactly those of retained gap certificate',
        deficit_radius='for valid 0<d0<=kappa-Ract, rd=min(rb,omega^2*d0/(2Lpsi));Psi>=omega^2*d0/2; oriented_unit_stress>=omega/3',
        cone_margin='gamma_B=omega^2*d0/(9kappa), gamma_frozen=gamma_B/2; no tilt or physical pulse is renormalized',
        freezing='same-coordinate tau_freeze<=Cfreeze ell^-6;ell^6>=72kappa*Cfreeze/(omega^2*d0), independently of background error',
        precise_joint_gate='E=Cstar(M,N)*2^(-2h*(ell-2));E<=rb;Ract+max(2Lpsi E/omega^2,72kappa Cfreeze/(omega^2 ell^6))<=kappa',
        rational_controls='replace exponent by its integer floor for a conservative error upper bound; exact sixth-root ceiling; later witness is conditional, not old first-band acceptance',
        thin_margin_family='fixed F=1,g=(-82/25,0),c^2=16/25,D=2624/625,u=3/4; unit directions n(z) and n(-z),0<z<=1/32; deficit<=3z and relative error=4z/(1+z^2)',
        diagnostics='96 exact rational perturbations,3 legacy constant/radius comparisons,4 fixed-tilt thin-margin controls,symbolic identities and complete pressure-shape-freezing composition; not actual profile samples',
        scope='conditional same-tilt value radius and joint quantitative obstruction; actual d0,cutoff,Cfreeze,retained scale and complete physical pulses/residual/blowup not verified')
@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        s=importlib.util.spec_from_file_location('same_tilt_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','deficit_bounds.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/p for p in ('harmonic_closure_foundation/audit.py','background_gap_certificate/gap_bounds.py',
        'leading_cone_margin/cone_bounds.py','cutoff_cone_scale_gate/gate_bounds.py','cutoff_shape_feasibility/shape_bounds.py',
        'actual_band_pulse_bounds/README.md','weighted_background_edges/README.md')]
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Same-tilt deficit record changed')
    if not all(expected['gates'].values()):raise ValueError('Same-tilt deficit check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; actual quantitative deficit, cutoff, Cfreeze and retained-scale acceptance remain unverified.')
if __name__=='__main__':main()
