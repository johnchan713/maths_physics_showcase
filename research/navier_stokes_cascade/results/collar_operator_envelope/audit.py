#!/usr/bin/env python3
"""Reproduce explicit actual-collar operator/source input, not a full proof."""
import argparse
from copy import deepcopy
from functools import lru_cache
import hashlib,importlib.util,json,sys
from pathlib import Path
from envelope import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
STATUS='explicit-actual-activation-collar-operator-and-full-source-envelope'
DEPENDENCIES=('leading_cone_margin','core_all_order_foundation','activated_collar_foundation',
    'first_correction_construction','first_seminorm_bound')
FLAGS=dict(actual_activation_used_instead_of_reference_or_protected_core_bound=True,
    original_squared_flat_cutoff_width_and_all_supports_unchanged=True,
    explicit_second_radial_and_common_complex_parameter_jet_envelopes_derived=True,
    regular_axis_radial_average_and_reciprocal_L_derivative_retained=True,
    full_radial_pressure_source_viscosity_and_shifted_axial_sources_retained=True,
    complete_local_pressure_incoming_axis_datum_preserved=True,
    explicit_collar_operator_and_Picard_strip_inputs_derived=True,
    exponent_separated_diagnostics_do_not_materialize_actual_giant_parameters=True,
    numerical_controls_are_manufactured_not_infinite_core_samples=True,
    actual_giant_parameter_values_numerically_evaluated=False,
    complete_future_pressure_integrals_and_global_first_coefficient_inputs_bounded=False,
    actual_existing_tilt_gate_verified=False,
    combined_scale_gate_or_historical_first_band_certified=False,
    independent_unconditional_core_and_global_construction_proved=False,
    unconditional_actual_pulses_complete_residual_or_blowup_proved=False,
    formal_certificate=False)
REMAINING=['quantify_global_leading_jets_full_future_W0_W1_shifted_axial_sources_and_actual_cutoff_repair_bump_envelopes',
    'verify_existing_unchanged_tilt_gate_or_sharpen_same_u_raw_cone_bound',
    'substitute_complete_Cstar_and_same_u_gap_radius_then_check_all_scale_freezing_geometry_and_inverse_gates',
    'retain_complete_wave_covariance_mean_axial_velocity_and_both_pressure_memories_in_full_residual_closure']

def protocol():
    return dict(parent_commit='1e910ccc6aa4da81eefb1d05e76e010943d0e91b',
        domain='real 0<=X<=Xa exp(tc/4), common complex eta strip; actual flat activation and unchanged core through Xa',
        inputs='S>=2, 0<rho<=1, Lambda>0, t=C^-20<=.01, ell=tc/4<t/20, construction viscosity=1',
        reference='x=4.1/Lambda,b=1+x,B=2(256S+9+S/Lambda)(1+Lambda)^2; delta=rho/(800S)',
        activation='j1=25S,j2=1600S/t+1000S^2,E=exp(25ell S); K given by README equation (4)',
        angular='H=K(1+delta^-1+2delta^-2) bounds radial orders 0..2 and eta orders 0..2 on 2delta',
        full_source='W=512(1+x)H^2,Z=128(1+x)^2H; f=(0,0,0,-xi*w0,-2CF,2eta X*w0/L-2CU)',
        operator='Ccol=4096(1+x+b)^3 H^2 bounds infinity row norms A0,A1 and vector source on 2delta',
        pressure='Pi=Pi_r(Xa)+int_0^y X F_a^2; Pcore+ell*x*K^2 value bound preserves incoming Pi_axis',
        Picard='sigma=2delta,Delta=d=delta/4,z=2Ccol*b,beta=e*z^2/Delta,M=(z+beta)exp(beta)',
        arithmetic='E<=2^ceil(50ell S); rational fixtures store the exponent separately, not actual core values',
        scope='conditional local collar input only; global pressure source integrals, fixed tilt, full scale acceptance and full residual remain open')

@lru_cache(None)
def provenance_tools():
    p=HERE.parent/'harmonic_closure_foundation';old=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('collar_operator_historical_bytes',p/'audit.py')
        m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
    finally:sys.path[:]=old
@lru_cache(None)
def history():return provenance_tools().verify_bindings(PROJECT,[HERE.parent/n/'evidence.json' for n in DEPENDENCIES])
def sources():
    files=[HERE/n for n in ('README.md','audit.py','envelope.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    files += [HERE.parent/n/f for n in DEPENDENCIES for f in ('README.md','evidence.json')]
    files += [HERE.parent/'harmonic_closure_foundation/audit.py',HERE.parent/'first_correction_construction/system.py',
        HERE.parent/'first_correction_construction/majorant.py',HERE.parent/'first_seminorm_bound/seminorms.py',
        HERE.parent/'physical_residual_budget/residual.py',HERE.parent/'physical_residual_budget/review.py']
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
    if json.dumps(expected,sort_keys=True)!=json.dumps(r,sort_keys=True):raise ValueError('Collar envelope record changed')
    if not all(expected['gates'].values()):raise ValueError('Collar envelope check failed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(f'{STATUS}: {len(r["gates"])} checks; global source/tilt/scale/residual acceptance not inferred.')
if __name__=='__main__':main()
