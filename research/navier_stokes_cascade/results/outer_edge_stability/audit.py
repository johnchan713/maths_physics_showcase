#!/usr/bin/env python3
"""Reproduce the direct outer-edge and first viscous-stress audit."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
from diagnostics import (algebra_review, edge_review, failed_shortcuts_review,
                         heat_review, physical_source_review, relative,
                         stress_review, uniform_bounds_review)

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = 'outer-heat-edge-and-first-viscous-stress-stability-audited'
FLAGS = dict(outer_heat_profile_growth_and_oriented_cone_argument_supplied=True,
             first_outer_viscous_stress_positive_angular_and_zero_axial=True,
             outer_T1_over_T0_ratio_bounded_and_smooth_at_flat_edge=True,
             nonnegative_cutoff_preserves_outer_background_stress_cone=True,
             actual_giant_schedule_or_global_flow_numerically_evaluated=False,
             global_forward_backward_moment_identification_independently_reverified=False,
             inner_core_collar_and_compact_cone_foundation_independently_verified=False,
             unconditional_uniform_pulse_estimate_verified=False,
             unconditional_complete_residual_improvement_verified=False,
             blowup_verified=False, formal_proof_assistant_certificate=False)
REVIEWS = dict(exact_axial_coefficient_algebra=algebra_review,
               whole_interval_heat_and_schedule_bounds=uniform_bounds_review,
               exact_heat_integral=heat_review, signed_and_positive_viscous_stress=stress_review,
               physical_axial_source_and_radial_stress_sign=physical_source_review,
               flat_edge_and_correction_ratio=edge_review, failed_shortcuts=failed_shortcuts_review)
REMAINING = ('verify_actual_inner_analytic_core_and_inner_stress_direction_factor_at_same_N_h',
             'independently_review_compact_modulation_and_repair_cone_certificate_with_pressure_memory',
             'verify_global_exact_moments_and_coefficient_support_identifying_the_outer_backward_stresses',
             'audit_reserved_mean_supports_localization_and_complete_residual_assembly')


def protocol():
    return dict(status=STATUS, parent_commit='f31516d84c8d46137e41c0882f5006e5d938d30a',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['4.2', '4.11', '5.3', '5.9', 'A.47-A.56'],
                same_fixed_N_and_h=True, proof_domain='1/2<=y=log(X/Xtail)<3, -1<=eta<=1, Xtail>=100',
                proof_h_range='0<h<1/100', actual_schedule='Md=64, T=exp(64)+10, h=exp(-8*T)',
                diagnostic_h='.005', diagnostic_Xtail=1000, diagnostic_amplitude='1e-10',
                diagnostic_precision=80, interval_precision=70,
                kernel='exact Gamma integral; numerical Tricomi U with independent positive integration',
                cutoff='chi(s)=1-sigma(2*s-1), 0<=chi<=1',
                first_stress='T1theta=alpha*r*K*fprime/2 + integral r_v^3*K_v*C*fprime/(2*r^2); T1z=0',
                edge_ratio_limit='2*eta^2*Xb/L^2', uniform_ratio_bound='2*Xb/.98^2+7',
                flat_deltas=['.3', '.03', '.001', '1e-7'], flat_eta=['-1', '0', '.5', '1'],
                physical_q=['1', '.4', '.2', '.7'], radial_difference_step='1e-8',
                global_identification_inputs=['leading exact total moments', 'positive-order exact moments and common supports'],
                diagnostics_are_not_continuum_proofs=True)


def inherited_hashes():
    folders = ('paper_profile_audit', 'outer_pressure_pilot', 'axial_stress_audit', 'intermediate_decay_audit',
               'pulse_moment_audit', 'pulse_stress_audit', 'post_pulse_stress_audit', 'heat_exterior_audit',
               'axis_matching_audit', 'axis_core_attachment', 'stress_realization_audit', 'joined_stress_construction',
               'angular_c2_transfer', 'compact_jet_envelope', 'loop_modulation_bounds', 'repair_state_bounds',
               'physical_residual_budget', 'first_correction_construction', 'higher_jet_budget',
               'annular_first_correction', 'background_summation', 'physical_wave_interface',
               'actual_band_pulse_bounds', 'compact_mean_corrections', 'residual_cycle_closure',
               'full_potential_realization', 'localized_force_extension', 'global_solution_comparison',
               'auxiliary_geometry_audit', 'geometric_pulse_replay', 'background_jet_foundation')
    checks = {}
    for folder in folders:
        record = json.loads((HERE.parent/folder/'evidence.json').read_text())
        if 'source_sha256' in record:
            hashes, root = record['source_sha256'], HERE.parent/folder
        elif 'source_hashes' in record:
            hashes, root = record['source_hashes'], PROJECT
        else:
            hashes, root = record.get('provenance', {}), None
        if not hashes:
            raise ValueError('Missing inherited ledger: '+folder)
        for path, expected in hashes.items():
            routed = root if root is not None else (PROJECT if path.startswith('results/') or path == 'GOAL.md' else HERE.parent)
            checks[folder+':'+path] = hashlib.sha256((routed/path).read_bytes()).hexdigest() == expected
    return checks


def source_hashes():
    names = ('README.md', '.gitignore', 'requirements.txt', 'protocol.json', 'tail.py',
             'diagnostics.py', 'audit.py', 'test_audit.py')
    paths = [HERE/name for name in names]+[PROJECT/'GOAL.md',
             HERE.parent/'background_jet_foundation/evidence.json', HERE.parent/'heat_exterior_audit/evidence.json',
             HERE.parent/'annular_first_correction/evidence.json']
    return {str(p.relative_to(PROJECT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def refinement_review(reviews):
    stress, edge = stress_review(True), edge_review(True)
    worst = mp.mpf(0)
    for a, b in zip(reviews['signed_and_positive_viscous_stress']['rows'], stress['rows']):
        for key in ('T0theta', 'T0z', 'T1positive', 'T1signed', 'correction_ratio'):
            worst = max(worst, relative(a[key], b[key]))
    for a, b in zip(reviews['flat_edge_and_correction_ratio']['rows'], edge['rows']):
        for key in ('b0', 'b1', 'bz', 'correction_ratio'):
            worst = max(worst, relative(a[key], b[key]))
    return dict(maximum_relative_refinement_error=float(worst),
                check=bool(worst < mp.mpf('1e-8') and all(stress['checks'].values()) and all(edge['checks'].values())))


def generate():
    if json.loads((HERE/'protocol.json').read_text()) != protocol():
        raise ValueError('Frozen outer-edge protocol changed')
    with mp.workdps(80):
        reviews = {name: fn() for name, fn in REVIEWS.items()}
        refinement = refinement_review(reviews)
    inherited = inherited_hashes()
    gates = {name: all(review['checks'].values()) for name, review in reviews.items()}
    gates['tighter_quadrature_agrees'] = refinement['check']
    gates['inherited_bytes_preserved'] = bool(inherited) and all(inherited.values())
    gates['global_and_inner_foundation_not_promoted'] = not any(FLAGS[name] for name in
        ('actual_giant_schedule_or_global_flow_numerically_evaluated',
         'global_forward_backward_moment_identification_independently_reverified',
         'inner_core_collar_and_compact_cone_foundation_independently_verified',
         'unconditional_uniform_pulse_estimate_verified',
         'unconditional_complete_residual_improvement_verified', 'blowup_verified'))
    return dict(status=STATUS, protocol=protocol(), flags=FLAGS, reviews=reviews, refinement=refinement,
                gates=gates, remaining_obligations=list(REMAINING),
                inherited_provenance=inherited, source_hashes=source_hashes())


def validate(record):
    if record.get('status') != STATUS or record.get('protocol') != protocol() or record.get('flags') != FLAGS:
        raise ValueError('Outer scope, sign, scale or assumptions changed')
    if record.get('remaining_obligations') != list(REMAINING):
        raise ValueError('Inner, compact or global identification obligation omitted')
    expected = set(REVIEWS)|{'tighter_quadrature_agrees', 'inherited_bytes_preserved', 'global_and_inner_foundation_not_promoted'}
    if set(record.get('gates', {})) != expected or not all(record['gates'].values()):
        raise ValueError('An outer-edge audit gate failed')
    if not record.get('inherited_provenance') or not all(record['inherited_provenance'].values()):
        raise ValueError('Inherited bytes changed')
    if set(record.get('reviews', {})) != set(REVIEWS):
        raise ValueError('Missing direct or failure-control review')
    if not all(r.get('checks') and all(r['checks'].values()) for r in record['reviews'].values()):
        raise ValueError('A required check failed')
    if not record.get('refinement', {}).get('check'):
        raise ValueError('Quadrature refinement failed')


def compare_records(actual, frozen, path='root'):
    if isinstance(actual, dict):
        if not isinstance(frozen, dict) or set(actual) != set(frozen):
            raise ValueError('Changed keys at '+path)
        for key in actual:
            compare_records(actual[key], frozen[key], path+'.'+key)
    elif isinstance(actual, list):
        if not isinstance(frozen, list) or len(actual) != len(frozen):
            raise ValueError('Changed sequence at '+path)
        for i, value in enumerate(actual):
            compare_records(value, frozen[i], path+'.'+str(i))
    elif isinstance(actual, float):
        if type(frozen) is not float or not math.isfinite(actual) or not math.isfinite(frozen):
            raise ValueError('Nonfinite or mistyped diagnostic at '+path)
        strict = path.endswith(('.T0theta', '.T0z', '.T1positive', '.T1signed', '.b0', '.b1', '.bz', '.directional_loss'))
        absolute = 1e-300 if strict else (1e-40 if 'source_normalized_error' in path or 'heat_ODE_error' in path else 1e-12)
        if not math.isclose(actual, frozen, rel_tol=3e-10, abs_tol=absolute):
            raise ValueError('Changed numerical result at '+path)
    elif type(actual) is not type(frozen) or actual != frozen:
        raise ValueError('Changed exact result at '+path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--verify-record', type=Path)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit('Output already exists; choose a fresh path')
    record = generate()
    validate(record)
    if args.verify_record:
        frozen = json.loads(args.verify_record.read_text())
        validate(frozen)
        compare_records(record, frozen)
    args.output.write_text(json.dumps(record, sort_keys=True, indent=2)+'\n')
    print(STATUS)
    print(f"{len(record['gates'])} gates pass; {len(record['inherited_provenance'])} inherited hashes agree.")
    print('Inner/compact/global identification, unconditional pulses and the complete residual remain unverified.')


if __name__ == '__main__':
    main()
