#!/usr/bin/env python3
"""Reproduce residual identities without promoting the open correction proof."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from residual import HERE, PROJECT
from bounds import ledger
from review import (BASE, CORRECTION, physical_review, omission_review,
                    derivative_gap_review, cutoff_review, exponent_review)

STATUS = 'physical-residual-derived-with-open-C3-and-correction-bounds'
FLAGS = dict(full_unlocalized_axisymmetric_residual_identity_derived=True,
             first_positive_order_equations_derived=True,
             physical_derivative_exponent_ledger_derived=True,
             inherited_C2_does_not_bound_full_residual=True,
             relies_on_inherited_analytic_construction=True,
             actual_post_modulation_C3_bound_verified=False,
             first_correction_for_actual_profile_solved=False,
             actual_high_frequency_profile_numerically_resolved=False,
             full_physical_residual_budget_verified=False,
             full_physical_stress_wave_construction_verified=False,
             full_PDE_corrections_verified=False, smooth_force_verified=False,
             blowup_verified=False, independent_foundation_review_completed=False,
             formal_proof_assistant_certificate=False)
RESOLVED = ['unlocalized_axisymmetric_residual_identity',
            'first_positive_order_coefficient_equations',
            'physical_derivative_exponent_ledger']
REMAINING = ['actual_post_modulation_radial_and_C3_angular_bounds',
             'positive_order_solvability_and_five_moment_extension',
             'full_wave_covariance_and_all_nonlinear_interactions',
             'all_order_PDE_correction_summation',
             'localization_and_smooth_forcing_through_T',
             'complete_field_blowup_lower_bound_and_energy',
             'independent_review_of_inherited_continuum_estimates']


def expected_protocol():
    return dict(status=STATUS, parent_commit='753501b7c205e6c9073e1c939bbc340dfbf65dc1',
                paper_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f',
                source_equations=['4.2', '4.3-4.7', '4.11-4.12', '5.1-5.6', '5.26', '5.41'],
                fixed_N='1+floor(H^32)',
                scope='Unlocalized axisymmetric base and finite coefficient identities; manufactured physical crosschecks; actual derivative constants, coefficient solvability, waves and summation remain open.',
                physical_sample_count=324, maximum_identity_error=2e-11,
                negative_control_minimum_gap=1e-4, maximum_finite_difference_error=1e-6,
                minimum_finite_difference_refinement=4,
                C2_counterexample_frequencies=[8, 32, 128, 512],
                phase_counterexample_frequencies=[8, 64, 512],
                diagnostic_cutoff='exp(-sigma); not compactly supported',
                required_actual_angular_order_for_full_base_residual=3,
                actual_profile_numerically_resolved=False)


def provenance():
    prior = json.loads((HERE.parent/'repair_state_bounds'/'evidence.json').read_text())
    return {path: hashlib.sha256((PROJECT/path).read_bytes()).hexdigest() == digest
            for path, digest in prior['source_hashes'].items()}


def hashes():
    paths = [HERE/name for name in ('README.md', 'protocol.json', 'requirements.txt',
                                   'residual.py', 'review.py', 'bounds.py', 'audit.py', 'test_audit.py')]
    paths += [HERE.parent/'repair_state_bounds'/'evidence.json',
              HERE.parent/'paper_profile_audit'/'profiles.py',
              HERE.parent/'paper_profile_audit'/'protocol.json']
    return {str(path.relative_to(PROJECT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def encode(value):
    if isinstance(value, dict):
        return {key: encode(v) for key, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, (float, np.floating)):
        if not math.isfinite(value):
            raise ValueError('Nonfinite diagnostic')
        return float(value)
    return value


def generate():
    protocol = json.loads((HERE/'protocol.json').read_text())
    if protocol != expected_protocol():
        raise ValueError('Frozen scope, source or derivative budget changed')
    reviews = dict(physical=physical_review(), omissions=omission_review(),
                   derivative_gap=derivative_gap_review(), cutoff=cutoff_review(),
                   exponents=exponent_review())
    inherited = provenance()
    constants = ledger()
    gates = dict(physical_coordinate_and_coefficient_identities=all(reviews['physical']['checks'].values()),
                 scalar_stencil_and_omission_controls=all(reviews['omissions']['checks'].values()),
                 C2_and_large_N_insufficiency_detected=all(reviews['derivative_gap']['checks'].values()),
                 divergence_preserving_cutoff_product_rule=all(reviews['cutoff']['checks'].values()),
                 derivative_loss_and_finite_order_limits=all(reviews['exponents']['checks'].values()),
                 conditional_operator_constant_arithmetic=all(constants['checks'].values()) and not constants['actual_K_evaluated'],
                 inherited_bytes_preserved=bool(inherited) and all(inherited.values()),
                 manufactured_fixtures_and_noncompact_cutoff_explicit=all(reviews[key]['manufactured'] for key in ('physical', 'omissions', 'derivative_gap', 'cutoff')) and not reviews['cutoff']['diagnostic_cutoff_is_compact'],
                 actual_PDE_and_blowup_unverified=not any(FLAGS[key] for key in ('actual_post_modulation_C3_bound_verified', 'first_correction_for_actual_profile_solved', 'full_physical_residual_budget_verified', 'full_PDE_corrections_verified', 'smooth_force_verified', 'blowup_verified')))
    return encode(dict(status=STATUS, flags=FLAGS, protocol=protocol,
                       fixtures=dict(base=BASE, correction=CORRECTION),
                       conditional_constant_ledger=constants,
                       source_hashes=hashes(), inherited_provenance=inherited,
                       resolved_obligations=RESOLVED, remaining_obligations=REMAINING,
                       reviews=reviews, gates=gates))


def validate(record):
    if record.get('status') != STATUS or record.get('flags') != FLAGS:
        raise ValueError('Residual identities cannot promote actual corrections or blow-up')
    if record.get('resolved_obligations') != RESOLVED or record.get('remaining_obligations') != REMAINING:
        raise ValueError('Unclosed mathematical obligations were changed')
    if not record.get('gates') or not all(record['gates'].values()):
        raise ValueError('A residual identity, negative control or provenance check failed')


def compare_records(actual, frozen, location='root'):
    """Only numerical roundoff may differ between supported NumPy/SciPy builds."""
    if isinstance(actual, bool) or isinstance(actual, (str, int)) or actual is None:
        if type(actual) is not type(frozen) or actual != frozen:
            raise ValueError('Frozen evidence differs at '+location)
    elif isinstance(actual, float):
        if not isinstance(frozen, (int, float)) or isinstance(frozen, bool) or not math.isclose(actual, frozen, rel_tol=1e-8, abs_tol=2e-10):
            raise ValueError('Frozen numerical evidence differs at '+location)
    elif isinstance(actual, dict):
        if not isinstance(frozen, dict) or set(actual) != set(frozen):
            raise ValueError('Frozen fields differ at '+location)
        for key in actual:
            compare_records(actual[key], frozen[key], location+'.'+key)
    elif isinstance(actual, list):
        if not isinstance(frozen, list) or len(actual) != len(frozen):
            raise ValueError('Frozen samples differ at '+location)
        for i, (a, b) in enumerate(zip(actual, frozen)):
            compare_records(a, b, location+'['+str(i)+']')
    else:
        raise TypeError('Unsupported evidence value')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--verify-record', type=Path)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit('Output already exists; choose a fresh path')
    result = generate()
    validate(result)
    if args.verify_record:
        frozen = json.loads(args.verify_record.read_text())
        validate(frozen)
        compare_records(result, frozen)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(STATUS)
    print('Residual and first-order identities pass; the actual C3 bound and correction construction remain open.')


if __name__ == '__main__':
    main()
