#!/usr/bin/env python3
"""Reproduce repair constants and conditional frequency inequalities, not a PDE proof."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import mpmath as mp
from bounds import HERE,PROJECT,lower,upper,scalar_checks,ledger
from review import bump_review,moment_pressure_review

STATUS='repair-state-bound-and-conditional-frequency-acceptance'
FLAGS=dict(actual_repair_error_constants_bounded=True,
           positive_repair_neighborhood_bounded=True,
           frequency_acceptance_given_inherited_estimates=True,
           conditional_leading_cone_argument_completed=True,
           relies_on_inherited_analytic_construction=True,
           actual_high_frequency_profile_numerically_resolved=False,
           full_physical_stress_wave_construction_verified=False,
           full_PDE_corrections_verified=False,smooth_force_verified=False,
           blowup_verified=False,independent_foundation_review_completed=False,
           formal_proof_assistant_certificate=False)
RESOLVED=['actual_repair_error_constants']
REMAINING=['independent_review_of_inherited_continuum_estimates',
           'full_physical_residual_and_correction_budget',
           'convergence_of_the_PDE_corrections',
           'admissible_smooth_forcing_through_the_singular_time',
           'blowup_lower_bound_for_the_complete_field']


def expected_protocol():
    return dict(status=STATUS,parent_commit='89e00da7d41e8ee0dfb77548e7c299308e73eb55',
                input_envelope='A=C^4096',
                coefficient_norm='max_j max(sup|c_j|,sup|c_j_eta|)',
                field_angular_norm='sup|f|+sup|f_eta|',
                patch='X=rhop*exp(y), 0<y<5',bump_width='0.1',
                lambda_domain='0<lambda<=0.01',Ccorr='A^128',
                coefficient_radius='A^(-128)',selected_N='1+floor(H^32)',
                acceptance_conditional_on_inherited_estimates=True,
                pressure_state_error_order=0,moment_angular_error_order=1,
                interval_digits=[80,110],diagnostic_digits=70)


def provenance():
    """Check both the last argument and the earlier hashes it depended on."""
    record=json.loads((HERE.parent/'loop_modulation_bounds'/'evidence.json').read_text())
    hashes=record['source_hashes']
    previous=json.loads((HERE.parent/'compact_jet_envelope'/'evidence.json').read_text())
    hashes={**previous['source_hashes'],**hashes}
    if not hashes:
        raise ValueError('Missing inherited provenance')
    return {p:hashlib.sha256((PROJECT/p).read_bytes()).hexdigest()==digest
            for p,digest in hashes.items()}


def hashes():
    paths=[HERE/p for p in ('README.md','protocol.json','requirements.txt','bounds.py',
                           'review.py','audit.py','test_audit.py')]
    for directory in ('loop_modulation_bounds','compact_jet_envelope','stress_realization_audit'):
        paths += [HERE.parent/directory/p for p in ('README.md','evidence.json','bounds.py')]
    paths += [HERE.parent/'stress_realization_audit'/'repair.py']
    return {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def encode(v):
    if hasattr(v,'_mpi_'):
        return dict(lower=mp.nstr(lower(v),120),upper=mp.nstr(upper(v),120))
    if isinstance(v,mp.mpf): return mp.nstr(v,65)
    if isinstance(v,Fraction): return str(v)
    if isinstance(v,dict): return {k:encode(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)): return [encode(x) for x in v]
    return v


def generate():
    protocol=json.loads((HERE/'protocol.json').read_text())
    if protocol!=expected_protocol():
        raise ValueError('Norms, neighborhood, frequency or scientific scope changed')
    scalars=[scalar_checks(d) for d in protocol['interval_digits']]
    reviews=dict(bumps=bump_review(70),moments_and_pressure=moment_pressure_review(70))
    old=provenance()
    arithmetic=ledger()
    prior=json.loads((HERE.parent/'loop_modulation_bounds'/'evidence.json').read_text())
    gates=dict(scalar_and_frequency_checks=all(all(r['checks'].values()) for r in scalars),
               bump_derivative_and_neighborhood_controls=all(reviews['bumps']['checks'].values()),
               original_moment_and_pressure_equations=all(reviews['moments_and_pressure']['checks'].values()),
               norm_conversion_and_state_arithmetic=all(arithmetic['checks'].values()),
               inherited_bytes_preserved=all(old.values()),
               preceding_repair_obligation_matches=set(prior['remaining_obligations'])==set(RESOLVED),
               manufactured_diagnostics_explicit=all(r['manufactured'] for r in reviews.values()),
               physical_PDE_and_blowup_unverified=not any(FLAGS[k] for k in
                   ('full_physical_stress_wave_construction_verified','full_PDE_corrections_verified',
                    'smooth_force_verified','blowup_verified','independent_foundation_review_completed')))
    return dict(status=STATUS,flags=FLAGS,protocol=protocol,source_hashes=hashes(),
                scalar_bounds=scalars,error_ledger=arithmetic,reviews=reviews,
                inherited_provenance=old,resolved_obligations=RESOLVED,
                remaining_obligations=REMAINING,gates=gates)


def validate(record):
    if record.get('status')!=STATUS or record.get('flags')!=FLAGS:
        raise ValueError('Conditional leading-profile acceptance cannot promote the PDE or blow-up')
    if record.get('resolved_obligations')!=RESOLVED or record.get('remaining_obligations')!=REMAINING:
        raise ValueError('The proof obligations changed')
    if not record.get('gates') or not all(record['gates'].values()):
        raise ValueError('A mathematical constant, equation or provenance check failed')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify-record',type=Path)
    args=parser.parse_args()
    result=encode(generate())
    validate(result)
    if args.verify_record and json.loads(args.verify_record.read_text())!=result:
        raise SystemExit('Frozen record differs; no output written')
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(STATUS)
    print('Repair neighborhood and finite N accepted conditional on inherited analytic estimates.')
    print('Full physical corrections, smooth forcing, independent review and blow-up remain unverified.')


if __name__=='__main__': main()
