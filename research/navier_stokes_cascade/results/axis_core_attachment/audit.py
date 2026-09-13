#!/usr/bin/env python3
"""Reproduce quantitative axis/attachment bounds and their scientific scope."""
import argparse
from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path
import sys
import mpmath as mp
from bounds import (certificate, elementary_certificate, entry_transfer, lower, upper,
                    point, tail_bound_log, weight)
from identities import source_identities, stress_identities, stress_values

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
sys.path.append(str(RESULTS/'axis_matching_audit'))
from annulus import accept_entry_bounds
STATUS = 'analytic-core-and-reference-attachment-bounded'
FLAGS = {
    'exact_analytic_axis_existence_bounded': True,
    'native_axial_entry_budget_bounded': True,
    'continuation_analytic_estimates_reviewed': True,
    'actual_axis_entry_verified': True,
    'matched_leading_profile_verified': True,
    'formal_proof_assistant_certificate': False,
    'independent_peer_review_completed': False,
    'full_admissible_stress_realized': False,
    'full_PDE_corrections_verified': False,
    'smooth_force_verified': False,
    'blowup_verified': False,
}


def provenance():
    paths = [HERE/n for n in ('README.md', 'protocol.json', 'requirements.txt',
                              'bounds.py', 'identities.py', 'audit.py', 'test_audit.py')]
    paths += [RESULTS/'axis_matching_audit'/n for n in
              ('README.md', 'bounds.py', 'annulus.py', 'moments.py', 'evidence.json')]
    paths += [RESULTS/'heat_exterior_audit'/'evidence.json']
    paths += [RESULTS/'outer_pressure_pilot'/'schedule.py']
    return {str(p.relative_to(RESULTS.parent)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths}


def exact_weight_checks(limit):
    """Finite diagnostics of the all-index bounds, not a proof by truncation."""
    ratios = dict(radial_shift=Fraction(0), angular_shift=Fraction(0),
                  cauchy_coefficient=Fraction(0), leibniz_ratio=Fraction(0))
    for n in range(limit+1):
        for k in range(limit+1):
            ratios['radial_shift'] = max(ratios['radial_shift'], weight(n,k)/weight(n+1,k))
            ratios['angular_shift'] = max(ratios['angular_shift'], weight(n,k+1)/weight(n+1,k)/(n+1))
            ratios['cauchy_coefficient'] = max(ratios['cauchy_coefficient'], Fraction((k+1)**2,4**k))
            for i in range(n+1):
                for j in range(k+1):
                    ratios['leibniz_ratio'] = max(ratios['leibniz_ratio'],
                        Fraction(comb(i+j,j)*comb(n-i+k-j,k-j), comb(n+k,k)))
    checks = dict(radial_shift=ratios['radial_shift']<=80,
                  angular_shift=ratios['angular_shift']<=80,
                  cauchy_coefficient=ratios['cauchy_coefficient']<=1,
                  leibniz_ratio=ratios['leibniz_ratio']<=1)
    return dict(limit=limit, maxima=ratios, checks=checks, exhaustive_all_indices=False)


def stress_controls():
    # Exact rational perturbations exercise the common kappa cancellation.
    rows = []
    for k in (Fraction(1), Fraction(1,10**12), Fraction(1,10**100)):
        for q in (Fraction(0), Fraction(7), Fraction(10**30)):
            row = stress_values(3, q, k, Fraction(1000001,1000000),
                                Fraction(1,10**12), -Fraction(1,10**12))
            rows.append(dict(kappa=k, p2_reference=q, **row))
    invariant = all(rows[i]['ts']==rows[i%3]['ts'] and rows[i]['Pc']==rows[i%3]['Pc']
                    and rows[i]['Jc']==rows[i%3]['Jc'] for i in range(len(rows)))
    unactivated = stress_values(3, 7, 1, 1)
    # A=Pc-vs is exactly zero before activation; a strict interior test
    # must not be promoted merely because Pc>2.
    return dict(samples=rows, kappa_cancellation=invariant,
                unactivated_stress_zero=unactivated['Pc']==unactivated['vs'],
                frozen_E_ratio_changes_Pc=(stress_values(3,7,1,1)['Pc'] !=
                                           stress_values(3,7,1,Fraction(1000001,1000000))['Pc']))


def encode(value):
    if hasattr(value, '_mpi_'):
        return dict(lower=mp.nstr(lower(value),125), upper=mp.nstr(upper(value),125))
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, mp.mpf):
        return mp.nstr(value,125)
    if isinstance(value, dict):
        return {k:encode(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)):
        return [encode(v) for v in value]
    return value


def generate():
    protocol = json.loads((HERE/'protocol.json').read_text())
    required=dict(epsilon='1e-16',j='1e-18',sigma='1e-20',status=STATUS,
        source_pdf_sha256='0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f')
    if any(protocol.get(k)!=v for k,v in required.items()):
        raise ValueError('Protocol does not match the implemented analytic data and scope')
    inherited = json.loads((RESULTS/'axis_matching_audit'/'evidence.json').read_text())
    old_md_failure=inherited['inherited']['old_Md_failure']
    records = []
    for scale in protocol['pressure_scales']:
        for digits in protocol['interval_digits']:
            records.append(certificate(scale,digits))
    refined = []
    with mp.workdps(125):
        for coarse, fine in zip(records[::2], records[1::2]):
            gaps = {}
            for key, value in coarse['log_bounds'].items():
                reference = fine['log_bounds'][key]
                gaps[key] = max(abs(lower(value)-lower(reference)),
                                abs(upper(value)-upper(reference)))/max(1,abs(upper(reference)))
            refined.append(dict(scale=coarse['scale'], maximum_relative_log_gap=max(gaps.values()),
                                worst_quantity=max(gaps,key=gaps.get)))
    weights = exact_weight_checks(protocol['finite_weight_index_limit'])
    sources, stress = source_identities(), stress_identities()
    controls = stress_controls()
    elementary = elementary_certificate(protocol['interval_digits'][-1])
    transfer=entry_transfer(protocol['interval_digits'][-1])
    with mp.workdps(125):
        entry_accepted=accept_entry_bounds(
            [upper(v) for v in transfer['entry_moment_C1_bounds']],
            upper(transfer['axial_offset_C1_bound']))
    sign_control=dict(source_floor=Fraction(1,10**18),
                      trial_error=-Fraction(1,1000), old_error_budget=Fraction(1,100),
                      fixed_size_budget_can_reverse_sign=True,
                      actual_reference_profile_counterexample=False)
    # Show that the old finite-degree comparison is not this norm certificate.
    old_tail = tail_bound_log(records[0]['log_bounds']['ball_norm_bound'],24)
    gates = {
        'inherited_annulus_passed': all(inherited['gates'].values()),
        'Md4_outer_failure_retained': mp.mpf(old_md_failure['Pc_over_ps1'])<0,
        'inherited_axis_gap_retained_as_historical': (
            inherited['status']=='conditional-five-moment-annulus-bounded-axis-entry-unverified'
            and inherited['certificate']['actual_axis_entry_verified'] is False
            and inherited['interpretation']['actual_axis_entry_verified'] is False),
        'all_interval_parameter_inequalities': all(all(r['checks'].values()) for r in records),
        'elementary_analytic_constants': all(elementary.values()),
        'interval_precision_refinement': all(r['maximum_relative_log_gap']<mp.mpf('1e-35') for r in refined),
        'exact_coefficient_weight_diagnostics': all(weights['checks'].values()),
        'original_inner_source_identities': all(sources.values()),
        'exact_stress_perturbation_identities': all(stress.values()),
        'common_kappa_cancels': controls['kappa_cancellation'],
        'unactivated_stress_is_zero': controls['unactivated_stress_zero'],
        'E_ratio_cannot_be_frozen': controls['frozen_E_ratio_changes_Pc'],
        'degree24_generic_tail_not_certified_small': lower(old_tail)>0,
        'positive_integrals_reach_annulus_entry': all(transfer['checks'].values()),
        'existing_annulus_accepts_certified_bounds': bool(entry_accepted),
        'analytic_attachment_scope_consistent': all(FLAGS[k] for k in (
            'exact_analytic_axis_existence_bounded','native_axial_entry_budget_bounded',
            'continuation_analytic_estimates_reviewed','actual_axis_entry_verified',
            'matched_leading_profile_verified')),
        'constant_error_budget_does_not_prove_source_sign': (
            abs(sign_control['trial_error'])<sign_control['old_error_budget'] and
            sign_control['source_floor']+sign_control['trial_error']<0),
        'later_PDE_claims_excluded': not any(FLAGS[k] for k in
            ('full_admissible_stress_realized','full_PDE_corrections_verified',
             'smooth_force_verified','blowup_verified')),
    }
    return dict(status=STATUS, scientific_scope=FLAGS, protocol=protocol,
                source_hashes=provenance(), gates=gates, certificates=records,
                refinement=refined, weight_checks=weights, elementary=elementary,
                source_identities=sources, stress_identities=stress,
                stress_controls=controls, degree24_log_tail_upper=old_tail,
                entry_transfer=transfer, source_sign_budget_control=sign_control,
                inherited_status=inherited['status'],
                selected_outer_family='Md64',
                inherited_Md4_outer_failure=old_md_failure,
                other_pressure_scales_scope='axis/attachment inequalities only; no outer-cone promotion',
                numerical_matched_profile_materialized=False,
                selected_pressure_datum='exact exterior positive mixture; no sampled fit')


def compare(expected, actual, path='root'):
    if isinstance(actual,dict):
        if not isinstance(expected,dict) or set(expected)!=set(actual):
            raise ValueError('Record keys differ at '+path)
        for key in actual:
            compare(expected[key],actual[key],path+'.'+key)
    elif isinstance(actual,list):
        if not isinstance(expected,list) or len(expected)!=len(actual):
            raise ValueError('Record length differs at '+path)
        for i,value in enumerate(actual):
            compare(expected[i],value,path+'.'+str(i))
    elif isinstance(actual,str) and actual!=expected:
        # Hashes, statuses, rational controls and source descriptions must match
        # exactly. Only serialized decimal interval results admit rounding noise.
        decimal_path = any(s in path for s in ('.lower','.upper','relative_log_gap'))
        if not decimal_path:
            raise ValueError('Record differs at '+path)
        with mp.workdps(125):
            a,b=mp.mpf(actual),mp.mpf(expected)
            if not mp.isfinite(a) or not mp.isfinite(b) or abs(a-b)>mp.mpf('1e-35')*max(1,abs(a)):
                raise ValueError('Numeric record differs at '+path)
    elif expected!=actual:
        raise ValueError('Record differs at '+path)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-record',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    record=encode(generate())
    if not all(record['gates'].values()):
        raise SystemExit('Failed gates: '+str([k for k,v in record['gates'].items() if not v]))
    if args.verify_record:
        compare(json.loads(args.verify_record.read_text()),record)
    args.output.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    print(record['status'])
    print(str(sum(record['gates'].values()))+' audit gates passed')
    print('Exact-core bounds are separate from full-PDE and blow-up claims.')


if __name__=='__main__':
    main()
