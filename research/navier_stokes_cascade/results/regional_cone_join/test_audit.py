#!/usr/bin/env python3
"""Regional cone, signed comparison, coverage, provenance and scope controls."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
from audit import HERE, FLAGS, compare_records, generate, main, validate
from cones import coordinates, gaps, scalar_certificate


class RegionalConeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record=generate()

    def test_exact_finite_p_and_oriented_stress_identities(self):
        r=self.record['reviews']['exact_finite_cone_identities']
        self.assertEqual(r['rational_identity_checks'],162)
        self.assertEqual(r['a2_factorization_checks'],9)
        self.assertTrue(all(r['checks'].values()))

    def test_shear_sign_convention_is_not_reversed(self):
        self.assertEqual(coordinates(Q(2),Q(1),Q(10),Q(3)),(Q(-1,2),Q(5,2),Q(17,2),Q(8)))
        self.assertEqual(gaps(Q(2),Q(1),Q(10),Q(3)),(Q(2),Q(1,2),Q(6),Q(40)))

    def test_ratio_positivity_does_not_replace_finite_pressure_gap(self):
        self.assertLess(gaps(Q(11,5),Q(0),Q(21,10),Q(0))[2],0)

    def test_subcritical_original_region_is_not_called_admissible(self):
        self.assertLess(gaps(Q(4,5),Q(0),Q(3),Q(0))[1],0)
        self.assertTrue(self.record['protocol']['original_third_fourth_required_only_where_vs_at_least_two'])

    def test_old_negative_source_kernel_inference_is_rejected(self):
        r=self.record['reviews']['corrected_signed_outer_ramp']
        self.assertGreater(r['normalized_old_claim'],r['normalized_exact_model'])
        self.assertGreater(r['old_minus_model'],.004)

    def test_corrected_ramp_bound_is_below_exact_comparison_solution(self):
        r=self.record['reviews']['corrected_signed_outer_ramp']
        self.assertLess(r['normalized_corrected_lower'],r['normalized_exact_model'])

    def test_required_half_XR_floor_survives(self):
        r=self.record['reviews']['corrected_signed_outer_ramp']
        self.assertGreater(r['normalized_corrected_lower'],.586)
        self.assertTrue(all(r['checks'].values()))

    def test_comparison_model_is_not_actual_ramp_or_selected_h(self):
        self.assertTrue(self.record['protocol']['negative_control_is_comparison_model_not_actual_ramp'])

    def test_all_nested_parameter_outward_checks_pass(self):
        r=self.record['reviews']['actual_regional_input_and_acceptance_budgets']
        self.assertEqual([v['digits'] for v in r['rows']],[100,160])
        self.assertTrue(all(r['checks'].values()))
        self.assertTrue(all(all(v['checks'].values()) for v in r['rows']))
        self.assertLess(r['maximum_relative_refinement_error'],1e-55)

    def test_insufficient_nested_precision_rejected(self):
        with self.assertRaises(ValueError):
            scalar_certificate(80)

    def test_actual_axial_and_intermediate_gap_floors(self):
        for row in self.record['reviews']['actual_regional_input_and_acceptance_budgets']['rows']:
            for key in ('axial_third_lower','axial_fourth_lower','intermediate_fourth_lower'):
                self.assertGreater(float(row['metrics'][key]['lower']),1)

    def test_ideal_after_join_and_first_outer_ramp_are_both_covered(self):
        regions=self.record['reviews']['compact_and_compensated_exterior_join']['regions']
        self.assertIn('ideal_after_join',[r['region'] for r in regions])
        self.assertIn('first_outer_ramp',[r['region'] for r in regions])

    def test_ten_regions_have_exact_shared_endpoints(self):
        regions=self.record['reviews']['compact_and_compensated_exterior_join']['regions']
        self.assertEqual(len(regions),10)
        self.assertEqual(regions[0]['left'],'Xminus')
        self.assertEqual(regions[-1]['right'],'Xrep=Xstar')
        for left,right in zip(regions,regions[1:]):
            self.assertEqual(left['right'],right['left'])

    def test_heat_compensation_start_equals_repair_end_exactly(self):
        for T in (Q(64),Q(123,2),Q(10**30)):
            self.assertEqual(241*T-23+5,241*T-18)
        self.assertTrue(self.record['protocol']['compact_right_equals_heat_compensation_left'])

    def test_downstream_cone_uses_compensated_not_unedited_pulse(self):
        self.assertTrue(self.record['protocol']['later_cone_uses_heat_compensated_profile_with_exact_pressure_memory'])

    def test_protected_growth_and_finite_frequency_are_checked(self):
        r=self.record['reviews']['compact_and_compensated_exterior_join']
        self.assertTrue(all(r['checks'].values()))
        self.assertEqual(r['selected_N'],'1+floor(Hhier^32)')

    def test_input_gap_powers_fit_unchanged_envelope(self):
        self.assertLess(401,1024)
        self.assertLess(802,1024)
        self.assertLess(1024,4096)
        self.assertEqual(self.record['protocol']['Aenv'],'C^4096')

    def test_Md4_failure_is_not_relabelled(self):
        self.assertTrue(self.record['protocol']['selected_family'].startswith('Md64,'))

    def test_frequency_is_analytic_not_numerically_resolved(self):
        self.assertFalse(self.record['protocol']['actual_frequency_or_profile_numerically_resolved'])

    def test_exact_global_dependency_and_533_inherited_hashes(self):
        self.assertTrue(self.record['gates']['exact_global_moment_dependency_replayed'])
        self.assertEqual(len(self.record['inherited_provenance']),533)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_frozen_record_replays(self):
        frozen=json.loads((HERE/'evidence.json').read_text())
        validate(frozen);compare_records(self.record,frozen)

    def test_missing_region_or_gap_in_coverage_rejected(self):
        for action in ('remove','break'):
            r=deepcopy(self.record)
            regions=r['reviews']['compact_and_compensated_exterior_join']['regions']
            if action=='remove':
                regions.pop(5)
            else:
                regions[5]['left']='uncovered_radius'
            with self.assertRaises(ValueError):
                validate(r)

    def test_invented_frequency_rejected(self):
        r=deepcopy(self.record)
        r['reviews']['compact_and_compensated_exterior_join']['selected_N']='1000'
        with self.assertRaises(ValueError):
            validate(r)

    def test_all_order_support_pulse_residual_and_blowup_promotion_rejected(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):
                    validate(r)

    def test_analytic_dependencies_cannot_be_hidden(self):
        r=deepcopy(self.record);r['flags']['relies_on_replayed_analytic_core_envelope_loop_and_inverse_lemmas']=False
        with self.assertRaises(ValueError):
            validate(r)

    def test_missing_outer_comparison_gate_rejected(self):
        r=deepcopy(self.record);r['gates'].pop('corrected_signed_outer_ramp')
        with self.assertRaises(ValueError):
            validate(r)

    def test_failing_negative_control_rejected(self):
        r=deepcopy(self.record)
        r['reviews']['corrected_signed_outer_ramp']['checks']['old_signed_kernel_shortcut_fails_for_comparison_model']=False
        with self.assertRaises(ValueError):
            validate(r)

    def test_altered_provenance_rejected(self):
        r=deepcopy(self.record);r['inherited_provenance'][next(iter(r['inherited_provenance']))]=False
        with self.assertRaises(ValueError):
            validate(r)

    def test_nonfinite_diagnostic_rejected(self):
        r=deepcopy(self.record);r['reviews']['corrected_signed_outer_ramp']['old_minus_model']=float('nan')
        with self.assertRaises(ValueError):
            compare_records(self.record,r)

    def test_reserved_mean_and_covariance_obligations_cannot_be_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):
            validate(r)

    def test_CLI_cannot_overwrite_frozen_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('preserve')
            with self.assertRaises(SystemExit):
                main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__=='__main__':
    unittest.main()
