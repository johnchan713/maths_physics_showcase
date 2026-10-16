#!/usr/bin/env python3
"""Scope, exact-cancellation, missing-pressure and frozen-budget regressions."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import mpmath as mp
from audit import HERE, compare_records, generate, main, validate
from continuation import budget, taper


class ContinuationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record=generate()

    def test_original_source_ODEs_with_pressure(self):
        r=self.record['reviews']['original_pressure_source_and_final_segment']
        self.assertLess(r['maximum_original_ODE_error'],1e-65)
        self.assertTrue(r['checks']['Hc_uses_local_U_not_radial_average'])

    def test_pressure_value_term_not_deleted(self):
        self.assertTrue(self.record['reviews']['original_pressure_source_and_final_segment']['checks']['pressure_value_term_cannot_be_dropped'])

    def test_pressure_eta_term_not_deleted(self):
        self.assertTrue(self.record['reviews']['original_pressure_source_and_final_segment']['checks']['pressure_eta_term_cannot_be_dropped'])

    def test_log_slope_keeps_half_offset(self):
        self.assertTrue(self.record['reviews']['original_pressure_source_and_final_segment']['checks']['ell_is_half_plus_E_logarithmic_slope'])

    def test_exact_taper_errors_and_common_kappa(self):
        self.assertTrue(all(self.record['reviews']['exact_taper_identities']['checks'].values()))

    def test_beta_and_beta_squared_are_distinct(self):
        r=taper(3,7,Q(1,10**40),Q(1,3),1)
        self.assertEqual(r['main_Pc'],3+Q(49,9))
        self.assertEqual(r['main_v'],Q(1,10**40)*(3+Q(49,27)))

    def test_zero_axial_endpoint_has_no_shear_ratio_defect(self):
        r=taper(3,10**30,Q(1,10**100),0,Q(3,2))
        self.assertEqual(r['t'],0)
        self.assertEqual(r['Pc'],3)
        self.assertEqual(r['v'],Q(3,10**100))

    def test_wrong_denominators_and_taper_range_rejected(self):
        for args in ((0,1,1,1,1),(3,1,0,1,1),(3,1,1,-1,1),(3,1,1,2,1),(3,1,1,1,0)):
            with self.assertRaises(ValueError):
                taper(*args)

    def test_source_absorption_and_positive_ODE_comparisons(self):
        self.assertTrue(all(self.record['reviews']['positive_source_ODE_comparisons']['checks'].values()))

    def test_source_budgets_outward_at_frozen_scales(self):
        r=self.record['reviews']['outward_frozen_budgets']
        self.assertTrue(all(r['checks'].values()))
        self.assertEqual([v['scale'] for v in r['rows']],['P16','Md4','Md64'])

    def test_each_refined_interval_gate_passes(self):
        for r in self.record['reviews']['outward_frozen_budgets']['rows']:
            self.assertTrue(all(r['checks'].values()))

    def test_parameter_precision_refinement(self):
        self.assertLess(self.record['reviews']['outward_frozen_budgets']['maximum_relative_refinement_error'],1e-55)
        with self.assertRaises(ValueError):
            budget('Md64',80)

    def test_small_tail_is_enclosed_not_deleted(self):
        self.assertEqual(self.record['protocol']['log_sum_negligible_tail_interval'],['0','1e-200'])
        self.assertEqual(self.record['protocol']['wide_log_sum_enclosure'],'max(lower)..max(upper)+log(n)')
        for r in self.record['reviews']['outward_frozen_budgets']['rows']:
            self.assertGreater(r['enclosures']['wide_intervals'],0)
            self.assertGreater(r['enclosures']['negligible_positive_tails'],0)

    def test_final_constant_shear_error_not_width_small(self):
        r=self.record['reviews']['original_pressure_source_and_final_segment']
        self.assertLess(r['post_taper_fixed_ratio_log_error'],-.03)
        self.assertGreater(r['error_divided_by_width'][1],r['error_divided_by_width'][0]*10**8)
        self.assertFalse(self.record['flags']['width_small_field_error_extended_after_final_taper'])

    def test_actual_continuation_margin_scope(self):
        self.assertEqual(self.record['protocol']['actual_v_upper'],'4/5')
        self.assertEqual(self.record['protocol']['actual_Pc_lower'],'11/5')
        self.assertEqual(self.record['protocol']['actual_a_lower'],'C^-22')
        self.assertTrue(self.record['flags']['actual_reduced_continuation_and_two_taper_relaxed_margins_derived'])

    def test_global_pressure_pulse_and_residual_not_promoted(self):
        for k in ('whole_compact_shape_join_outer_foundation_independently_verified',
                  'completed_global_backward_pressure_normalization_verified',
                  'unconditional_uniform_actual_pulse_estimate_verified',
                  'unconditional_complete_residual_improvement_verified','blowup_verified'):
            r=deepcopy(self.record)
            r['flags'][k]=True
            with self.assertRaises(ValueError):
                validate(r)

    def test_omitted_obligation_or_pressure_check_rejected(self):
        r=deepcopy(self.record)
        r['remaining_obligations'].pop()
        with self.assertRaises(ValueError):
            validate(r)
        r=deepcopy(self.record)
        r['reviews']['original_pressure_source_and_final_segment']['checks']['pressure_eta_term_cannot_be_dropped']=False
        with self.assertRaises(ValueError):
            validate(r)

    def test_evidence_replay_and_source_hashes(self):
        validate(self.record)
        compare_records(self.record,json.loads((HERE/'evidence.json').read_text()))
        self.assertEqual(len(self.record['inherited_provenance']),497)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_nonfinite_and_exact_record_tampering_rejected(self):
        r=deepcopy(self.record)
        r['reviews']['original_pressure_source_and_final_segment']['maximum_original_ODE_error']=float('nan')
        with self.assertRaises(ValueError):
            compare_records(self.record,r)
        r=deepcopy(self.record)
        r['protocol']['width_and_kappa0']='C^-21'
        with self.assertRaises(ValueError):
            validate(r)

    def test_output_overwrite_protected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json'
            p.write_text('preserve')
            with self.assertRaises(SystemExit):
                main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__=='__main__':
    unittest.main()
