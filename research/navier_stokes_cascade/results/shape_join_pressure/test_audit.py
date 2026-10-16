#!/usr/bin/env python3
"""C2 input, exact pressure preservation, and scientific-scope regressions."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
from audit import HERE, compare_records, generate, main, validate
from shape import budget,Jet2


class ShapeJoinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record=generate()

    def test_shape_source_barrier_has_rational_slack(self):
        self.assertTrue(self.record['reviews']['shape_source_and_C2_root_scalars']['checks']['source_shape_margin'])

    def test_source_barrier_crossing_sign(self):
        self.assertTrue(self.record['reviews']['shape_source_and_C2_root_scalars']['checks']['barrier_crossing_points_upward'])

    def test_shape_logarithmic_slope_and_width(self):
        self.assertTrue(self.record['reviews']['shape_source_and_C2_root_scalars']['checks']['shape_slope_budget_at_B10'])

    def test_frozen_parameter_bounds_pass_at_both_precisions(self):
        r=self.record['reviews']['outward_actual_shape_join_and_Cp_root_budgets']
        self.assertTrue(all(r['checks'].values()))
        self.assertLess(r['maximum_relative_refinement_error'],1e-55)
        for row in r['rows']:
            self.assertTrue(all(row['checks'].values()))

    def test_actual_C2_entry_below_p255epsilon(self):
        for row in self.record['reviews']['outward_actual_shape_join_and_Cp_root_budgets']['rows']:
            self.assertLess(row['entry_upper'],2.55e-17)

    def test_positive_bump_inverse_not_numerical_newton_certificate(self):
        for row in self.record['reviews']['outward_actual_shape_join_and_Cp_root_budgets']['rows']:
            self.assertLess(row['inverse_upper'],1000)
            self.assertTrue(row['checks']['C2_quadratic_bound_below100'])

    def test_same_C2_root_and_factor_two_for_raw_curvature(self):
        r=self.record['reviews']['shape_source_and_C2_root_scalars']
        self.assertEqual(Q(r['coefficient_C2_bound']),Q(572,10**16))
        self.assertTrue(r['checks']['same_root_inside_C1_ball'])
        self.assertTrue(r['checks']['raw_second_derivative_is_twice_factorial_bound'])

    def test_pressure_and_energy_normalizations_keep_C2_jets(self):
        r=self.record['reviews']['complete_C2_normalization_and_energy_cancellation']
        self.assertTrue(all(r['checks'].values()))
        self.assertLess(r['maximum_normalized_jet_error'],1e-65)

    def test_inverse_f_squared_norm_is_20_not_12(self):
        e=Jet2(Q(1),Q(1),Q(0));r=(1+e*e)*(1+e*e)
        self.assertEqual(abs(r.v)+abs(r.d)+abs(r.dd)/2,20)

    def test_complete_pressure_ledger_identity(self):
        self.assertTrue(self.record['reviews']['exact_pressure_ledger_and_renormalization_controls']['checks']['forward_equals_backward_when_exact_total_restored'])

    def test_quadratic_pressure_defect_cannot_be_deleted(self):
        r=self.record['reviews']['exact_pressure_ledger_and_renormalization_controls']
        self.assertTrue(r['checks']['quadratic_pressure_remainder_nonzero'])
        self.assertGreater(Q(r['quadratic_pressure_residual']),0)

    def test_one_angle_pressure_match_does_not_restore_derivative(self):
        self.assertTrue(self.record['reviews']['exact_pressure_ledger_and_renormalization_controls']['checks']['one_angle_zero_not_functional_zero'])

    def test_raw_I_is_not_a_finite_zero_total(self):
        r=self.record['reviews']['exact_pressure_ledger_and_renormalization_controls']
        self.assertEqual(r['toy_raw_I'],['2','4','8','16'])
        self.assertTrue(r['checks']['renormalized_I_can_vanish'])

    def test_outer_angular_and_heat_Cp_roots_are_explicit_dependencies(self):
        self.assertTrue(self.record['flags']['exact_outer_Cp_correction_maps_are_explicit_dependencies'])
        for row in self.record['reviews']['outward_actual_shape_join_and_Cp_root_budgets']['rows']:
            self.assertTrue(row['checks']['exact_angular_Cp_root_contraction'])
            self.assertTrue(row['checks']['exact_heat_Cp_root_contraction'])

    def test_no_giant_actual_profile_or_Md4_cone_promotion(self):
        self.assertFalse(self.record['protocol']['manufactured_fixtures_are_actual_profile'])
        self.assertTrue(self.record['reviews']['outward_actual_shape_join_and_Cp_root_budgets']['checks']['outer_Md4_cone_failure_not_relabelled'])

    def test_insufficient_precision_rejected(self):
        with self.assertRaises(ValueError):
            budget('Md64',80)

    def test_scope_cannot_be_promoted_to_all_moments_pulses_or_blowup(self):
        for key in ('all_global_moment_and_signed_stress_identities_verified','whole_leading_profile_independently_verified',
                    'coefficient_and_reserved_mean_supports_verified','unconditional_uniform_actual_pulse_estimate_verified',
                    'unconditional_complete_residual_improvement_verified','blowup_verified'):
            r=deepcopy(self.record);r['flags'][key]=True
            with self.assertRaises(ValueError):
                validate(r)

    def test_pressure_scope_is_leading_not_full_corrected_pressure(self):
        self.assertIn('not full corrected pressure',self.record['protocol']['backward_pressure_scope'])
        self.assertTrue(self.record['protocol']['compact_modulation_input_hypotheses_remain_dependencies'])

    def test_omitted_global_obligation_rejected(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop()
        with self.assertRaises(ValueError):
            validate(r)

    def test_omitted_pressure_gate_rejected(self):
        r=deepcopy(self.record)
        r['reviews']['exact_pressure_ledger_and_renormalization_controls']['checks']['quadratic_pressure_remainder_nonzero']=False
        with self.assertRaises(ValueError):
            validate(r)

    def test_frozen_record_and_508_inherited_hashes(self):
        validate(self.record);compare_records(self.record,json.loads((HERE/'evidence.json').read_text()))
        self.assertEqual(len(self.record['inherited_provenance']),508)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_nonfinite_and_changed_exact_records_rejected(self):
        r=deepcopy(self.record)
        r['reviews']['complete_C2_normalization_and_energy_cancellation']['maximum_normalized_jet_error']=float('nan')
        with self.assertRaises(ValueError):
            compare_records(self.record,r)
        r=deepcopy(self.record);r['protocol']['inverse_f_squared_C2_norm']=12
        with self.assertRaises(ValueError):
            validate(r)

    def test_CLI_cannot_overwrite_existing_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('preserve')
            with self.assertRaises(SystemExit):
                main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__=='__main__':
    unittest.main()
