#!/usr/bin/env python3
"""Tests target derivative losses, implicit inverses, and inappropriate promotion."""
from copy import deepcopy
from fractions import Fraction as F
import json
import unittest
from bounds import HERE, scalar_checks, derivative_ledger, error_ledger, transformed_target_bound
from review import root_review,phase_review,modulation_review
from audit import STATUS,FLAGS,RESOLVED,INHERITED,REMAINING,validate,provenance,expected_protocol


class AnalyticAndDiagnosticTests(unittest.TestCase):
    def test_outward_constants_at_both_precisions(self):
        for digits in (80,110):
            r = scalar_checks(digits)
            self.assertTrue(all(r['checks'].values()),r['checks'])
            self.assertFalse(r['G_H_and_N_materialized'])

    def test_root_differentiation_at_zero_and_nonzero_pressure(self):
        r = root_review()
        self.assertTrue(all(r['checks'].values()),r['checks'])

    def test_moving_phase_chain_rule_and_failure_control(self):
        r = phase_review()
        self.assertTrue(all(r['checks'].values()),r['checks'])

    def test_actual_field_differentiation_of_manufactured_modulation(self):
        r = modulation_review()
        self.assertTrue(all(r['checks'].values()),r['checks'])

    def test_lambda_loss_survives_transformed_matrix_conditioning(self):
        small = transformed_target_bound(F(1,10**12),7,F(1,10**3))
        tiny = transformed_target_bound(F(1,10**12),7,F(1,10**30))
        self.assertEqual(tiny/small,10**27)
        self.assertGreater(tiny,1)

    def test_target_norm_controls_both_normalized_rows_before_subtraction(self):
        error, multiplier, lam = F(1,100),F(3),F(1,1000)
        # Opposite signed row errors attain the triangle bound exactly.
        self.assertEqual((multiplier*error-(-multiplier*error))/lam,
                         transformed_target_bound(error,multiplier,lam))

    def test_target_rejects_invalid_lambda_and_negative_norms(self):
        for args in ((1,1,0),(1,1,-1),(1,1,2),(-1,1,1),(1,0,1)):
            with self.assertRaises(ValueError):
                transformed_target_bound(*args)

    def test_derivative_budget_stops_at_state_values(self):
        ledger = error_ledger()
        self.assertEqual(ledger['angular_moment_error_order'],1)
        self.assertEqual(ledger['pressure_state_error_order'],0)
        self.assertEqual(ledger['p1_coefficient'],96)
        self.assertEqual(ledger['p2_coefficient'],196)
        self.assertTrue(derivative_ledger()['inverse_phase_derivative_included'])

    def test_inherited_source_bytes_and_protocol(self):
        self.assertTrue(all(provenance().values()))
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()),expected_protocol())


class ScopeTests(unittest.TestCase):
    def record(self):
        return deepcopy(dict(status=STATUS,flags=FLAGS,resolved_obligations=RESOLVED,
                             inherited_resolved_obligations=INHERITED,
                             remaining_obligations=REMAINING,gates={'valid':True}))

    def test_only_modulation_obligation_is_newly_closed(self):
        validate(self.record())
        self.assertEqual(RESOLVED,['actual_modulation_error_constants'])
        self.assertEqual(set(REMAINING),{'actual_repair_error_constants'})

    def test_repair_frequency_and_blowup_cannot_be_promoted(self):
        for key in ('actual_repair_error_constants_bounded',
                    'actual_joined_profile_frequency_selected',
                    'full_admissible_stress_realized','full_PDE_corrections_verified',
                    'smooth_force_verified','blowup_verified'):
            record = self.record()
            record['flags'][key] = True
            with self.assertRaises(ValueError):
                validate(record)

    def test_missing_obligation_and_failed_gate_rejected(self):
        record = self.record()
        record['remaining_obligations'] = {}
        with self.assertRaises(ValueError):
            validate(record)
        record = self.record()
        record['gates']['valid'] = False
        with self.assertRaises(ValueError):
            validate(record)


if __name__ == '__main__':
    unittest.main()
