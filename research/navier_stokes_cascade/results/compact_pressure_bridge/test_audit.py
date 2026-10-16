#!/usr/bin/env python3
"""Meaningful missing-term and scientific-scope controls."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import mpmath as mp
from audit import HERE, compare_records, generate, main, validate
from bridge import activation_budget, memory, pulse, sigma


class BridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = generate()

    def test_full_activation_not_only_first_quarter(self):
        rows = self.record['reviews']['positive_whole_interval_operator']['rows']
        self.assertEqual([r['s'] for r in rows],[.2,.5,.9,1.0])
        self.assertTrue(all(self.record['reviews']['positive_whole_interval_operator']['checks'].values()))

    def test_exact_scalar_growth_and_cone_budgets(self):
        self.assertTrue(all(self.record['reviews']['whole_activation_scalar_budget']['checks'].values()))

    def test_final_activation_need_only_relaxed_cone(self):
        self.assertLess(activation_budget(2)['exit_v_upper'],1)

    def test_out_of_range_activation_rejected(self):
        with self.assertRaises(ValueError):
            activation_budget(1)
        with self.assertRaises(ValueError):
            sigma('1.1')
        self.assertEqual(sigma(0),0)
        self.assertEqual(sigma(1),1)

    def test_exact_memory_against_original_forward_map(self):
        self.assertTrue(self.record['reviews']['exact_pressure_memory']['checks']['full_forward_subtraction_exact'])

    def test_four_rows_do_not_restore_pressure(self):
        zeros = [(Q(0),Q(0))]*4
        result = memory(Q(2),Q(1),Q(0),Q(2),Q(0),Q(1,100),zeros+[(Q(1,1000),Q(1,1000))])
        self.assertEqual(result['p1'],0)
        self.assertEqual(result['p2'],Q(-1,500))

    def test_pressure_memory_amplifies_even_on_equal_fields(self):
        self.assertTrue(self.record['reviews']['exact_pressure_memory']['checks']['pressure_memory_amplifies_with_radius'])

    def test_arbitrarily_small_pressure_memory_breaks_zero_stress_edge(self):
        self.assertTrue(self.record['reviews']['exact_pressure_memory']['checks']['pressure_error_survives_a_zero_stress_edge'])

    def test_nonzero_U_memory_term_retained(self):
        self.assertTrue(self.record['reviews']['exact_pressure_memory']['checks']['minus_U_delta_W_is_required'])

    def test_all_five_rows_with_derivatives_restore_exactly(self):
        result = memory(Q(2),Q(1),Q(4),Q(2),Q(1,3),Q(1,100),[(Q(0),Q(0))]*5)
        self.assertTrue(all(v == 0 for v in result.values()))

    def test_agreement_at_one_angle_does_not_restore_derivative(self):
        mu = [(Q(0),Q(0))]*4+[(Q(0),Q(1,1000))]
        self.assertNotEqual(memory(Q(2),Q(1),0,Q(2),0,Q(1,100),mu)['p2'],0)

    def test_both_pressure_terms_present(self):
        self.assertTrue(self.record['reviews']['exact_pressure_memory']['checks']['pressure_value_and_eta_terms_separately_exercised'])

    def test_rotated_stress_orientation_and_gap_identities(self):
        self.assertTrue(all(self.record['reviews']['compact_pulse_value_identities']['checks'].values()))

    def test_unoriented_or_nonpositive_fourth_gap_rejected(self):
        for args in ((0,0,4,0,1),(3,0,2,0,1),(3,0,6,100,1)):
            with self.assertRaises(ValueError):
                pulse(*map(Q,args))

    def test_negative_orientation_not_accepted_by_squared_ratio(self):
        with self.assertRaises(ValueError):
            pulse(Q(3),Q(0),Q(0),Q(0),Q(1))

    def test_global_promotion_rejected(self):
        for key in ('whole_leading_profile_independently_verified',
                    'post_activation_continuation_margins_independently_verified',
                    'unconditional_uniform_actual_pulse_estimate_verified',
                    'completed_global_backward_pressure_normalization_verified','blowup_verified'):
            r = deepcopy(self.record)
            r['flags'][key] = True
            with self.assertRaises(ValueError):
                validate(r)

    def test_missing_obligation_or_gate_rejected(self):
        r = deepcopy(self.record)
        r['remaining_obligations'].pop()
        with self.assertRaises(ValueError):
            validate(r)
        r = deepcopy(self.record)
        r['gates']['exact_pressure_memory'] = False
        with self.assertRaises(ValueError):
            validate(r)

    def test_inherited_bytes_preserved(self):
        self.assertGreater(len(self.record['inherited_provenance']),474)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_record_replay_and_output_overwrite_protection(self):
        frozen = json.loads((HERE/'evidence.json').read_text())
        validate(frozen)
        compare_records(self.record,frozen)
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'record.json'
            path.write_text('keep')
            with self.assertRaises(SystemExit):
                main(['--output',str(path)])
            self.assertEqual(path.read_text(),'keep')


if __name__ == '__main__':
    unittest.main()
