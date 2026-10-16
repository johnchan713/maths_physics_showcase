#!/usr/bin/env python3
"""Support, pressure memory, potential cutoff, provenance and scope controls."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
from audit import HERE,FLAGS,compare_records,generate,main,validate
from supports import support_edges,radial_flux,cut_potential_radial_flux,mean_pressure


class SupportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record=generate()

    def test_actual_unchanged_patch_geometry(self):
        r=self.record['reviews']['actual_reserved_geometry']
        self.assertTrue(all(r['checks'].values()))
        self.assertEqual(r['stage_intervals'],[[-25,-20],[-20,-15],[-14,-9],[-8,-3]])

    def test_five_disjoint_bumps_precede_velocity_edge_and_mean_patch(self):
        r=self.record['reviews']['actual_reserved_geometry']
        self.assertEqual(len(r['positive_bump_supports']),5)
        self.assertEqual(r['positive_bump_supports'][-1],['49/20','51/20'])
        self.assertLess(Q(51,20),3);self.assertLess(3,6)

    def test_first_pressure_is_not_supported_with_velocity(self):
        self.assertEqual(support_edges(1),dict(order=1,velocity='Xvel',pressure='Xv',stress='Xb'))

    def test_second_pressure_and_stress_keep_axial_source(self):
        self.assertEqual(support_edges(2),dict(order=2,velocity='Xvel',pressure='Xv',stress='Xv'))

    def test_arbitrary_later_order_has_sharper_common_support(self):
        for n in (3,4,100,10**6):
            self.assertEqual(set(support_edges(n).values())-{'Xvel'},{n})

    def test_invalid_formal_orders_rejected(self):
        for n in (0,-1,1.5,True,'2'):
            with self.assertRaises(ValueError):
                support_edges(n)

    def test_exact_rational_patch_and_curl_identities(self):
        r=self.record['reviews']['exact_support_and_cutoff_identities']
        self.assertEqual(r['rational_identity_checks'],90)
        self.assertTrue(all(r['checks'].values()))

    def test_zero_U_without_zero_M_is_not_zero_radial_velocity(self):
        self.assertNotEqual(radial_flux(Q(0),Q(1),Q(0),Q(10),Q(1,2),Q(1,200),1),0)

    def test_potential_cutoff_product_rule_sign_and_factor(self):
        self.assertEqual(cut_potential_radial_flux(Q(5,32),Q(-27,32),Q(0),Q(1),Q(1,2),Q(1,200)),Q(225,266))

    def test_zero_primitive_kills_cutoff_transition_at_any_order(self):
        self.assertEqual(cut_potential_radial_flux(Q(2,7),Q(-31,17),Q(0),Q(0),Q(-3,4),Q(1,200)),0)

    def test_first_pressure_changes_on_mean_patch(self):
        self.assertEqual(mean_pressure(1,Q(20),Q(3),Q(7))-mean_pressure(1,Q(10),Q(3),Q(7)),Q(9,80))

    def test_second_pressure_constant_is_not_forced_zero(self):
        self.assertEqual(mean_pressure(2,Q(20),Q(3),Q(7)),7)

    def test_higher_mean_pressure_is_zero_after_exact_total_moment(self):
        self.assertEqual(mean_pressure(3,Q(20),Q(3),Q(7)),0)

    def test_current_order_map_and_complete_pressure_target(self):
        r=self.record['reviews']['finite_formal_induction']
        self.assertEqual(r['formal_grade_checks'],36)
        self.assertTrue(all(r['checks'].values()))
        self.assertEqual(r['raw_normalization_costs'],['lambda^-1','e(eta)^-1'])

    def test_continuum_inverse_has_outward_strict_bounds(self):
        r=self.record['reviews']['continuum_linear_inverse']
        self.assertGreater(r['E_inverse_upper'],787.56)
        self.assertLess(r['E_inverse_upper'],1000)
        self.assertLess(r['U_inverse_upper'],2)
        self.assertTrue(all(r['checks'].values()))

    def test_shifted_viscosity_matches_distinct_full_series_assembly(self):
        r=self.record['reviews']['shifted_pressure_source_controls']
        self.assertLess(r['maximum_source_equation_error'],1e-11)
        self.assertTrue(r['checks']['inherited_recursion_matches_independent_shifted_operator'])

    def test_dropping_radial_axial_viscosity_is_detected(self):
        self.assertGreater(self.record['reviews']['shifted_pressure_source_controls']['maximum_omitted_viscosity_error'],1)

    def test_order3_and_above_have_no_pressure_source_after_velocity_support(self):
        self.assertEqual(self.record['reviews']['shifted_pressure_source_controls']['zero_high_order_source_checks'],60)

    def test_nonzero_second_pressure_control_is_labelled_manufactured(self):
        r=self.record['reviews']['shifted_pressure_source_controls']
        self.assertEqual(r['toy_second_pressure_plateau'],'153/10')
        self.assertTrue(r['checks']['manufactured_profiles_not_actual_selected_h_or_actual_coefficients'])
        self.assertTrue(self.record['protocol']['actual_second_pressure_constant_not_evaluated'])

    def test_actual_h_not_replaced_by_diagnostic_h(self):
        self.assertEqual(self.record['protocol']['manufactured_h'],.005)
        self.assertTrue(self.record['protocol']['fixed_parameters_and_all_historical_evidence_unchanged'])

    def test_finite_order_constants_not_uniform_in_order(self):
        self.assertTrue(self.record['protocol']['finite_coefficients_not_uniform_in_order'])
        self.assertFalse(FLAGS['full_weighted_background_cutoff_inputs_verified'])

    def test_parent_replayed_and_all_551_inherited_hashes_match(self):
        self.assertTrue(self.record['gates']['parent_leading_join_replayed'])
        self.assertEqual(len(self.record['inherited_provenance']),551)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_frozen_record_reproduces(self):
        frozen=json.loads((HERE/'evidence.json').read_text())
        validate(frozen);compare_records(self.record,frozen)

    def test_promoted_pulse_mean_residual_or_blowup_claim_rejected(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):
                    validate(r)

    def test_analytic_dependency_cannot_be_removed(self):
        r=deepcopy(self.record);r['flags']['relies_on_analytic_collar_and_continuum_inverse_lemmas']=False
        with self.assertRaises(ValueError):
            validate(r)

    def test_invented_pressure_zero_on_mean_patch_rejected(self):
        r=deepcopy(self.record);r['protocol']['mean_pressure_not_equal_to_leading_pressure']=False
        with self.assertRaises(ValueError):
            validate(r)

    def test_exceptional_second_pressure_edge_cannot_be_shortened(self):
        r=deepcopy(self.record)
        r['reviews']['exact_support_and_cutoff_identities']['support_edges'][1]['pressure']='Xvel'
        with self.assertRaises(ValueError):
            validate(r)

    def test_failing_source_control_rejected(self):
        r=deepcopy(self.record)
        r['reviews']['shifted_pressure_source_controls']['checks']['omitted_radial_axial_viscosity_changes_second_pressure_source']=False
        with self.assertRaises(ValueError):
            validate(r)

    def test_missing_gate_or_review_rejected(self):
        for section in ('gates','reviews'):
            r=deepcopy(self.record);r[section].pop(next(iter(r[section])))
            with self.assertRaises(ValueError):
                validate(r)

    def test_missing_remaining_obligation_rejected(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop()
        with self.assertRaises(ValueError):
            validate(r)

    def test_changed_frozen_provenance_rejected(self):
        r=deepcopy(self.record);r['inherited_provenance'][next(iter(r['inherited_provenance']))]=False
        with self.assertRaises(ValueError):
            validate(r)

    def test_changed_current_source_hash_rejected(self):
        r=deepcopy(self.record);r['source_hashes'][next(iter(r['source_hashes']))]='0'*64
        with self.assertRaises(ValueError):
            validate(r)

    def test_nonfinite_diagnostic_rejected(self):
        r=deepcopy(self.record);r['reviews']['shifted_pressure_source_controls']['maximum_source_equation_error']=float('nan')
        with self.assertRaises(ValueError):
            compare_records(self.record,r)

    def test_CLI_does_not_overwrite_frozen_record(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('preserve')
            with self.assertRaises(SystemExit):
                main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__=='__main__':
    unittest.main()
