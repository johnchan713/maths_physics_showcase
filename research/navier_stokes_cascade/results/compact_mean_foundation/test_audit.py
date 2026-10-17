#!/usr/bin/env python3
"""Exact compact-mean controls and scoped evidence mutation checks."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from mean_maps import blocks,solve,determinant,vandermonde,flat_order,matrix_checks,pressure_checks,inverse_checks,normalization_checks
from audit import HERE,FLAGS,generate,validate,compare_records,parent_review,protocol,main


class CompactMeanFoundationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_live_parent_and_unchanged_positive_lambda(self):
        self.assertEqual(protocol()['parent_commit'],'372574a3e9f08bba9d31c351ce04298cd709a0ab')
        self.assertTrue(FLAGS['actual_positive_lambda_not_retuned'])

    def test_all_240_exact_signed_moment_rows(self):
        r=matrix_checks();self.assertEqual(r['moment_rows'],240);self.assertTrue(r['rows'])

    def test_angular_matrix_remains_nonsingular_at_zero_lambda(self):
        self.assertNotEqual(determinant(blocks(2,1)[0]),0)
        self.assertTrue(matrix_checks()['angular_limit_nonsingular'])

    def test_axial_matrix_is_singular_at_zero_lambda(self):
        self.assertEqual(determinant(blocks(2,1)[1]),0)
        with self.assertRaises(ValueError):solve(blocks(2,1)[1],[0,1])

    def test_axial_scaled_inverse_has_exact_one_over_gap_loss(self):
        self.assertTrue(matrix_checks()['sharp_axial_scaled_identity'])
        for q in (F(1,2),F(255,256),F(65535,65536)):
            x=solve(blocks(2,q,3)[1],[0,-7])
            self.assertEqual((1-q)*x[1],F(7,6));self.assertEqual(x[0],-2*x[1])

    def test_axial_inverse_diverges_for_fixed_nonzero_target(self):
        coarse=solve(blocks(2,F(1,2))[1],[0,-1])[1]
        fine=solve(blocks(2,F(255,256))[1],[0,-1])[1]
        self.assertEqual(fine/coarse,128)

    def test_zero_target_remains_zero_near_singular_limit(self):
        self.assertEqual(solve(blocks(2,F(65535,65536))[1],[0,0]),[0,0])

    def test_invalid_geometry_or_moment_rows_rejected(self):
        for args in ((1,F(1,2)),(2,0),(2,F(3,2)),(2,F(1,2),0),(2,F(1,2),1,(1,0,1,1,1))):
            with self.assertRaises(ValueError):blocks(*args)

    def test_vandermonde_orientation_not_absolute_determinant(self):
        nodes=[F(4),F(1,8),F(1,2)]
        self.assertEqual(determinant(vandermonde(nodes)),(nodes[1]-nodes[0])*(nodes[2]-nodes[0])*(nodes[2]-nodes[1]))

    def test_all_576_exact_quadratic_norm_divisor_controls(self):
        r=inverse_checks();self.assertEqual(r['quadratic_integer_norm_controls'],576);self.assertTrue(r['divisors_nonzero'])

    def test_shifted_primitive_all_21_jet_controls(self):
        r=inverse_checks();self.assertEqual(r['shifted_primitive_jet_controls'],21);self.assertTrue(r['shifted_path_identity'])

    def test_unshifted_radial_inverse_is_detectably_wrong(self):
        self.assertTrue(inverse_checks()['omitting_shift_detected'])

    def test_18_strict_margin_choices_keep_one_domain(self):
        r=inverse_checks();self.assertEqual(r['strict_flatness_choices'],18);self.assertTrue(r['strict_flatness_margins'])

    def test_flat_order_requires_strict_not_equal_exponent(self):
        self.assertEqual(flat_order(1,2,0),100001)
        self.assertEqual(flat_order(5,0,0),1)
        with self.assertRaises(ValueError):flat_order(1,2,0,0)

    def test_inverse_requires_zero_radial_Haar_moment(self):
        self.assertIn('zero radial Haar moment required',protocol()['compact_remainder'])

    def test_temporal_inverse_requires_zero_Haar_mode(self):
        self.assertIn('zero Haar input',protocol()['temporal_inverse'])

    def test_lower_frequency_bound_not_only_upper_bound(self):
        self.assertIn('M^(-1)<',protocol()['exact_floor_relation'])

    def test_s_dependent_reserved_swirl_not_frozen_to_q_equals_Q(self):
        self.assertIn('s^(lambda-h)',protocol()['physical_reserved_swirl'])
        self.assertTrue(normalization_checks()['nonconstant_s_swirl_power'])

    def test_all_30_physical_row_normalizations(self):
        r=normalization_checks();self.assertEqual(r['normalization_rows'],30);self.assertTrue(r['all_five_rows'])

    def test_unit_bump_and_compact_pressure_exact(self):
        r=pressure_checks();self.assertTrue(r['unit_bump']);self.assertTrue(r['compact_pressure'])

    def test_complete_radial_pressure_moment(self):
        r=pressure_checks();self.assertTrue(r['pressure_moment']);self.assertTrue(r['pressure_moment_omission_nonzero'])

    def test_moving_pressure_bump_derivative_cannot_be_dropped(self):
        r=pressure_checks();self.assertTrue(r['moving_pressure_moment']);self.assertTrue(r['moving_product_omission_nonzero'])

    def test_complete_covariance_all_six_entries(self):
        r=pressure_checks();self.assertEqual(r['tensor_entries'],6);self.assertTrue(r['axial_covariances_nonzero'])

    def test_initial_full_P_and_Jz_identities(self):
        self.assertTrue(pressure_checks()['complete_initial_defects'])

    def test_radial_axial_flux_not_dropped(self):
        self.assertTrue(pressure_checks()['radial_axial_flux_omission_nonzero'])

    def test_full_background_pressure_and_flat_remainders_retained(self):
        self.assertTrue(FLAGS['actual_background_Pi1_Pi2_and_flat_A0_A1_not_discarded'])

    def test_solved_linear_targets_not_zero_nonlinear_defects(self):
        self.assertTrue(FLAGS['five_map_recomputed_nonlinear_defects_not_identically_zero'])

    def test_no_unconditional_source_residual_or_blowup_promotion(self):
        for k,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):validate(r)

    def test_all_874_inherited_and_15_direct_bindings(self):
        self.assertEqual(len(self.record['inherited_provenance']),874)
        self.assertEqual(len(self.record['source_hashes']),15)

    def test_frozen_evidence_and_JSON_roundtrip(self):
        f=json.loads((HERE/'evidence.json').read_text());validate(f);compare_records(self.record,f)
        validate(json.loads(json.dumps(self.record)))

    def test_provenance_key_substitution_rejected_at_same_count(self):
        r=generate();r['inherited_provenance'].pop(next(iter(r['inherited_provenance'])));r['inherited_provenance']['fake:key']=True
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('fake:key',parent_review()[1])

    def test_cached_checks_are_not_mutable_evidence_aliases(self):
        r=generate();r['reviews']['exact_signed_five_map']['rows']=False
        with self.assertRaises(ValueError):validate(r)
        self.assertTrue(matrix_checks()['rows'])

    def test_protocol_equation_or_remaining_gate_mutation_rejected(self):
        for key in ('protocol','remaining_obligations','gates','source_hashes'):
            r=deepcopy(self.record)
            if isinstance(r[key],dict):r[key].pop(next(iter(r[key])))
            else:r[key].pop()
            with self.assertRaises(ValueError):validate(r)

    def test_cli_will_not_overwrite_any_existing_output(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('keep')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')


if __name__=='__main__':unittest.main()
