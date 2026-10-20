"""Tame fourth jets, all normalized rows, primitive gate and provenance controls."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import target_bounds as b

class FourthTargetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_exponential_fourth_all_partitions(self):self.assertTrue(self.r['symbolic']['checks']['exponential_fourth_has_all_partitions'])
    def test_swirl_square_exponential_is_different(self):self.assertTrue(self.r['symbolic']['checks']['squared_swirl_exponential_fourth_has_all_partitions'])
    def test_mixed_lower_jets_not_removed(self):self.assertTrue(self.r['symbolic']['checks']['omitting_mixed_lower_jets_detected'])
    def test_tame_first_exponential_lower_bound(self):self.assertTrue(self.r['checks']['first_exponential_lower_quartic_terms_fit_three'])
    def test_tame_squared_exponential_lower_bound(self):self.assertTrue(self.r['checks']['squared_exponential_lower_quartic_terms_fit_96'])
    def test_tame_product_b_exponential_lower_bound(self):self.assertTrue(self.r['checks']['mixed_b_exponential_fourth_lower_terms_fit23'])
    def test_all_five_density_affine_rows(self):
        d=b.density_affines(3);self.assertEqual(d,dict(U=(1,0,1),E=(36,9,0),U2=(18,0,12),E2=(1026,162,0),UE=(342,54,27)))
    def test_density_highest_derivatives_only_linear(self):
        d=b.density_affines(3)
        for row in d.values():self.assertEqual(b.value(row,64,4,6)-b.value(row,64,2,3),b.value(row,0,2,3))
    def test_native_fourth_envelope_not_dropped(self):self.assertGreater(b.density_affines(5)['E2'][1],b.density_affines(3)['E2'][1])
    def test_raw_M_J_I_moment_geometry(self):
        m=b.moment_affines(3,3);self.assertEqual(m['M'],(3,0,3));self.assertEqual(m['J'],(6156,972,486));self.assertEqual(m['I'],(648,162,0))
    def test_raw_S_retains_both_squares(self):self.assertEqual(b.moment_affines(3,3)['S'],(1593,243,36))
    def test_raw_pressure_has_inverse_X_and_half(self):self.assertEqual(b.moment_affines(3,3)['Cp'],(4617,729,0))
    def test_all_five_normalizers_distinct(self):
        q=b.target_affines(3,3,2,F(1,4),3,F(1,200));self.assertEqual(q['normalizer_C4'],dict(M=F(5,6),J=F(25,36),I=F(5,12),S=F(25,18),Cp=F(25,9)))
    def test_reciprocal_e_has_fourth_norm5(self):self.assertTrue(self.r['symbolic']['checks']['native_reciprocal_e_normalizer_C4_is_five'])
    def test_reciprocal_e_squared_has_fourth_norm25(self):self.assertTrue(self.r['symbolic']['checks']['native_reciprocal_e_squared_normalizer_C4_is25'])
    def test_inverse_lambda_transformed_row(self):
        q=b.target_affines(3,3,2,F(1,4),3,F(1,200));self.assertEqual(q['transformed_rows'][1],(855500,135000,68000))
    def test_no_zero_lambda_limit_is_used(self):
        q=b.target_affines(3,3,2,F(1,4),3,F(1,200));qq=b.target_affines(3,3,2,F(1,4),3,F(1,400));self.assertEqual(qq['transformed_rows'][1],tuple(2*v for v in q['transformed_rows'][1]))
    def test_componentwise_max_encloses_all_rows(self):self.assertTrue(self.r['checks']['componentwise_max_bounds_every_row'])
    def test_explicit_Z4_affine_constants(self):self.assertEqual(b.target_affines(3,3,2,F(1,4),3,F(1,200))['Z4_affine'],(855500,135000,68000))
    def test_Z4_no_longer_independent_unknown(self):self.assertTrue(self.r['checks']['target_Z4_not_an_additional_independent_unknown'])
    def test_B4_full_norm_includes_lower_jets(self):self.assertTrue(self.r['checks']['axial_primitive_B4_retains_lower_jets'])
    def test_same_pressure_axis_datum_only_cancels_in_difference(self):self.assertTrue(self.r['symbolic']['checks']['same_pressure_axis_datum_cancels_only_in_difference'])
    def test_changed_pressure_gauge_has_nonzero_fourth_jet(self):self.assertTrue(self.r['symbolic']['checks']['changed_axis_gauge_would_add_a_fourth_jet'])
    def test_moving_inverse_phase_cannot_be_frozen(self):self.assertEqual(self.r['symbolic']['moving_phase_fourth_control'],'9/512')
    def test_fifteen_full_target_samples(self):
        d=self.r['manufactured'];self.assertEqual(len(d['rows']),15);self.assertTrue(d['checks']['all_fifteen_full_target_controls']);self.assertTrue(d['manufactured_not_selected_loop_or_PDE']);self.assertTrue(d['floating_not_interval_certificate'])
    def test_original_row_preconditioner_independent_controls(self):self.assertTrue(all(r['checks']['original_preconditioner_fourth_jet_matches'] for r in self.r['manufactured']['rows']))
    def test_full_normalizer_derivative_Leibniz_controls(self):self.assertTrue(all(r['checks']['independent_normalizer_Leibniz_convolution'] for r in self.r['manufactured']['rows']))
    def test_manufactured_native_and_lower_inputs_genuinely_bounded(self):self.assertTrue(self.r['manufactured']['checks']['actual_manufactured_lower_jet_inputs_satisfy_H']);self.assertTrue(self.r['manufactured']['checks']['actual_manufactured_native_C4_inputs_satisfy_Ahi'])
    def test_manufactured_eta_endpoints_and_phases(self):
        rows=self.r['manufactured']['rows'];self.assertEqual({r['eta'] for r in rows},{-1.,-.4,0.,.6,1.});self.assertEqual(len({r['phase'] for r in rows}),3)
    def test_exact_substituted_gate_boundary(self):self.assertTrue(self.r['checks']['primitive_gate_is_exact_at_substituted_boundary'])
    def test_gate_matches_parent_substituted_majorant(self):
        for a4 in (0,2,1000):
            for b4 in (0,3,1000):
                for K in (0,1,F(10**20)):
                    q=b.primitive_gate(3,3,2,F(1,4),3,F(1,200),64,4096,F(1,10**9),a4,b4,K);self.assertEqual(q['accepts'],q['substituted_majorant_accepts'])
    def test_highest_jet_need_not_be_small_for_density_bounds(self):
        q=b.primitive_gate(3,3,2,F(1,4),3,F(1,200),64,4096,F(1,10**9),10**30,10**40,10**100);self.assertTrue(q['accepts'])
    def test_negative_headroom_rejects_this_majorant(self):
        q=b.primitive_gate(3,3,2,F(1,4),3,F(1,200),64,4096,0,0,0,0);self.assertLess(q['rhs'],0);self.assertFalse(q['accepts'])
    def test_highest_primitive_gate_coefficients(self):
        q=b.primitive_gate(3,3,2,F(1,4),3,F(1,200),64,4096,0,2,3,0);self.assertEqual(q['a4_coefficient'],3200000*3*135000);self.assertEqual(q['b4_coefficient'],1+3200000*3*68000)
    def test_rational_sqrt_enclosure_rejects_underbound(self):
        with self.assertRaises(ValueError):b.target_affines(3,3,2,F(1,5),3,F(1,200))
    def test_sqrt_enclosure_monotonicity(self):
        a=b.target_affines(3,3,2,F(1,4),3,F(1,200));c=b.target_affines(3,3,2,F(1,2),3,F(1,200));self.assertTrue(all(x<=y for x,y in zip(a['Z4_affine'],c['Z4_affine'])))
    def test_invalid_normalizer_inputs(self):
        for rp,sr,es,lam in ((0,1,1,F(1,200)),(2,0,3,F(1,200)),(2,F(1,4),0,F(1,200)),(2,F(1,4),3,0),(2,F(1,4),3,F(1,50))):
            with self.assertRaises(ValueError):b.target_affines(3,3,rp,sr,es,lam)
    def test_invalid_geometry_native_envelope(self):
        with self.assertRaises(ValueError):b.density_affines(0)
        with self.assertRaises(ValueError):b.moment_affines(0,3)
    def test_invalid_tame_lower_jet_hypothesis(self):
        for H,N in ((0,1),(2,1)):
            with self.assertRaises(ValueError):b.primitive_gate(3,3,2,F(1,4),3,F(1,200),H,N,0,0,0,1)
    def test_float_bool_inputs_rejected(self):
        for v in (3.,True):
            with self.assertRaises(TypeError):b.density_affines(v)
    def test_parent_protocol_and_frozen_record(self):self.assertEqual(audit.protocol()['parent_commit'],'4f093b36c10abea2fb4c8e6b31bb522418df7963');audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_source_bindings(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_unproved_scope_promotion_rejected(self):
        for key,value in audit.FLAGS.items():
            if not value:
                r=deepcopy(self.record);r['flags'][key]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_numeric_zero_is_not_false(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_normalizers_gauge_pressure_or_moving_phase_not_erased(self):
        for key in ('normalizers','transform','raw','inputs','pressure'):
            r=deepcopy(self.record);r['protocol'][key]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_remaining_actual_fourth_primitive_obligation_not_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_evidence_and_frozen_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
