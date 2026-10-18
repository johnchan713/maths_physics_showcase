"""Exact rational regression and conditional-scope tests; no actual profile evaluation."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import gap_bounds as b

class BackgroundGapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record=audit.generate();cls.c=b.certificate(F(1,2),2,2,5,4,F(1,2),F(1,2),1)
        cls.rows=cls.record['reviews']['explicit_background_perturbations']['rows']
    def controls(self,key):self.assertTrue(all(r['checks'][key] for r in self.rows),key)
    def test_exact_fixture_radius(self):self.assertEqual(self.c['radius'],F(1,68992))
    def test_all_six_safety_constraints_charged(self):
        self.assertEqual(len(self.c['constraints']),6)
        for x in self.c['constraints'].values():self.assertLessEqual(self.c['radius'],x)
    def test_positive_explicit_floors(self):self.assertTrue(all(v>0 for v in self.c['floors'].values()))
    def test_same_tilt_margin_strict_between_zero_and_one(self):self.assertTrue(0<self.c['floors']['cone_margin']<1)
    def test_error_controls_fit_euclidean_norms(self):self.controls('input_bounds')
    def test_leading_data_satisfy_all_original_inputs(self):self.controls('leading_input_bounds')
    def test_exact_full_discriminant_expansion(self):self.controls('exact_full_discriminant_expansion')
    def test_discriminant_expansion_budget(self):self.controls('discriminant_difference')
    def test_product_difference_budget(self):self.controls('product_difference')
    def test_exact_quotient_difference_identity(self):self.controls('exact_quotient_identity')
    def test_quotient_lipschitz_budget(self):self.controls('quotient_difference')
    def test_orientation_difference_budget(self):self.controls('H_difference')
    def test_transverse_direction_difference_budget(self):self.controls('J_difference')
    def test_cone_polynomial_difference_budget(self):self.controls('cone_polynomial_difference')
    def test_actual_swirl_floor(self):self.controls('swirl_floor')
    def test_actual_shear_floor(self):self.controls('angular_shear_floor')
    def test_actual_discriminant_floor(self):self.controls('discriminant_floor')
    def test_actual_spectral_floor(self):self.controls('lambda_floor')
    def test_actual_mixing_floor(self):self.controls('c_floor')
    def test_actual_normal_theta_floor(self):self.controls('normal_floor')
    def test_actual_orientation_floor(self):self.controls('orientation')
    def test_actual_cone_polynomial_floor(self):self.controls('cone_floor')
    def test_same_fixed_tilt_has_explicit_margin(self):self.controls('same_tilt_margin')
    def test_all_thirty_two_signed_perturbation_controls(self):self.assertEqual(len(self.rows),32)
    def test_twenty_seven_independent_reference_frames(self):
        r=b.frame_review();self.assertEqual(len(r['rows']),27);self.assertTrue(all(r['checks'].values()))
    def test_positive_swirl_and_shear_can_lose_growth(self):self.assertTrue(b.obstruction_review()['checks']['positive_swirl_and_shear_not_sufficient'])
    def test_tiny_error_not_a_uniform_near_neutral_certificate(self):
        r=b.obstruction_review()['near_neutral_growth_rows'];self.assertLess(F(r[-1]['absolute_error']),F(1,10**20));self.assertTrue(all(x['sign_failure'] for x in r))
    def test_strict_growth_not_sufficient_for_prescribed_cone(self):self.assertTrue(b.obstruction_review()['checks']['strict_growth_does_not_supply_same_tilt_cone'])
    def test_absolute_stress_error_is_not_relative_input(self):self.assertTrue(b.obstruction_review()['checks']['unweighted_absolute_stress_error_not_relative_control'])
    def test_allowed_diagnostic_h_selector_and_all_later_grid_points(self):self.assertTrue(all(b.tail_review()['checks'].values()))
    def test_whole_four_Q_chart_pays_two_extra_bands(self):
        r=b.tail_review();self.assertEqual(r['sufficient_first_fixture_grid_band'],r['first_fixture_grid_point']+2)
    def test_small_tilt_keeps_positive_radius_without_retuning(self):
        c=b.certificate(F(1,2),2,2,5,4,F(1,2),F(1,2),F(1,10**6));self.assertGreater(c['radius'],0);self.assertEqual(c['inputs']['u'],F(1,10**6))
    def test_arbitrarily_small_positive_floors_not_rounded_to_zero(self):
        c=b.certificate(F(1,10**100),2,F(1,10**100),5,F(1,10**200),F(1,2),F(1,2),1);self.assertGreater(c['radius'],0)
    def test_vanishing_cone_margin_is_rejected(self):
        with self.assertRaises(ValueError):b.certificate(1,2,1,5,1,F(1,2),0,1)
    def test_zero_tilt_is_rejected(self):
        with self.assertRaises(ValueError):b.certificate(1,2,1,5,1,F(1,2),F(1,2),0)
    def test_unordered_or_invalid_inputs_rejected(self):
        for args in ((2,1,1,5,1,1,1,1),(1,2,6,5,1,1,1,1),(0,2,1,5,1,1,1,1),(1,2,1,5,1,2,1,1)):
            with self.assertRaises(ValueError):b.certificate(*args)
    def test_float_and_bool_cannot_masquerade_as_exact_parameters(self):
        for x in (1.,True):
            with self.assertRaises(TypeError):b.certificate(x,2,1,5,1,F(1,2),F(1,2),1)
    def test_inconsistent_exact_norm_rejected(self):
        with self.assertRaises(ValueError):b.shear(1,-4,0,3)
    def test_neutral_discriminant_has_zero_spectral_and_mixing_gap(self):
        r=b.shear(1,-2,0,2);self.assertEqual(r['lambda2'],0);self.assertEqual(r['c2'],0)
    def test_protocol_and_parent(self):
        self.assertEqual(audit.protocol()['parent_commit'],'165b75afee625f8ec597e117e4986452907e7a3c');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()))
    def test_frozen_exact_record_reproduces(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_recursive_historical_bytes_unchanged(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_all_unsupported_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_derived_radius_tamper_rejected(self):
        r=deepcopy(self.record);r['reviews']['explicit_background_perturbations']['certificate']['radius']='1'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_remaining_actual_input_gate_cannot_be_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_bool_integer_scope_tamper_rejected(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_cli_record_and_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
