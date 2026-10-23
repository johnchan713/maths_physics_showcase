"""Leading directional enclosure, domain, hierarchy and selector checks."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review,fixture
from bounds import (collar_recipe,compact_majorants,leading_enclosures,
    reference_enclosures,native_first_jets,comparison,selector_margin,certify_same_band)

class LeadingDirectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()
    def control(self,key):
        vals=[v for k,v in self.r['checks'].items() if key in k]
        self.assertTrue(vals,key);self.assertTrue(all(vals),key)
    def test_first_phase_pullback(self):self.control('pullback_first_radial_has_P_phi_without_N')
    def test_angular_phase_pullback(self):self.control('pullback_angular_retains_inverse_N')
    def test_complete_swirl_product(self):self.control('exponential_swirl_first_product_rule')
    def test_axial_additive_pullback(self):self.control('axial_first_derivative_has_same_pullback')
    def test_second_radial_frequency_loss_not_erased(self):self.control('second_radial_still_has_N')
    def test_dropped_fast_phase_mutation(self):self.control('omitting_phase_chain_term_detected')
    def test_regular_swirl_half_power(self):self.control('regular_swirl_conversion_retains_half_power');self.control('using_A_instead_of_A_plus_half_detected')
    def test_pointwise_log_jet_dual_bound_sharp(self):self.control('pointwise_dual_norm_two_is_sharp')
    def test_radius_majorant_has_no_square_root(self):self.control('sqrt_radius_bound_has_positive_remainder')
    def test_joint_norm_counts_pullback_terms(self):self.control('joint_J1_budget')
    def test_compact_hierarchy_and_same_repair(self):self.control('first_E_below_H4');self.control('first_U_below_H2');self.control('leading_repair_below_one')
    def test_direction_uses_regional_pointwise_sums(self):self.control('axis_direction_below_pointwise_dual_bound');self.control('regional_max_of_sums_not_global_sum_of_suprema')
    def test_exact_axis_eta_envelope(self):self.control('axis_sharp_envelope_dominates')
    def test_whole_hull_ordinary_X_inverse_radius_and_all_keys(self):
        self.control('whole_hull_ordinary_X_keeps_inverse_radius');self.control('whole_hull_has_all_two_profile_primitive_keys');self.control('whole_hull_has_original_Xlo')
        with self.assertRaises(ValueError):native_first_jets(1,1,1,16,0)
    def test_axis_recipe_retains_width_and_exponential(self):self.control('narrower_width_increases_same_collar_prefactor');self.control('axis_recipe_no_giant_exponential_evaluation')
    def test_reference_growth_shear_and_tilt_outward(self):
        self.control('reference_gmin_outward');self.control('reference_growth_floor_outward');self.control('reference_growth_ceiling_outward');self.control('original_squared_tilt_not_changed')
    def test_existing_selector_margin_pass_and_failure(self):self.control('margin_has_expected_outcome');self.control('same_selector_polynomial_independent_formula')
    def test_true_radius_and_slow_gate(self):self.control('original_radius_factor_retained');self.control('original_selector_supplies_slow')
    def test_same_comparison_all_later_S(self):self.control('later_S_improves_same_comparison')
    def test_no_selector_reselection_or_actual_margin_claim(self):self.control('no_scale_or_Lreq_reselection');self.control('no_actual_margin_claim')
    def test_first_coefficient_retained_in_complete_tail(self):self.control('full_first_order_tail_includes_first_coefficient');self.control('first_term_cannot_be_removed')
    def test_wrong_hierarchy_and_collar_rejected(self):
        for args in ((16,16,16**32+1),(1,16,16**32),(1,8,8**32+1)):
            with self.assertRaises(ValueError):compact_majorants(*args)
        with self.assertRaises(ValueError):collar_recipe(2,1,8,Q(1,100),Q(1,100))
        with self.assertRaises(ValueError):leading_enclosures(Q(1,100),1,1,1,16)
    def test_invalid_reference_selector_and_radius_rejected(self):
        with self.assertRaises(ValueError):reference_enclosures(0,1,1,3,1,1,1)
        with self.assertRaises(ValueError):selector_margin(fixture(),1,25,Q(1,100),1,1)
        with self.assertRaises(ValueError):selector_margin(fixture(),1,25,Q(1,6057072),Q(1,2),1)
    def test_same_band_claim_refused(self):
        with self.assertRaisesRegex(NotImplementedError,'Actual same-Lreq'):certify_same_band()
        self.control('no_full_residual_or_blowup_promotion')
    def test_all_controls_ten_records_and_frozen_roundtrip(self):
        self.assertTrue(all(self.r['checks'].values()));self.assertTrue(all(self.record['gates'].values()))
        self.assertEqual(len(self.record['dependency_evidence_sha256']),10)
        audit.validate(json.loads(json.dumps(self.record)))
        p=Path(__file__).with_name('evidence.json')
        if p.exists():audit.validate(json.loads(p.read_text()))
    def test_claim_scope_and_protocol_mutations_rejected(self):
        for key in ('actual_giant_parameters_or_same_Lreq_margin_numerically_verified','new_selector_Lreq_tilt_width_frequency_or_support_selected','unconditional_pulse_complete_residual_infinite_field_or_blowup_proved'):
            r=deepcopy(self.record);r['flags'][key]=True
            with self.assertRaises(ValueError):audit.validate(r)
        r=deepcopy(self.record);r['protocol']['selector']='increase Lreq'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_output_and_frozen_record_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--output',str(p)])
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
