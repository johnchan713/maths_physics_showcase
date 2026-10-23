"""Sharp compact finite-jet operator tests and non-promotion controls."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from bounds import pell,select_mode,finite_jet_obstruction,epsilon_target,compare_request,actual_frontier
from controls import review

class CompactSharpnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()
    def control(self,key):
        values=[v for k,v in self.r['checks'].items() if key in k]
        self.assertTrue(values,key);self.assertTrue(all(values),key)
    def test_actual_quadratic_integer_norm_and_nonzero_divisors(self):self.control('nonzero_integer_norm_and_actual_vr_divisor')
    def test_nonzero_real_Fourier_pair_and_l1_frequency(self):self.control('real_pair_nonzero_k_and_l1_size')
    def test_Pell_selection_at_arbitrary_unmodified_carriers(self):self.control('same_mode_selection_without_retuning_M')
    def test_radial_source_frequency_cost_matches_torus_size(self):self.control('full_radial_phase_frequency_at_most_P')
    def test_exact_full_shift_cancels_the_source_radial_phase(self):self.control('complete_shift_cancels_radial_source_oscillation')
    def test_unshifted_primitive_is_not_the_same_operator(self):self.control('unshifted_primitive_does_not_give_same_integral')
    def test_wrong_drift_orientation_mutation_detected(self):self.control('wrong_drift_sign_is_detected')
    def test_real_source_has_zero_Haar_at_every_radius(self):self.control('real_cosine_pair_has_zero_Haar_and_moment')
    def test_independent_unit_mass_bump_integrals_and_endpoint_jets(self):self.control('unit_mass_and_compact_finite_jet_extension')
    def test_original_f_e_inverse_multiplier_is_retained(self):
        self.control('actual_kernel_f_e_inverse_factor_retained');self.control('dropping_actual_f_e_inverse_changes_kernel_mass')
    def test_independent_integral_detects_missing_shift(self):self.control('omitting_shift_changes_integral_not_just_formula')
    def test_zero_Haar_compensation_does_not_remove_compact_output(self):self.control('zero_Haar_compensation_leaves_nonzero_compact_output')
    def test_raw_product_and_all_lower_order_Leibniz_terms(self):self.control('full_Leibniz_cost_dominates_all_lower_jets')
    def test_finite_jet_fixture_does_not_certify_all_fixed_jets(self):self.control('finite_jet_fixture_not_uniform_all_jet_source')
    def test_two_p_minus_one_source_derivatives_are_insufficient(self):self.control('one_fewer_than_two_p_is_obstructed')
    def test_equality_boundary_is_left_open(self):self.control('equality_not_disproved_or_promoted_to_upper_bound')
    def test_original_upper_rule_only_three_above_M_power_boundary(self):self.control('inherited_upper_rule_is_three_derivatives_above_boundary')
    def test_independent_rational_resonance_ratio_diverges(self):self.control('Pell_subsequence_ratio_grows_exactly_like_P')
    def test_real_noninteger_M_gain_boundary(self):self.control('noninteger_gain_also_has_sharp_power_boundary')
    def test_exact_parent_normalized_and_physical_route_jet_counts(self):self.control('exact_parent_source_jet_and_sharp_route_threshold')
    def test_strict_margin_gives_five_extra_jets_at_integer_target(self):self.control('strict_rule_is_five_above_integer_epsilon_boundary')
    def test_epsilon_deficit_and_unresolved_equality(self):
        self.control('one_jet_below_boundary_has_strict_epsilon_deficit');self.control('equality_is_not_resolved_by_comparison')
    def test_parent_physical_loss_is_not_promoted_to_necessary_loss(self):self.control('same_physical_route_not_new_loss_lower_bound')
    def test_all_three_strict_rule_excess_cases(self):
        self.control('strict_rule_excess3');self.control('strict_rule_excess4');self.control('strict_rule_excess5')
    def test_actual_current_pressure_and_all_covariance_gate_remains_open(self):
        self.control('actual_full_current_pressure_source_inequalities_not_evaluated');self.control('scalar_obstruction_not_canonical_or_blowup_claim')
    def test_large_orders_only_exact_integers_no_giant_factorials(self):
        r=compare_request(0,7,physical_power=11,h=Q(1,1000))
        self.assertGreater(r['retained_raw_source_order'],10**8)
        self.assertIn(r['excess_over_necessary_order'],(3,4,5))
        self.assertEqual(select_mode(Q(10)**100)['norm']**2,1)
    def test_finite_jet_counterexample_distinct_from_actual_input_failure(self):
        f=actual_frontier();self.assertFalse(f['actual_current_source_high_jet_norms_evaluated'])
        self.assertFalse(f['actual_input_validity_residual_decay_or_blowup_certified'])
        self.assertFalse(finite_jet_obstruction(0,1,1)['canonical_source_or_all_jet_counterexample'])
        self.assertFalse(epsilon_target(0,0,0,1)['all_jet_or_actual_source_failure'])
    def test_invalid_parameters_rejected(self):
        for fn,args in ((pell,(0,)),(pell,(-1,)),(select_mode,(Q(1,2),)),
                        (finite_jet_obstruction,(0,0,-1)),(epsilon_target,(0,0,2,1)),
                        (epsilon_target,(0,0,0,1,0)),(compare_request,(2,0))):
            with self.assertRaises(ValueError):fn(*args)
        with self.assertRaises(TypeError):select_mode(1.0)
    def test_four_historical_records_and_fresh_record_reproduce(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),4)
        self.assertTrue(all(self.record['gates'].values()))
        self.assertEqual(self.record['protocol']['new_canonical_analytical_assumptions'],[])
        audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_scope_and_input_provenance_mutations_rejected(self):
        for field in ('flags','dependency_evidence_sha256','remaining_obligations','protocol'):
            bad=deepcopy(self.record)
            if field=='flags':bad[field]['current_source_high_jet_support_pressure_input_inequalities_validated']=True
            elif field=='dependency_evidence_sha256':bad[field][next(iter(bad[field]))]='wrong'
            elif field=='remaining_obligations':bad[field].pop()
            else:bad[field]['all_jet_caveat']='canonical counterexample'
            with self.assertRaises(ValueError):audit.validate(bad)
    def test_historical_output_never_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')

if __name__=='__main__':unittest.main()
