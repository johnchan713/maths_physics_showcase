"""Native transport identities and supplied-scale floor regression checks."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review,fixture
from bounds import (directional_envelope,leading_direction,normal_gate,uniform_slow_bound,
                    comparison,certify_actual_inputs)

class NativeFloorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()
    def control(self,key):
        values=[v for k,v in self.r['checks'].items() if key in k]
        self.assertTrue(values,key);self.assertTrue(all(values),key)
    def test_implicit_chart_derivatives(self):self.control('implicit_s_Z_keeps_full_denominator');self.control('eta_Z_has_exact_one_minus_eta_square')
    def test_full_native_Euler_phase_derivative(self):self.control('direct_implicit_chain_matches_native_Euler_direction')
    def test_frozen_q_cutoff_mutation_detected(self):self.control('freezing_q_in_Borel_coefficient_loses_Euler_transport');self.control('actual_q_dependent_cutoff_Euler_term_is_nonzero')
    def test_actual_heat_exterior_cancellation(self):self.control('actual_heat_exterior_direction_cancels_for_full_heat_function')
    def test_heat_cancellation_not_zero_individual_jets(self):self.control('individual_heat_profile_first_jets_are_not_zero');self.control('heat_cancellation_needs_both_native_radial_and_eta_terms')
    def test_sharp_eta_envelope_factorizations(self):self.control('interior_vertex_bound_has_exact_square_remainder');self.control('endpoint_bound_has_exact_positive_factorization')
    def test_sharp_eta_extrema_and_looser_sum(self):self.control('sharp_eta_envelope_attained');self.control('sharp_eta_bound_no_larger_than_sum')
    def test_full_error_independently_dominated_by_retained_powers(self):self.control('independent_full_error_is_bounded_with_distinct_powers')
    def test_independent_primitive_direction_and_diagonal_arithmetic(self):self.control('native_profile_and_full_diagonal_formula_agree')
    def test_retained_powers_can_pass_while_collapsed_program_fails(self):self.control('retained_powers_can_pass_when_collapsed_same_program_fails')
    def test_no_universal_dominance_of_unrelated_profile_norms_claimed(self):self.assertIn('different leading norms',self.record['protocol']['comparison_scope'])
    def test_margin_improves_for_later_S_with_same_inputs(self):self.control('same_bounds_margin_improves_on_later_bands')
    def test_native_norm_domain_matches_original_hull(self):self.control('native_norm_domain_is_exact_original_enlarged_chart')
    def test_original_lower_band_uniform_slow_bound(self):self.control('original_lower_band_implies_whole_slow_gate')
    def test_rounded_carrier_retained_no_band_reselection(self):self.control('rounded_carrier_not_replaced_and_no_new_band')
    def test_missing_band_identification_not_repaired(self):self.control('missing_original_band_gate_is_reported_not_repaired')
    def test_fixed_exponential_series_constant(self):self.control('sixteenth_exponential_series_term_has_positive_remaining_coefficients')
    def test_actual_enclosures_and_blowup_claim_refused(self):
        self.control('no_actual_input_full_residual_or_blowup_claim')
        with self.assertRaisesRegex(NotImplementedError,'Actual leading E0/U0'):certify_actual_inputs()
    def test_supplied_scale_only_and_slow_failure(self):
        r=normal_gate(fixture(),Q(2,5),1,Q(1),Q(1))
        self.assertFalse(r['gates']['original_slow_gate']);self.assertFalse(r['new_scale_or_length_selected'])
    def test_missing_native_jet_rejected(self):
        with self.assertRaises(ValueError):comparison(Q(1,200),dict(value=1),dict(value=1,radial=1,angular=1),Rmin=1,Rmax=2,tilt_upper=1,lambda_max=1,gmin=1,Llow=1,Lhigh=2,Dbox=6)
    def test_bad_growth_and_norms_rejected(self):
        with self.assertRaises(ValueError):directional_envelope(-1,1)
        with self.assertRaises(ValueError):leading_direction(Q(1,100),1,value=1,radial=1,angular=1)
        with self.assertRaises(ValueError):normal_gate(fixture(),2,100,Q(1,1000000),Q(1,1000000))
        with self.assertRaises(ValueError):uniform_slow_bound(Q(0),16)
    def test_all_controls_and_six_records(self):
        self.assertTrue(all(self.r['checks'].values()));self.assertTrue(all(self.record['gates'].values()))
        self.assertEqual(len(self.record['dependency_evidence_sha256']),6)
    def test_round_trip_and_frozen_record(self):
        audit.validate(json.loads(json.dumps(self.record)))
        p=Path(__file__).with_name('evidence.json')
        if p.exists():audit.validate(json.loads(p.read_text()))
    def test_scope_and_scale_mutations_rejected(self):
        for key in ('actual_profile_suprema_diagonal_identity_or_same_band_floor_certified','nearest_nonzero_carrier_replaced_or_new_scale_selected','full_pressure_memory_or_radial_cutoff_curl_dropped'):
            r=deepcopy(self.record);r['flags'][key]=not r['flags'][key]
            with self.assertRaises(ValueError):audit.validate(r)
        r=deepcopy(self.record);r['protocol']['native_identity']='freeze q'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_output_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--output',str(p)])
            audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
