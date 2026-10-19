"""Exact normalized-deficit and independent-freezing checks, not actual acceptance."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import deficit_bounds as b

class DeficitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record=audit.generate();cls.r=cls.record['reviews'];cls.data=b.constants(F(1,2),2,2,5,4,F(1,2),F(3,4))
        cls.rows={k:tuple(map(F,v)) for k,v in cls.r['join']['complete_rows'].items()}
    def symbolic(self,k):self.assertTrue(self.r['symbolic']['checks'][k])
    def controls(self,k):self.assertTrue(all(r['checks'][k] for r in self.r['perturbations']['rows']))
    def test_original_kappa_exact(self):self.assertEqual(b.tilt_slack(F(3,4),F(1,4))['kappa'],F(9,25))
    def test_raw_slack_is_not_relative_squared_margin(self):self.assertEqual(b.tilt_slack(F(3,4),F(1,4))['slack'],F(11,100))
    def test_exact_original_tilt_boundary_uncertified(self):self.assertFalse(b.tilt_slack(F(3,4),F(9,25))['admissible_by_supplied_bound'])
    def test_raw_supremum_above_kappa_uncertified(self):self.assertFalse(b.tilt_slack(F(3,4),F(1,2))['admissible_by_supplied_bound'])
    def test_zero_raw_ratio_maximum_deficit(self):self.assertEqual(b.tilt_slack(F(3,4),0)['slack'],F(9,25))
    def test_coarse_gate_failure_does_not_erase_original_strictness(self):self.assertTrue(self.r['obstruction']['checks']['coarse_gate_failure_does_not_remove_original_admissibility'])
    def test_all_five_noncone_constraints_retained(self):self.assertEqual(len(self.data['base_constraints']),5)
    def test_normalized_Lipschitz_matches_old_exactly(self):
        self.assertTrue(all(r['checks']['normalized_Lipschitz_exactly_matches_legacy'] for r in self.r['perturbations']['legacy_comparisons']))
    def test_cone_radius_linear_in_raw_deficit(self):
        r=b.certificate(self.data,F(1,100));self.assertEqual(r['cone_radius'],self.data['inputs']['omega']**2*F(1,100)/(2*self.data['Lpsi']))
    def test_new_radius_never_smaller_than_legacy(self):self.assertTrue(all(r['checks']['new_radius_never_smaller'] for r in self.r['perturbations']['legacy_comparisons']))
    def test_new_direct_cone_term_strictly_improves_half_gamma_route(self):self.symbolic('direct_deficit_improves_gamma_half_route');self.assertTrue(all(r['checks']['direct_cone_radius_strictly_larger'] for r in self.r['perturbations']['legacy_comparisons']))
    def test_orientation_floor_is_omega_over_three(self):self.assertEqual(b.certificate(self.data,F(1,100))['floors']['oriented_unit_stress'],F(1,6))
    def test_margin_nine_and_frozen_half(self):
        r=b.certificate(self.data,F(1,100));self.assertEqual(r['floors']['unfrozen_cone_margin'],F(1,1296));self.assertEqual(r['floors']['frozen_cone_margin'],F(1,2592));self.symbolic('direct_margin_uses_nine_from_unit_relative_target_bound')
    def test_maximum_deficit_margin_still_below_one(self):self.assertLessEqual(b.certificate(self.data,self.data['kappa'])['floors']['unfrozen_cone_margin'],F(1,9))
    def test_all_ninety_six_exact_controls(self):self.assertEqual(len(self.r['perturbations']['rows']),96);self.assertTrue(self.r['perturbations']['checks']['all_ninety_six_exact_perturbation_controls'])
    def test_complete_normalized_polynomial_identity(self):self.symbolic('exact_normalized_cone_polynomial');self.controls('normalized_cone_polynomial_identity')
    def test_normalized_Lipschitz_controls(self):self.controls('normalized_cone_Lipschitz')
    def test_positive_cone_floor_controls(self):self.controls('normalized_positive_cone_floor')
    def test_spectral_growth_and_mixing_controls(self):self.controls('spectral_and_swirl_floors')
    def test_relative_stress_norm_not_reset(self):self.controls('relative_target_norm_not_reset')
    def test_oriented_unit_target_floor_controls(self):self.controls('oriented_unit_target_floor')
    def test_same_u_cone_margin_controls(self):self.controls('same_unchanged_tilt_new_margin')
    def test_four_thin_margin_controls_keep_fixed_spectral_data(self):self.assertEqual(len(self.r['obstruction']['rows']),4);self.assertTrue(self.r['obstruction']['checks']['all_four_fixed_tilt_thin_margin_controls'])
    def test_strictly_admissible_gaps_and_radii_can_shrink(self):self.assertTrue(self.r['obstruction']['checks']['strictly_decreasing_slack_and_cone_radius'])
    def test_thin_margin_exact_unit_and_cone_factorization(self):self.symbolic('thin_margin_exact_unit_target');self.symbolic('thin_margin_exact_cone_factorization')
    def test_original_gap_identity_includes_actual_supremum(self):self.symbolic('exact_leading_raw_supremum_deficit');self.symbolic('inverse_squared_tilt_coarse_margin_is_raw_slack_over_kappa')
    def test_dyadic_scale_whole_chart_shift(self):self.assertEqual(b.dyadic_scale(F(1,200),102)['integer_floor'],1);self.assertEqual(b.dyadic_scale(F(1,200),2)['integer_floor'],0)
    def test_nonintegral_exponent_rounded_down_conservatively(self):
        r=b.dyadic_scale(F(1,200),103);self.assertEqual(r['exact_exponent'],F(101,100));self.assertEqual(r['factor'],2)
    def test_sixth_root_exact_boundaries_and_strict_above(self):
        for n in (1,2,17,10**50):
            self.assertEqual(b.sixth_root_ceiling(n**6),n);self.assertEqual(b.sixth_root_ceiling(F(n**6)+F(1,10**20)),n+1)
    def test_sixth_root_zero_and_fraction(self):self.assertEqual(b.sixth_root_ceiling(0),0);self.assertEqual(b.sixth_root_ceiling(F(1,1000)),1)
    def test_background_deficit_exact_requirement(self):
        r=b.required_deficit(self.data,F(1,10**20),0,100);self.assertEqual(r['required_deficit'],2*self.data['Lpsi']*F(1,10**20)/self.data['inputs']['omega']**2);self.symbolic('required_background_deficit_exact')
    def test_freezing_deficit_has_seventy_two(self):
        r=b.required_deficit(self.data,0,7,100);self.assertEqual(r['required_deficit'],72*self.data['kappa']*7/(self.data['inputs']['omega']**2*100**6));self.symbolic('freezing_condition_has_seventy_two')
    def test_freezing_not_implied_by_tiny_background_error(self):self.assertTrue(self.r['join']['checks']['freezing_gate_cannot_be_erased_even_when_background_error_tiny'])
    def test_full_pressure_and_viscosity_rows_retained(self):self.assertTrue(self.r['join']['checks']['same_full_pressure_and_viscosity_rows_retained'])
    def test_joint_later_witness_checks_both_costs(self):self.assertTrue(self.r['join']['checks']['whole_chart_and_independent_freezing_witness'])
    def test_sharp_auxiliary_shape_program_not_actual_cutoff(self):self.assertTrue(self.r['join']['checks']['sharp_shape_program_receives_new_radius_but_does_not_select_old_cutoff']);self.assertFalse(audit.FLAGS['auxiliary_sharp_shape_program_identifies_unchanged_cutoff'])
    def test_independent_freezing_can_fail_when_shape_program_passes(self):
        ell=self.r['join']['joint_witness']['sufficient_joint_ell'];c=b.certificate(self.data,F(1,100));C=c['floors']['unfrozen_cone_margin']*ell**6
        r=b.auxiliary_shape_join(self.data,self.rows,F(1,100),F(1,200),ell,C);self.assertFalse(r['freezing_gate']);self.assertTrue(r['auxiliary_shape_program']['possible'])
    def test_fixed_scale_reports_capacity_not_old_norm_verification(self):
        r=b.fixed_scale(self.data,self.rows,1,4,F(1,200),100,0);self.assertFalse(r['shape_budgets_are_possible_for_some_smooth_cutoff']);self.assertEqual(r['shape_capacity'],F(3,4))
    def test_large_required_deficit_cannot_fit_nonnegative_raw_supremum(self):
        r=b.required_deficit(self.data,1,1,2);self.assertFalse(r['deficit_requirement_within_possible_range']);self.assertLess(r['raw_squared_supremum_required_upper'],0)
    def test_later_witness_rejects_impossible_shape_budgets(self):
        with self.assertRaises(ValueError):b.later_witness(self.data,self.rows,1,4,F(1,100),F(1,200),0)
    def test_invalid_tilt_or_raw_supremum(self):
        for u,R in ((0,0),(1,-1),(1,2)):
            with self.assertRaises(ValueError):b.tilt_slack(u,R)
    def test_invalid_leading_bounds(self):
        for args in ((0,2,2,5,4,F(1,2),1),(1,0,2,5,4,F(1,2),1),(1,2,6,5,4,F(1,2),1),(1,2,2,5,4,2,1)):
            with self.assertRaises(ValueError):b.constants(*args)
    def test_invalid_deficit_lower_bound(self):
        for d in (0,-1,self.data['kappa']+1):
            with self.assertRaises(ValueError):b.certificate(self.data,d)
    def test_invalid_scale_h_ell_and_negative_costs(self):
        for h,ell in ((0,100),(F(1,100),100),(F(1,200),1),(F(1,200),True)):
            with self.assertRaises(ValueError):b.dyadic_scale(h,ell)
        with self.assertRaises(ValueError):b.required_deficit(self.data,-1,1,100)
        with self.assertRaises(ValueError):b.sixth_root_ceiling(-1)
    def test_exact_inputs_reject_floats_and_bools(self):
        for value in (True,.1):
            with self.assertRaises(TypeError):b.tilt_slack(value,0)
            with self.assertRaises(TypeError):b.sixth_root_ceiling(value)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'abfaf6c133a1b24a88c3eb55f54ac0546178a770');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()))
    def test_frozen_record_reproduces(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_sources_preserved(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_all_unsupported_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_scope_type_mutation_rejected(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_freezing_obligation_cannot_be_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(1)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_false_claim_original_choice_passes_coarse_gate_rejected(self):
        r=deepcopy(self.record);r['flags']['original_admissibility_implies_the_stronger_coarse_numeric_tilt_test']=True
        with self.assertRaises(ValueError):audit.validate(r)
    def test_derived_margin_tampering_rejected(self):
        r=deepcopy(self.record);r['reviews']['join']['auxiliary_shape_join']['certificate']['floors']['unfrozen_cone_margin']='1'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_cli_fresh_replay_and_refusal_to_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
