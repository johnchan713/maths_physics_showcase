"""Sharp capacity/infimum/elimination checks; auxiliary shapes are not the old field."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import json,tempfile,unittest
from pathlib import Path
import audit
import shape_bounds as b

class ShapeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record=audit.generate();cls.r=cls.record['reviews'];cls.rows={k:tuple(map(F,v)) for k,v in cls.r['gate']['toy_rows'].items()}
    def capacity_check(self,k):self.assertTrue(self.r['capacity']['checks'][k])
    def gate_check(self,k):self.assertTrue(self.r['gate']['checks'][k])
    def symbolic(self,k):self.assertTrue(self.r['symbolic']['checks'][k])
    def test_trapezoid_exact_integral(self):self.symbolic('exact_clipped_trapezoid_integral');self.assertEqual(b.capacity(F(3,2),8),F(39,32))
    def test_triangle_exact_integral(self):self.symbolic('exact_triangular_integral');self.assertEqual(b.capacity(8,12),3)
    def test_branch_boundary_agreement(self):self.symbolic('branches_agree_at_M_equals_N_over_two');self.assertEqual(b.capacity(4,8),2)
    def test_all_twenty_five_exact_integrations(self):self.assertEqual(len(self.r['capacity']['exact_rows']),25);self.capacity_check('all_twenty_five_exact_capacity_controls')
    def test_old_corner_capacity_three_quarters(self):self.capacity_check('rectangular_lower_corner_unrealizable');self.symbolic('rectangular_lower_corner_capacity_is_three_quarters')
    def test_M_one_never_enough(self):
        for N in (1,4,8,10**100):self.assertFalse(b.shape_budget(1,N)['smooth_shape_exists'])
    def test_N_four_never_enough(self):
        for M in (F(1,2),1,2,3,10**100):self.assertFalse(b.shape_budget(M,4)['smooth_shape_exists'])
    def test_exact_triangular_frontier_equality_rejected(self):self.assertEqual(b.capacity(2,4),1);self.assertFalse(b.shape_budget(2,4)['smooth_shape_exists'])
    def test_exact_trapezoid_frontier_equality_rejected(self):self.assertFalse(b.shape_budget(F(3,2),F(9,2))['smooth_shape_exists']);self.symbolic('strict_frontier_from_area_one')
    def test_all_four_strict_frontier_padding_controls(self):self.assertEqual(len(self.r['capacity']['near_frontier']),4);self.capacity_check('all_four_strict_frontier_and_smoothing_controls')
    def test_huge_loose_upper_enclosures_not_exact_norms(self):self.capacity_check('loose_first_upper_bound_does_not_force_second_above_its_square');self.assertTrue(b.shape_budget(1000,1216)['smooth_shape_exists'])
    def test_frontier_t_parameter_identity(self):self.symbolic('t_parameter_frontier')
    def test_arbitrarily_close_to_M_one_requires_large_N(self):
        e=F(1,10**100);M=1+e;N=M*M/e
        self.assertFalse(b.shape_budget(M,N)['smooth_shape_exists']);self.assertTrue(b.shape_budget(M,N+1)['smooth_shape_exists'])
    def test_smoothing_padding_exact_formula(self):
        M,N=F(3,2),F(5);A=b.capacity(M,N);w=b.smoothing_witness(M,N)
        self.assertEqual(w['epsilon'],(A-1)/(16*A));self.assertEqual(w['pre_normalization_mass'],(1-2*w['epsilon'])**2*A)
    def test_smoothing_support_strictly_inside_same_interval(self):
        w=b.smoothing_witness(32,1216);lo,hi=w['slope_support'];self.assertGreater(lo,0);self.assertLess(hi,1);self.assertEqual(lo,1-hi)
    def test_smoothing_both_derivative_bounds_strict(self):
        w=b.smoothing_witness(32,1216);self.assertLess(w['first_norm_upper'],32);self.assertLess(w['second_norm_upper'],1216);self.assertTrue(all(w['checks'].values()))
    def test_smoothing_normalized_enclosures_jointly_feasible(self):
        w=b.smoothing_witness(F(3,2),F(9,2)+F(1,2**60));self.assertGreater(b.capacity(w['first_norm_upper'],w['second_norm_upper']),1)
    def test_smoothing_refuses_nonstrict_capacity(self):
        for M,N in ((1,4),(2,4),(F(3,2),F(9,2))):
            with self.assertRaises(ValueError):b.smoothing_witness(M,N)
    def test_exact_weighted_stationary_identity(self):self.symbolic('weighted_infimum_at_stationary_t')
    def test_weighted_high_M_comparison(self):self.symbolic('higher_M_comparison_square')
    def test_all_fifteen_weighted_approximation_controls(self):self.assertEqual(len(self.r['weighted']['rows']),15);self.assertTrue(self.r['weighted']['checks']['all_fifteen_exact_weighted_infimum_controls'])
    def test_unattained_stress_floor_boundary_rejected(self):self.assertFalse(b.stress_floor_test(0,3,1,9)['passes_necessary_strict_stress_floor'])
    def test_positive_slack_above_weighted_floor_accepted(self):self.assertTrue(b.stress_floor_test(0,3,1,F(9)+F(1,10**100))['passes_necessary_strict_stress_floor'])
    def test_zero_first_weight_still_strict_second_floor(self):self.assertFalse(b.stress_floor_test(0,0,1,4)['passes_necessary_strict_stress_floor']);self.assertTrue(b.stress_floor_test(0,0,1,F(401,100))['passes_necessary_strict_stress_floor'])
    def test_negative_shift_not_rescued_by_squaring(self):self.assertFalse(b.stress_floor_test(100,3,1,1)['passes_necessary_strict_stress_floor'])
    def test_old_relaxed_corner_pass_is_not_smooth_feasibility(self):self.gate_check('old_corner_pass_can_coexist_with_no_smooth_shape_certificate')
    def test_exact_weighted_infimum_no_shape(self):self.gate_check('exact_unattainable_stress_infimum_rejected')
    def test_strictly_above_infimum_can_have_shape(self):self.gate_check('strict_above_infimum_can_work')
    def test_exact_quadratic_square_completion(self):self.symbolic('square_complete_for_quadratic_minimum')
    def test_exact_capped_minimum_identity(self):self.symbolic('capped_minimum_difference')
    def test_derivative_cap_equality_cannot_use_nonsmooth_limit(self):self.gate_check('derivative_cap_equality_has_no_smooth_shape')
    def test_derivative_cap_positive_slack(self):self.gate_check('derivative_cap_strict_surplus_has_smooth_shape')
    def test_rational_witness_has_all_original_rows(self):
        r=b.shape_gate(self.rows,10);self.assertTrue(r['possible']);self.assertTrue(all(r['checks'].values()));self.assertEqual(r['M'],F(13,8));self.assertEqual(r['N'],F(187,40))
    def test_constant_row_failure_not_hidden_by_quadratic(self):
        rows=dict(self.rows,value=(100,0,0));self.assertEqual(b.shape_gate(rows,10)['reason'],'constant_row_exceeds_budget')
    def test_first_norm_cap_at_one_excludes_shape(self):
        rows=dict(self.rows,derivative=(0,10,0));self.assertEqual(b.shape_gate(rows,10)['reason'],'first_shape_budget_at_most_one')
    def test_nonpositive_stress_surplus_excludes_shape(self):
        rows=dict(self.rows,stress=(100,3,1));self.assertEqual(b.shape_gate(rows,10)['reason'],'nonpositive_stress_surplus')
    def test_no_square_roots_or_float_parameters_in_witness(self):
        r=b.shape_gate(self.rows,10)
        for k in ('M','N','quadratic_minimizer','quadratic_minimum','N_frontier','N_stress_cap'):self.assertIsInstance(r[k],F)
    def test_complete_corrected_width_full_memory_composition(self):self.gate_check('complete_corrected_width_pressure_and_viscosity_rows_composed');self.gate_check('no_source_term_or_old_field_replaced_by_elimination')
    def test_160_gate_variants_and_4000_auxiliary_candidate_controls(self):
        old=b.gate_module();cases=candidates=0
        for A0,A1,A2,D0,D1,B in product((0,5),(0,3),(1,2),(0,7),(1,3),(4,8,10,20,100)):
            rows=dict(self.rows,derivative=(D0,D1,0),stress=(A0,A1,A2));r=b.shape_gate(rows,B);cases+=1
            if r['possible']:self.assertTrue(all(r['checks'].values()));self.assertTrue(b.shape_budget(r['M'],r['N'])['smooth_shape_exists'])
            for M,N in product((F(5,4),F(3,2),F(7,4),F(2),F(3)),(F(41,10),F(9,2),F(5),F(8),F(32))):
                candidates+=1
                if b.capacity(M,N)>1 and old.shape_region(rows,B,M,N)['accepted']:self.assertTrue(r['possible'])
        self.assertEqual((cases,candidates),(160,4000))
    def test_invalid_capacity_domain(self):
        for M,N in ((0,1),(1,0),(-1,4)):
            with self.assertRaises(ValueError):b.capacity(M,N)
    def test_nonexact_capacity_inputs_rejected(self):
        for value in (True,.1):
            with self.assertRaises(TypeError):b.capacity(value,4)
    def test_invalid_affine_row_schema(self):
        for rows in ({},dict(self.rows,stress=(1,2)),dict(self.rows,stress=(-1,0,1)),dict(self.rows,value=(1,1,0)),dict(self.rows,derivative=(0,0,0)),dict(self.rows,stress=(0,1,0))):
            with self.assertRaises(ValueError):b.shape_gate(rows,10)
    def test_invalid_budget_and_weights(self):
        with self.assertRaises(ValueError):b.shape_gate(self.rows,0)
        with self.assertRaises(ValueError):b.stress_floor_test(0,-1,1,10)
        with self.assertRaises(ValueError):b.stress_floor_test(0,1,0,10)
    def test_nonexact_gate_and_weight_inputs_rejected(self):
        with self.assertRaises(TypeError):b.shape_gate(self.rows,10.)
        with self.assertRaises(TypeError):b.stress_floor_test(0,1,True,10)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'355dce1e76bf0e6d97ec7da84aa06e0fddff5522');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()))
    def test_frozen_record_exactly_reproduced(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_source_bytes_unchanged(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_all_unsupported_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_scope_type_mutation_rejected(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_strict_frontier_cannot_be_changed_to_closed(self):
        r=deepcopy(self.record);r['reviews']['gate']['capped_boundary']['possible']=True
        with self.assertRaises(ValueError):audit.validate(r)
    def test_no_promotion_of_auxiliary_shape_to_actual_cutoff(self):
        self.assertFalse(audit.FLAGS['auxiliary_feasible_shape_or_budget_is_the_unchanged_actual_cutoff']);self.assertFalse(audit.FLAGS['actual_cutoff_shape_existing_tilt_or_retained_scale_acceptance_verified'])
    def test_actual_shape_obligation_cannot_be_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_cli_fresh_replay_and_refusal_to_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
