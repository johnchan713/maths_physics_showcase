"""Sharp affine-program algebra and coupled scale controls; no actual cutoff acceptance."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import sympy as s
import audit
import cost_floor as b

class SharpCurrentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record=audit.generate();cls.r=cls.record['reviews']
        cls.data=b.current().deficit().constants(F(1,2),2,2,6,4,F(1,2),F(3,4))
        cls.cost=b.current().current_costs(cls.data,0,0,6,F(1,20000));cls.slack=F(1,100)
    def symbolic(self,k):self.assertTrue(self.r['symbolic']['checks'][k])
    def test_exact_quadratic_sign_guards_opposite_sign_squaring(self):
        self.assertLess(b.qsign(b.quad(-3,2,2)),0);self.assertGreater(b.qsign(b.quad(3,-2,2)),0)
        self.assertGreater(b.qsign(b.quad(-1,1,2)),0);self.assertLess(b.qsign(b.quad(1,-1,2)),0)
    def test_quadratic_sign_agrees_with_independent_symbolic_comparisons(self):
        for a in (-3,-1,0,1,3):
            for c in (-3,-1,0,1,3):
                for q in (F(1,2),F(2),F(7)):
                    v=b.quad(a,c,q);e=b.as_symbolic(v);self.assertEqual(b.qsign(v),int(bool(e>0))-int(bool(e<0)))
    def test_perfect_rational_roots_and_zero_radicand(self):self.assertEqual(b.quad(1,2,F(9,4)),b.quad(4));self.assertEqual(b.quad(7,123,0),b.quad(7))
    def test_quadratic_comparison_and_affine_map(self):
        x=b.quad(1,2,2);self.assertGreater(b.qcmp(x,b.quad(3)),0);self.assertEqual(b.linear(x,3,2),b.quad(5,4,2))
    def test_different_irrational_fields_are_not_silently_compared(self):
        with self.assertRaises(ValueError):b.qcmp(b.quad(0,1,2),b.quad(0,1,3))
    def test_negative_radicand_rejected(self):
        with self.assertRaises(ValueError):b.quad(0,1,-1)
    def test_all_ten_closed_optimizer_regimes(self):self.assertEqual(len(self.r['optimizer']['regimes']),10);self.assertTrue(self.r['optimizer']['checks']['all_ten_exact_optimizer_regimes'])
    def test_all_414_independent_budget_decisions(self):self.assertEqual(len(self.r['optimizer']['budget_controls']),414);self.assertTrue(all(all(r['checks'].values()) for r in self.r['optimizer']['budget_controls']))
    def test_stress_stationary_full_infimum(self):
        r=b.full_floor(b.toy());self.assertEqual(r['active_variable_branch'],'stress_stationary');self.assertEqual(r['closure_t'],b.quad(F(1,2)));self.assertEqual(r['infimum'],b.quad(9))
    def test_derivative_row_raises_full_floor_above_stress_floor(self):
        r=b.full_floor(b.toy(D1=3,A1=0));self.assertEqual(r['active_variable_branch'],'derivative_stress_crossing');self.assertEqual(r['closure_t'],b.quad(F(1,2)));self.assertEqual(r['infimum'],b.quad(F(9,2)))
        self.assertTrue(b.shape().stress_floor_test(0,0,1,F(9,2))['passes_necessary_strict_stress_floor']);self.assertFalse(b.decision(b.toy(D1=3,A1=0),F(9,2))['possible'])
    def test_linear_crossing_instead_of_division_by_zero(self):
        r=b.full_floor(b.toy(D0=3,D1=1,A1=0));self.assertEqual(r['closure_t'],b.quad(F(1,2)));self.assertEqual(r['infimum'],b.quad(F(9,2)))
    def test_irrational_stationary_cost_exact(self):
        r=b.full_floor(b.toy(A1=1));self.assertEqual(s.simplify(b.as_symbolic(r['infimum'])-3-2*s.sqrt(2)),0)
    def test_negative_quadratic_crossing_exact(self):
        r=b.full_floor(b.toy(D0=1,D1=7));self.assertEqual(s.simplify(b.as_symbolic(r['infimum'])-s.Rational(9,2)-7*s.sqrt(21)/6),0)
    def test_positive_quadratic_crossing_selects_first_root(self):
        r=b.full_floor(b.toy(D0=10));self.assertEqual(s.simplify(b.as_symbolic(r['closure_t'])-1+s.sqrt(6)/3),0)
        self.assertEqual(s.simplify(b.as_symbolic(r['infimum'])-12+s.sqrt(6)/3),0)
    def test_stationary_equal_derivative_branch(self):
        r=b.full_floor(b.toy(D0=F(9,2),D1=3));self.assertEqual(r['stationary_branch_comparison'],0);self.assertEqual(r['active_variable_branch'],'stress_stationary');self.assertEqual(r['infimum'],b.quad(9))
    def test_endpoint_stationary_M_two_still_strict(self):
        r=b.full_floor(b.toy(A1=0));self.assertEqual(r['closure_t'],b.quad(1));self.assertEqual(r['infimum'],b.quad(4));self.assertFalse(b.decision(b.toy(A1=0),4)['possible'])
    def test_variable_minimum_is_unattained(self):self.assertFalse(b.full_floor(b.toy())['infimum_attained_by_some_smooth_shape']);self.assertFalse(b.decision(b.toy(),9)['possible'])
    def test_dominating_constant_floor_is_attained(self):
        r=b.full_floor(b.toy(C0=10));self.assertTrue(r['infimum_attained_by_some_smooth_shape']);self.assertTrue(b.decision(b.toy(C0=10),10)['possible'])
        witness=b.shape().shape_gate(b.toy(C0=10),10);self.assertTrue(all(witness['checks'].values()));self.assertEqual(max(witness['row_values'].values()),10)
    def test_tied_constant_is_unattained(self):self.assertFalse(b.full_floor(b.toy(C0=9))['infimum_attained_by_some_smooth_shape']);self.assertFalse(b.decision(b.toy(C0=9),9)['possible'])
    def test_positive_slack_above_unattained_floor_passes(self):
        for e in (F(1,10),F(1,10**8),F(1,10**40)):self.assertTrue(b.decision(b.toy(),9+e)['possible'])
    def test_unit_row_is_included_in_constant_floor(self):
        r=b.full_floor(b.toy(C0=0));self.assertEqual(r['constant_floor'],1)
    def test_exact_crossing_polynomial_and_root_identity(self):self.symbolic('exact_derivative_stress_crossing_polynomial');self.symbolic('crossing_root_rationalization')
    def test_exact_stationary_minimum_identity(self):self.symbolic('stationary_stress_minimum')
    def test_mixed_value_cost_is_retained_in_budget_headroom(self):self.symbolic('mixed_budget_headroom_solved_without_dropping_value_cost')
    def test_raw_deficit_keeps_mesh_and_value_terms(self):self.symbolic('freezing_raw_deficit_keeps_mesh_and_value_terms')
    def test_smooth_second_derivative_slack_increases_stress_cost(self):self.symbolic('stress_cost_strictly_increases_with_second_enclosure')
    def test_zero_or_negative_budget_cannot_accept_shape(self):self.assertFalse(b.decision(b.toy(),0)['possible']);self.assertFalse(b.decision(b.toy(),-1)['possible'])
    def test_zero_mesh_headroom_excludes_all_shapes(self):
        c=dict(self.cost);gamma=b.current().deficit().certificate(self.data,self.slack)['floors']['unfrozen_cone_margin'];ell=10002;c['Cmesh']=gamma*ell**6/8
        r=b.joint_test(self.data,c,b.toy(),self.slack,10**50,ell);self.assertEqual(r['budget']['mesh_headroom'],0);self.assertFalse(r['possible'])
    def test_negative_mesh_headroom_cannot_be_paid_by_large_decay(self):
        c=dict(self.cost,Cmesh=F(10**50));r=b.joint_test(self.data,c,b.toy(),self.slack,10**100,10002)
        self.assertLess(r['budget']['mesh_headroom'],0);self.assertFalse(r['possible'])
    def test_chart_gate_is_independent_of_shape_feasibility(self):
        c=dict(self.cost,chart_distance_threshold=F(10**100));r=b.joint_test(self.data,c,b.toy(),self.slack,10**50,10002)
        self.assertTrue(r['shape_decision']['possible']);self.assertFalse(r['possible'])
    def test_old_separate_tests_can_pass_when_coupled_shape_test_fails(self):self.assertTrue(self.r['coupled']['checks']['old_background_shape_and_mesh_tests_pass_but_mixed_test_fails'])
    def test_mixed_boundary_nine_fails_and_ten_passes(self):self.assertFalse(self.r['coupled']['boundary']['possible']);self.assertTrue(self.r['coupled']['strictly_after']['possible'])
    def test_dominating_constant_can_pass_coupled_equality(self):self.assertTrue(self.r['coupled']['attainable_constant_boundary']['possible'])
    def test_required_deficit_floor_is_unattained_at_boundary(self):self.assertTrue(self.r['coupled']['checks']['raw_deficit_floor_equals_supplied_slack_but_is_unattained'])
    def test_complete_cost_floor_below_every_rational_feasible_upper_cost(self):
        rows=b.toy(D0=1,D1=7);floor=b.full_floor(rows)['infimum']
        for M,N in ((F(3,2),5),(2,5),(1000,1216)):
            self.assertGreater(b.shape().capacity(M,N),1);cost=max(b.shape().affine_value(v,M,N) for v in rows.values());self.assertGreater(b.qcmp(b.quad(cost),floor),0)
    def test_complete_global_jet_pressure_composition(self):self.assertTrue(all(self.r['full_composition']['checks'].values()))
    def test_auxiliary_smooth_witness_passes_original_mixed_gate(self):self.assertTrue(self.r['full_composition']['checks']['rational_smooth_witness_passes_original_current_frame_gate'])
    def test_pressure_and_shifted_viscosity_inputs_not_erased(self):self.assertTrue(self.r['full_composition']['checks']['full_pressure_and_shifted_viscosity_costs_retained'])
    def test_dyadic_rounding_is_conservative_with_whole_chart_shift(self):self.assertTrue(self.r['full_composition']['checks']['exact_integer_floor_whole_chart_shift_retained'])
    def test_invalid_rows_rejected(self):
        for rows in (dict(b.toy(),derivative=(0,0,0)),dict(b.toy(),stress=(0,3,0)),dict(b.toy(),stress=(0,-1,1))):
            with self.assertRaises(ValueError):b.full_floor(rows)
    def test_invalid_scale_and_band_rejected(self):
        for S,ell in ((0,10),(-1,10),(1,1),(1,True)):
            with self.assertRaises(ValueError):b.mixed_budget(self.data,self.cost,self.slack,S,ell)
            with self.assertRaises(ValueError):b.scale_floor(self.data,self.cost,b.toy(),S,ell)
    def test_invalid_mixed_costs_rejected(self):
        for key,v in (('profile_normalization',1),('Cpert',0),('Cmesh',-1),('chart_distance_threshold',-1)):
            c=dict(self.cost);c[key]=v
            with self.assertRaises(ValueError):b.mixed_budget(self.data,c,self.slack,10,10002)
            with self.assertRaises(ValueError):b.scale_floor(self.data,c,b.toy(),10,10002)
    def test_float_and_bool_inputs_rejected(self):
        for v in (True,.1):
            with self.assertRaises(TypeError):b.quad(v)
            with self.assertRaises(TypeError):b.decision(b.toy(),v)
            with self.assertRaises(TypeError):b.mixed_budget(self.data,self.cost,self.slack,v,10002)
            c=dict(self.cost,Cpert=v)
            with self.assertRaises(TypeError):b.checked_cost(c)
    def test_invalid_raw_deficit_rejected(self):
        for v in (0,-1,self.data['kappa']+1):
            with self.assertRaises(ValueError):b.mixed_budget(self.data,self.cost,v,10,10002)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'bbfe73f645de9184bc4cd4e52bd8cb8d64a099c2');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()))
    def test_frozen_record_exact_reproduction(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_bindings_preserved(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_all_unsupported_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_formal_flag_not_replaceable_by_numeric_zero(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_attainment_claim_tampering_rejected(self):
        r=deepcopy(self.record);r['reviews']['optimizer']['regimes'][0]['optimization']['infimum_attained_by_some_smooth_shape']=True
        with self.assertRaises(ValueError):audit.validate(r)
    def test_mixed_value_freezing_cost_cannot_be_erased(self):
        r=deepcopy(self.record);r['protocol']['mixed_budget']='mesh condition only'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_old_cutoff_verification_obligation_cannot_be_erased(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_cli_replay_and_refusal_to_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
