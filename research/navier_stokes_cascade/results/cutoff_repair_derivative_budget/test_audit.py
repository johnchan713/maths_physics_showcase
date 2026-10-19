"""Derivative algebra, actual bump geometry, conditional scope and provenance."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import derivatives as w

class DerivativeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_first_cutoff_chain(self):self.assertTrue(self.r['symbolic']['checks']['cutoff_first_chain'])
    def test_second_cutoff_log_Jacobian(self):self.assertTrue(self.r['symbolic']['checks']['cutoff_second_chain_includes_log_jacobian'])
    def test_omitted_Jacobian_rejected(self):self.assertTrue(self.r['symbolic']['checks']['omitted_log_jacobian_detected'])
    def test_logistic_third_chain(self):self.assertTrue(self.r['symbolic']['checks']['logistic_third_chain'])
    def test_exact_center_step_jets(self):self.assertEqual(self.r['symbolic']['center_sigma_jets'],['8','0','-3712'])
    def test_full_bump_to_F_second(self):self.assertTrue(self.r['symbolic']['checks']['bump_to_F_second_with_both_inverse_radius_terms'])
    def test_step_first_enclosure(self):self.assertEqual(w.step_bounds()[0],32)
    def test_step_second_enclosure(self):self.assertEqual(w.step_bounds()[1],1216)
    def test_step_third_enclosure(self):self.assertLess(self.r['exact']['step_enclosure'][2],w.step_bounds()[2])
    def test_actual_bump_value(self):self.assertEqual(w.bump_bounds()[0],320)
    def test_actual_bump_first(self):self.assertEqual(w.bump_bounds()[1],121600)
    def test_actual_bump_second(self):self.assertEqual(w.bump_bounds()[2],220000000)
    def test_unit_mass_scaling_not_unscaled_annular_bump(self):self.assertEqual(w.bump_bounds(F(1,5)),tuple(v/F(1,5)**(j+1) for j,v in enumerate(w.step_bounds())))
    def test_all_five_supports(self):self.assertEqual(self.r['exact']['supports'],[['9/20','11/20'],['19/20','21/20'],['29/20','31/20'],['39/20','41/20'],['49/20','51/20']])
    def test_disjoint_supports(self):self.assertTrue(self.r['exact']['checks']['five_bumps_pairwise_disjoint'])
    def test_reserved_mean_unchanged(self):self.assertTrue(self.r['exact']['checks']['reserved_mean_I4_not_touched'])
    def test_cutoff_necessary_lower_not_upper(self):self.assertEqual(self.r['exact']['necessary_normalized_shape_bounds'],['1','4']);self.assertGreater(w.compressed_cutoff(1000)['first_norm_lower'],32)
    def test_compressed_center_first(self):self.assertEqual(w.compressed_cutoff(7)['center_slope'],-56)
    def test_compressed_second_lower(self):self.assertEqual(w.compressed_cutoff(7)['second_norm_lower'],784)
    def test_compressed_center_third(self):self.assertEqual(w.compressed_cutoff(2)['center_third'],29696)
    def test_arbitrarily_large_exact_compression(self):self.assertEqual(w.compressed_cutoff(10**30)['first_norm_lower'],8*10**30)
    def test_compressed_transition_endpoints(self):self.assertEqual(w.compressed_cutoff(1)['support_left'],0);self.assertEqual(w.compressed_cutoff(1)['support_right'],1)
    def test_cutoff_conditional_radial_cost(self):self.assertEqual(w.cutoff_bounds(2,F(1,4),3,9)['radial'],6)
    def test_cutoff_conditional_second_cost(self):self.assertEqual(w.cutoff_bounds(2,F(1,4),3,9)['radial2'],39)
    def test_width_halving_first_cost(self):self.assertEqual(w.cutoff_bounds(2,F(1,8),3,9)['radial'],2*w.cutoff_bounds(2,F(1,4),3,9)['radial'])
    def test_second_has_both_width_powers(self):self.assertGreater(w.cutoff_bounds(2,F(1,8),3,9)['radial2'],2*w.cutoff_bounds(2,F(1,4),3,9)['radial2'])
    def test_rho_inverse_costs(self):
        o=w.seminorms().extended_jets(dict(value=F(0),radial=F(0),radial2=F(0)),1,1,0,0,F(1,4),F(1,2),w.bump_bounds())
        self.assertEqual(o['value'],640);self.assertEqual(o['radial'],8*121760);self.assertEqual(o['radial2'],32*220243440)
    def test_conditional_composition(self):self.assertTrue(self.r['propagation']['checks']['composition_equals_retained_exact_jet_interface'])
    def test_pressure_memory_rows_not_erased(self):self.assertTrue(self.r['propagation']['checks']['full_pressure_memory_changes_both_normalized_rows'])
    def test_shape_cost_strictly_changes_derivative_budget(self):self.assertTrue(self.r['propagation']['checks']['larger_unknown_shape_norms_increase_derivative_budget'])
    def test_shape_cost_not_value_cost(self):self.assertTrue(self.r['propagation']['checks']['cutoff_derivatives_do_not_change_value_bound'])
    def test_no_e_amplitude_multiplier(self):self.assertTrue(self.r['propagation']['checks']['positive_order_bumps_are_additive_not_multiplied_by_e'])
    def test_invalid_compressions(self):
        for n in (0,-1,True,F(2),2.0):
            with self.assertRaises(ValueError):w.compressed_cutoff(n)
    def test_invalid_cutoff_budgets(self):
        for v in ((0,1,1,4),(1,0,1,4),(1,1,F(1,2),4),(1,1,1,3)):
            with self.assertRaises(ValueError):w.cutoff_bounds(*v)
    def test_invalid_bump_widths(self):
        for v in (0,-1):
            with self.assertRaises(ValueError):w.bump_bounds(v)
    def test_float_bool_inputs_rejected(self):
        for v in (.1,True):
            with self.assertRaises(TypeError):w.bump_bounds(v)
    def test_invalid_kept_radius(self):
        with self.assertRaises(ValueError):w.propagate(dict(),1,2,3,1,32,1216,1)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'ec57a0079d884edb2b1be30a587f98455628b1ef');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()))
    def test_historical_bindings(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_frozen_record(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_no_false_claim_promoted(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_scope_type_mutation(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_budget_tamper(self):
        r=deepcopy(self.record);r['reviews']['exact']['bump_bounds'][2]='1'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_remaining_gap_not_erased(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_replay_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
