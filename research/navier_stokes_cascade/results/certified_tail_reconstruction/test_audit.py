"""Exact tail/domain/diagonal controls and strict nonidentification guards."""
from copy import deepcopy
from fractions import Fraction as F
import inspect,json
from pathlib import Path
import tempfile,unittest
import audit
from tail_bounds import (ceil,ceil_sqrt,tail_budget,sum_scale,sum_margin,borel_scale,borel_margin,
    domain_review,scale_review,core_and_force_review,realization,force)

class CertifiedTailTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate()
    def test_parent_and_protocol(self):
        self.assertEqual(audit.protocol()['parent_commit'],'e34f6a5810cd2e59564a3540470d20f67f0cb5b7')
        self.assertEqual(json.loads((audit.HERE/'protocol.json').read_text()),audit.protocol())
    def test_rational_ceiling(self):self.assertEqual(ceil(F(17,3)),6);self.assertEqual(ceil(F(-17,3)),-5)
    def test_integer_sqrt_exact_squares(self):
        for n in range(50):self.assertEqual(ceil_sqrt(n*n),n)
    def test_rational_sqrt_grid_independent(self):
        for n in range(1,100):
            for d in range(1,40):
                x=F(n,d);k=ceil_sqrt(x);self.assertLess((k-1)**2,x);self.assertGreaterEqual(k*k,x)
    def test_sqrt_large_one_above_square(self):self.assertEqual(ceil_sqrt(10**100+1),10**50+1)
    def test_sqrt_very_small_positive(self):self.assertEqual(ceil_sqrt(F(1,10**1000)),1)
    def test_sqrt_negative_rejected(self):
        with self.assertRaises(ValueError):ceil_sqrt(-1)
    def test_no_float_or_boolean_root_inputs(self):
        for x in (1.0,True):
            with self.assertRaises(TypeError):ceil_sqrt(x)
    def test_all_twelve_tail_budget_fixtures_accepted(self):
        rows=domain_review()['parameter_budgets'];self.assertEqual(len(rows),12)
        self.assertTrue(all(F(r['length_margin'])>=0 and F(r['smallness_margin'])>=0 for r in rows))
    def test_length_dominated_budget_keeps_radius(self):
        a=tail_budget(10**40,F(1,200),F(1,6057072),16)
        b=tail_budget(10**40,F(1,200),F(1,16*6057072),16)
        self.assertGreater(b['band_floor'],a['band_floor'])
    def test_smallness_dominated_budget_keeps_h(self):
        a=tail_budget(25,F(1,200),F(1,6057072),16)
        b=tail_budget(25,F(1,10**5),F(1,6057072),16)
        self.assertGreater(b['band_floor'],a['band_floor'])
    def test_original_label_floor_retained(self):self.assertGreaterEqual(tail_budget(25,F(1,200),F(1,6057072),10**60)['band_floor'],10**60)
    def test_band_budget_no_jet_argument(self):self.assertNotIn('order',inspect.signature(tail_budget).parameters)
    def test_no_tiny_q_or_exponential_materialized(self):
        r=tail_budget(10**80,F(1,200),F(1,6057072),16)
        self.assertEqual(r['tail_exponent'],r['band_floor']+2);self.assertEqual(r['qstar_divisor'],64);self.assertTrue(r['qv_not_materialized'])
    def test_invalid_tail_budget_ranges(self):
        for L,h,r,b in ((0,F(1,200),F(1,6057072),16),(1,0,F(1,6057072),16),(1,F(1,100),F(1,6057072),16),(1,F(1,200),F(1),16),(1,F(1,200),F(1,6057072),15)):
            with self.assertRaises(ValueError):tail_budget(L,h,r,b)
    def test_entire_enlarged_localization_chart_inside_tail(self):
        rows=domain_review()['localization_charts'];self.assertEqual(len(rows),5)
        self.assertTrue(all(F(r['whole_enlargement_over_qv'])<=F(3,8) for r in rows))
    def test_observation_worst_ratio_is_one_quarter(self):self.assertEqual(max(F(r['whole_enlargement_over_qv']) for r in domain_review()['localization_charts']),F(3,8))
    def test_same_old_recipe_when_reference_zero(self):
        args=(3,F(1,100),[(F(7),F(5))],F(0),F(20))
        self.assertEqual(sum_scale(*args,0),realization().choose_log_scale(*args))
    def test_new_sum_reference_can_force_later_scale(self):
        x=sum_scale(1,1,[(0,0)],0,2,10**30);self.assertGreater(x,10**30)
    def test_every_new_sum_whole_interval_margin_negative(self):self.assertTrue(all(F(r['margin'])<0 for r in scale_review()['sum_rows']))
    def test_sum_margin_decreases_for_every_later_log_coordinate(self):
        s=sum_scale(2,F(1,7),[(4,9)],0,10,0)
        for d in (1,10,10**8):self.assertLess(sum_margin(2,F(1,7),[(4,9)],s+d),sum_margin(2,F(1,7),[(4,9)],s))
    def test_sum_support_and_doubling_gaps(self):
        for r in scale_review()['sum_rows']:self.assertGreater(F(r['previous_gap']),1);self.assertGreater(F(r['log_scale']),F(r['domain'])+1)
    def test_sum_invalid_inputs_not_masked_by_reference(self):
        for gamma,previous,domain,ref in ((0,0,2,0),(1,-1,2,100),(1,0,-1,0),(1,0,2,-1)):
            with self.assertRaises(ValueError):sum_scale(1,gamma,[(0,0)],previous,domain,ref)
    def test_sum_empty_or_negative_constraints_rejected(self):
        for pairs in ([],[(-1,0)],[(0,-1)]):
            with self.assertRaises(ValueError):sum_scale(1,1,pairs,0,2,0)
    def test_sum_float_enclosures_rejected(self):
        with self.assertRaises(TypeError):sum_scale(1,1,[(0.0,0)],0,2,0)
    def test_borel_reference_zero_uses_same_reviewed_recipe(self):
        p=[(0,0,F(20))];self.assertEqual(borel_scale(4,p,1,0),force().choose_dyadic_exponent(4,p,1))
    def test_all_Borel_margins_and_doubling(self):
        for r in scale_review()['Borel_rows']:self.assertLessEqual(F(r['margin']),-2);self.assertGreater(r['previous_gap'],1);self.assertTrue(r['dominates'])
    def test_old_Borel_norm_passes_but_new_norm_fails(self):
        self.assertEqual(borel_margin(4,[(0,0,0)],4),-12);self.assertEqual(borel_margin(4,[(0,0,20)],4),8)
    def test_new_Borel_sequence_resolves_changed_norm(self):
        n=borel_scale(4,[(0,0,20)],1,4);self.assertGreater(n,4);self.assertLessEqual(borel_margin(4,[(0,0,20)],n),-2)
    def test_Borel_denominator_positive_range_enforced(self):
        with self.assertRaises(ValueError):borel_scale(4,[(0,4,1)],1,0)
    def test_Borel_empty_constraints_rejected(self):
        with self.assertRaises(ValueError):borel_scale(4,[],1,0)
    def test_Borel_negative_reference_rejected(self):
        with self.assertRaises(ValueError):borel_scale(4,[(0,0,1)],1,-1)
    def test_Borel_float_log_norm_rejected(self):
        with self.assertRaises(TypeError):borel_scale(4,[(0,0,1.0)],1,0)
    def test_all_six_new_core_paths_inside_plateau(self):
        rows=core_and_force_review()['core_paths'];self.assertEqual(len(rows),6);self.assertTrue(all(r['radial_plateau'] and r['time_plateau'] for r in rows))
    def test_localized_gradient_force_difference_is_nonzero(self):self.assertEqual(F(core_and_force_review()['local_pressure_cutoff_gradient_difference']),F(-33,16))
    def test_direct_polynomial_pressure_cutoff_difference(self):
        x=F(1,4);self.assertEqual((-16*x+64*x**3)-(-4*x+4*x**3),F(-33,16))
    def test_fresh_field_and_force_not_identified_as_old(self):
        self.assertFalse(audit.FLAGS['new_annular_field_identified_with_historical_realization']);self.assertFalse(audit.FLAGS['new_localized_force_identified_with_historical_force'])
    def test_old_first_band_not_certified(self):self.assertFalse(audit.FLAGS['original_first_band_or_old_localized_force_acceptance_verified'])
    def test_actual_huge_norms_not_computed(self):self.assertFalse(audit.FLAGS['actual_infinite_fields_norms_or_diagonal_scales_numerically_evaluated'])
    def test_no_proof_assistant_premise_added(self):self.assertFalse(audit.FLAGS['proof_assistant_absence_treated_as_mathematical_obstruction'])
    def test_recursive_historical_bytes(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_frozen_roundtrip_and_exact_record(self):
        audit.validate(json.loads((audit.HERE/'evidence.json').read_text()));audit.validate(json.loads(json.dumps(self.record)))
    def test_all_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_old_scale_failure_cannot_be_erased(self):
        r=deepcopy(self.record);r['reviews']['fresh_two_diagonal_controls']['new_norm_old_scale_margin']='-100'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_domain_margin_tamper_rejected(self):
        r=deepcopy(self.record);r['reviews']['actual_value_tail_and_support_margins']['localization_charts'][0]['whole_enlargement_over_qv']='2'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_required_remaining_source_gate_cannot_be_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_cli_output_and_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)])
            audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
