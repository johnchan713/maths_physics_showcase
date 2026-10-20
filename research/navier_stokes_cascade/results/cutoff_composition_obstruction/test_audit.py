"""Regression and scope guards for the local composition nonimplication."""
from copy import deepcopy
import json,tempfile,unittest
from pathlib import Path
import audit
import composition as c

class CompositionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.s=cls.record['reviews']['symbolic'];cls.rows=cls.record['reviews']['numerical']['rows']
    def test_step_center_half(self):self.assertTrue(self.s['checks']['sigma_center_half'])
    def test_step_center_derivative_eight(self):self.assertTrue(self.s['checks']['sigma_center_slope_eight'])
    def test_target_center_exact(self):self.assertTrue(self.s['checks']['exact_s_center'])
    def test_target_first_exact_and_unbounded(self):self.assertTrue(self.s['checks']['exact_s_first']);self.assertIn('n',self.s['first_target'])
    def test_positive_bracket_has_no_hidden_cancellation(self):self.assertTrue(self.s['checks']['exact_positive_bracket'])
    def test_full_inverse_derivative_not_fixed_theta(self):self.assertTrue(self.s['checks']['exact_moving_inverse_derivative'])
    def test_frozen_phase_is_different(self):self.assertTrue(self.s['checks']['frozen_angle_wrong_by_nonunit_factor'])
    def test_lift_density_exact(self):self.assertTrue(self.s['checks']['exact_lift_density'])
    def test_lift_mean_one_exact(self):self.assertTrue(self.s['checks']['exact_mean_one_lift'])
    def test_primitive_source_exact(self):self.assertTrue(self.s['checks']['exact_anchored_primitive'])
    def test_primitive_inverse_identity(self):self.assertTrue(self.s['checks']['primitive_inverse_identity'])
    def test_full_fourth_chain(self):self.assertTrue(self.s['checks']['exact_fourth_chain_all_terms'])
    def test_fourth_principal_coefficient(self):self.assertTrue(self.s['checks']['highest_shape_coefficient'])
    def test_missing_mixed_term_is_rejected(self):self.assertTrue(self.s['checks']['missing_cross_term_rejected'])
    def test_quadrature_cases_cover_two_scales_and_three_compressions(self):self.assertEqual([(r['delta'],r['n']) for r in self.rows],[(d,n) for d in ('0.0625','0.00390625') for n in (1,4,16)])
    def test_independent_moving_phase_quadratures(self):self.assertTrue(all(r['checks']['independent_quadrature_inverse_derivative'] for r in self.rows))
    def test_all_frozen_phase_controls_fail(self):self.assertTrue(all(r['checks']['frozen_phase_rejected'] for r in self.rows))
    def test_primitive_lower_bound_checks(self):self.assertTrue(all(r['checks']['primitive_oscillation_lower_bound'] for r in self.rows))
    def test_positive_inverse_domain(self):self.assertTrue(all(r['checks']['positive_phase_density'] and r['checks']['fixed_external_phase_in_period'] for r in self.rows))
    def test_numerics_not_actual_or_interval(self):self.assertIs(self.record['reviews']['numerical']['manufactured_floating_not_interval_or_actual_profile'],True)
    def test_tolerance_not_frozen_roundoff(self):self.assertTrue(all(r['digits']==80 and r['relative_error_limit']=='1e-20' and 'relative_error' not in r for r in self.rows))
    def test_invalid_compressions(self):
        for n in (0,-1,True,1.5):
            with self.assertRaises(ValueError):c.validate_parameters('0.5',n)
    def test_invalid_delta(self):
        for d in (0,-1,2,True):
            with self.assertRaises(ValueError):c.validate_parameters(d,1)
    def test_uniform_native_and_primitive_inequalities(self):self.assertTrue(all(self.record['reviews']['inequalities']['checks'].values()))
    def test_actual_hierarchy_not_materialized(self):self.assertIs(self.record['reviews']['inequalities']['actual_G_H_N_or_deltaL_materialized'],False)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'117160fffb9824fb26831b72ab8ce03638eba082');self.assertEqual(json.loads((audit.HERE/'protocol.json').read_text()),audit.protocol())
    def test_all_historical_source_bindings(self):self.assertIs(self.record['inherited_provenance']['all_bindings_verified'],True);self.assertGreater(self.record['inherited_provenance']['binding_count'],2000)
    def test_parent_frozen_record_retained(self):self.assertIn('results/fourth_loop_shape_gate/evidence.json',self.record['inherited_provenance']['frozen_records'])
    def test_history_manifest_fixed_not_a_future_checkpoint_glob(self):
        self.assertEqual(len(audit.HISTORICAL_RECORDS),128);self.assertEqual(len(set(audit.HISTORICAL_RECORDS)),128)
        self.assertNotIn(audit.HERE.name,audit.HISTORICAL_RECORDS)
        self.assertEqual(set(self.record['inherited_provenance']['frozen_records']),{'results/'+n+'/evidence.json' for n in audit.HISTORICAL_RECORDS})
    def test_exact_frozen_evidence_reproduction(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_no_actual_profile_or_H_failure_promotion(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_zero_cannot_replace_false(self):
        r=deepcopy(self.record);r['flags']['actual_primitive_H_bound_disproved']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_composition_scope_cannot_be_removed(self):
        for k in ('family','cutoff','bounds','moving','fourth','memory','scope'):
            r=deepcopy(self.record);r['protocol'][k]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_actual_obligation_not_discharged(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_historical_hash_mutation_rejected(self):
        r=deepcopy(self.record);r['inherited_provenance']['frozen_records']['results/fourth_loop_shape_gate/evidence.json']='0'*64
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_output_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
