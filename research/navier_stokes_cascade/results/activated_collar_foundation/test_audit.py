#!/usr/bin/env python3
"""Actual-collar widths, complete memory, higher jets, scope and provenance."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import mpmath as mp
from activation_bounds import gains,falling_operator,domain_review,formula_review,width_review,collar_review
from activation_checks import RJet,state,native,diagnostic_stats,old_collar
from audit import HERE,FLAGS,generate,validate,compare_records,parent_review,main


class ActivatedCollarTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_exact_current_parent_and_unchanged_parameters(self):
        self.assertEqual(self.record['protocol']['parent_commit'],'b1c66663d020a3920aafc58cb46aae421aeda8c1')
        self.assertTrue(self.record['protocol']['no_parameter_or_support_change'])

    def test_actual_collar_precedes_reference_taper_and_modulation(self):
        self.assertTrue(all(domain_review()['checks'].values()))

    def test_common_strip_is_strictly_inside_core_zero_free_domain(self):
        self.assertEqual(domain_review()['final_eta_strip'],'delta=rho/(800*S)')
        self.assertTrue(domain_review()['checks']['four_delta_inside_parent_zero_free_strip'])

    def test_fixed_width_gains_include_zero_and_arbitrary_orders(self):
        for r in (0,1,2,20,1000):
            self.assertEqual(gains(r)['field'],1-r);self.assertEqual(gains(r)['moment'],2-r)

    def test_invalid_derivative_orders_rejected(self):
        for r in (-1,True,Q(1,2)):
            with self.assertRaises(ValueError):gains(r)
            with self.assertRaises(ValueError):falling_operator(r)

    def test_all_order_log_radial_conversion_on_independent_monomials(self):
        for r in range(13):
            for a in (Q(0),Q(1),Q(4),Q(17,3),Q(-2)):
                expected=Q(1)
                for j in range(r):expected*=a-j
                self.assertEqual(sum(c*a**j for j,c in enumerate(falling_operator(r))),expected)

    def test_fourth_log_conversion_retains_lower_orders(self):
        self.assertEqual(falling_operator(4),[0,-6,11,-6,1])
        self.assertNotEqual(sum(c*Q(2)**j for j,c in enumerate(falling_operator(4))),Q(2)**4)

    def test_exact_exponential_and_Volterra_actual_fields_not_seed_replacement(self):
        self.assertIn('F_a=F_r*exp(J)',formula_review()['actual_fields'])
        self.assertIn('U_a=U_r-K',formula_review()['actual_fields'])
        self.assertTrue(all(formula_review()['checks'].values()))

    def test_actual_pressure_density_not_reference_density(self):
        self.assertIn('F_a(v,eta)^2',formula_review()['actual_pressure'])
        self.assertGreater(diagnostic_stats()['using_reference_pressure_increment_gap'],mp.mpf('1e-20'))

    def test_all_five_full_moment_densities_retained(self):
        self.assertEqual(set(formula_review()['full_moment_densities']),{'M','I','J','S','Pi'})
        self.assertIn('U^2-X*F^2',formula_review()['full_moment_densities']['S'])

    def test_small_excess_quadrature_refines_without_order_one_normalization(self):
        self.assertLess(diagnostic_stats()['maximum_refinement_error'],mp.mpf('1e-12'))

    def test_independent_historical_cancellation_safe_formula_matches(self):
        self.assertLess(diagnostic_stats()['maximum_old_formula_error'],mp.mpf('1e-12'))

    def test_all_42_complex_eta_Cauchy_checks_through_order_six(self):
        r=diagnostic_stats();self.assertEqual(r['parameter_checks'],42)
        self.assertLess(r['maximum_parameter_Cauchy_error'],mp.mpf('1e-40'))

    def test_complete_actual_radial_shear_jets_through_order_four(self):
        self.assertLess(diagnostic_stats()['maximum_radial_flow_error'],mp.mpf('1e-50'))

    def test_fourth_radial_derivative_is_not_uniform_in_width(self):
        self.assertGreater(diagnostic_stats()['fourth_derivative_growth'],10**5)

    def test_predicted_field_width_gain_tracks_all_three_widths(self):
        self.assertLess(diagnostic_stats()['scaled_field_ratio'],2)

    def test_predicted_pressure_width_gain_tracks_all_three_widths(self):
        self.assertLess(diagnostic_stats()['scaled_pressure_ratio'],2)

    def test_zero_incoming_pressure_shortcut_detected(self):
        self.assertGreater(diagnostic_stats()['omitted_incoming_pressure_gap'],mp.mpf('.1'))

    def test_Cp_and_Pi_keep_distinct_incoming_values_and_same_radial_density(self):
        with mp.workdps(75):
            eta=mp.mpf('.3');a=state('.01','.24',eta)
            expected=-mp.mpf('.73')+mp.mpf('.19')*eta+mp.mpf('.15')*eta**2
            self.assertLess(abs((a['Pi']-a['Cp']).c[0]-expected),mp.mpf('1e-60'))
            for r in range(1,5):self.assertEqual((a['Pi']-a['Cp']).c[r],0)

    def test_pressure_memory_retains_eta_dependent_axis_datum(self):
        with mp.workdps(75):
            def f(e):
                a=state('.01','.24',e,work_digits=75);return (a['Pi']-a['Cp']).c[0]
            self.assertLess(abs(mp.diff(f,mp.mpf('.3'))-mp.mpf('.28')),mp.mpf('1e-60'))

    def test_omitted_Cp_incoming_trace_is_also_nonzero(self):
        with mp.workdps(75):
            a=state('.01','.24',mp.mpf('.3'));b=state('.01','.24',mp.mpf('.3'),drop_incoming=True)
            self.assertLess(abs(a['Cp'].c[0]-b['Cp'].c[0]-mp.mpf('.033')),mp.mpf('1e-60'))

    def test_zero_kappa0_allowed_without_inverse_kappa(self):
        with mp.workdps(75):
            a=state('.01','.24',mp.mpf('.3'),mp.mpf(0))
            self.assertTrue(all(mp.isfinite(v) for n in ('F','U','M','Pi') for v in a[n].c))

    def test_closed_eta_endpoints_and_complex_parameters_are_allowed(self):
        with mp.workdps(75):
            for eta in (mp.mpf(-1),mp.mpf(1),mp.mpc('.3','.05')):
                a=state('.01','.24',eta)
                self.assertNotEqual(a['F'].c[0],0)
                self.assertTrue(all(mp.isfinite(v) for v in a['Pi'].c))

    def test_invalid_fixture_width_or_rescaled_s_rejected(self):
        for width,s in (('0','.24'),('.1','.24'),('.01','0'),('.01','.5')):
            with self.assertRaises(ValueError):state(width,s,mp.mpf('.3'))

    def test_radial_jet_inverse_exp_and_derivative_budget(self):
        z=RJet([1,1]);inverse=z.inverse();identity=z*inverse
        self.assertEqual(identity.c[0],1);self.assertTrue(all(v==0 for v in identity.c[1:]))
        self.assertEqual(RJet([0,1]).exp().derivative_value(4),1)
        with self.assertRaises(ValueError):RJet(0).inverse()
        with self.assertRaises(ValueError):z.derivative_value(5)

    def test_single_common_collar_certificate_has_all_five_regular_fields(self):
        self.assertEqual(collar_review()['H'],['F0','U0','Pi0','M0/X','V0/X'])
        self.assertTrue(all(collar_review()['checks'].values()))

    def test_radial_analyticity_and_uniform_width_constants_not_claimed(self):
        self.assertFalse(FLAGS['actual_profile_radially_analytic_at_flat_join'])
        self.assertFalse(FLAGS['radial_constants_uniform_as_width_or_C_vary'])
        self.assertTrue(width_review()['checks']['actual_radial_costs_include_C_to_twenty_r'])

    def test_771_inherited_checks_and_14_current_source_bindings(self):
        self.assertEqual(len(self.record['inherited_provenance']),771)
        self.assertTrue(all(self.record['inherited_provenance'].values()));self.assertEqual(len(self.record['source_hashes']),14)

    def test_frozen_record_and_JSON_roundtrip(self):
        f=json.loads((HERE/'evidence.json').read_text());validate(f);compare_records(self.record,f)
        validate(json.loads(json.dumps(self.record)))

    def test_raw_roundoff_values_not_frozen_as_exact_diagnostic_data(self):
        r=self.record['reviews']['higher_jet_complex_parameter_and_full_memory_controls']
        self.assertTrue(r['raw_roundoff_scores_tested_not_refrozen'])
        self.assertNotIn('maximum_parameter_Cauchy_error',r)

    def test_unconditional_promotion_or_scope_changes_rejected(self):
        for n,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][n]=True
                with self.assertRaises(ValueError):validate(r)

    def test_actual_pressure_formula_cannot_be_changed_behind_passed_gates(self):
        r=deepcopy(self.record);r['reviews']['complete_actual_forward_formulas']['actual_pressure']='Pi_a=0'
        with self.assertRaises(ValueError):validate(r)

    def test_generated_record_does_not_alias_cached_provenance(self):
        r=generate();d=r['inherited_provenance'];d.pop(next(iter(d)));d['invented:binding']=True
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('invented:binding',parent_review()[1])

    def test_global_obligation_or_gate_cannot_be_erased(self):
        for section in ('remaining_obligations','gates','source_hashes'):
            r=deepcopy(self.record)
            if section=='remaining_obligations':r[section].pop()
            else:r[section].pop(next(iter(r[section])))
            with self.assertRaises(ValueError):validate(r)

    def test_CLI_preserves_existing_output(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('preserve')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__=='__main__':unittest.main()
