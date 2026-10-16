#!/usr/bin/env python3
"""Current comparison, exact forced-pressure distinction, limits and scope."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
from comparison import exponent_review,fixed_interval_bound,core_threshold
from polynomial_checks import Poly,gradient,laplacian,forced_difference_review
from audit import HERE,FLAGS,generate,validate,compare_records,main,parent_review,general_review


class CurrentComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def checks(self,n):return self.record['reviews'][n]['checks']
    def historical(self,n):return general_review()['controls'][n]

    def test_current_parent_and_not_historical_field(self):
        self.assertEqual(self.record['protocol']['parent_commit'],'236056afbfd2d20e1839733708fdc53371835bc4')
        self.assertTrue(general_review()['historical_field_not_substituted_for_current'])

    def test_same_force_datum_complete_pressure_and_energy_interface(self):
        self.assertTrue(all(self.checks('current_hypothesis_interface').values()))
        self.assertIn('force_divergence_in_exact_pressure_equation_retained',parent_review()[0]['flags'])

    def test_full_polynomial_derivatives_are_exact(self):
        x=Poly({(1,0,0):1});y=Poly({(0,1,0):1})
        self.assertEqual(laplacian((x*x+y*y)*(x*x+y*y)),16*(x*x+y*y))
        self.assertEqual(gradient(3*x*x)[0],6*x)

    def test_both_velocities_and_time_jets_solenoidal(self):
        c=forced_difference_review()['checks']
        self.assertTrue(c['both_velocities_solenoidal']);self.assertTrue(c['both_instantaneous_time_jets_solenoidal'])

    def test_common_force_divergence_nonzero(self):
        self.assertTrue(forced_difference_review()['checks']['common_force_has_nonzero_divergence'])
        self.assertEqual(forced_difference_review()['exact_force_divergence'],'4')

    def test_both_individual_forced_pressure_identities(self):
        c=forced_difference_review()['checks']
        self.assertTrue(c['individual_forced_pressure_u']);self.assertTrue(c['individual_forced_pressure_v'])

    def test_unforced_formula_for_individual_pressure_fails(self):
        self.assertTrue(forced_difference_review()['checks']['individual_unforced_pressure_shortcut_fails'])

    def test_common_force_cancels_exactly_in_difference_pressure(self):
        self.assertTrue(forced_difference_review()['checks']['force_cancels_in_difference_pressure'])

    def test_subtracting_force_divergence_twice_fails(self):
        self.assertTrue(forced_difference_review()['checks']['retaining_force_divergence_again_in_difference_is_wrong'])

    def test_full_advective_difference_and_nontrivial_pressure(self):
        c=forced_difference_review()['checks']
        for n in ('full_advective_difference_equation','difference_pressure_gradient_nonzero','difference_stress_double_divergence_nonzero'):self.assertTrue(c[n])

    def test_polynomials_not_claimed_as_finite_energy_trajectories(self):
        self.assertTrue(forced_difference_review()['no_global_finite_energy_trajectory_claim'])

    def test_harmonic_pressure_and_negative_Sobolev_controls(self):
        self.assertTrue(all(self.historical('pressure_recovery_and_negative_Sobolev').values()))

    def test_independent_fourier_pressure_and_difference_algebra(self):
        self.assertTrue(all(self.historical('complete_difference_pressure_algebra').values()))

    def test_complete_weighted_energy_including_fluxes_reproduced(self):
        self.assertTrue(all(self.historical('complete_weighted_energy_identity').values()))

    def test_kernel_cancellation_far_cap_and_independent_symbols(self):
        self.assertTrue(all(self.historical('commutator_kernel_and_independent_symbols').values()))

    def test_all_failed_shortcut_controls_reproduced(self):
        self.assertTrue(all(self.historical('sharp_failed_shortcuts').values()))

    def test_exact_kernel_radial_integral_and_norm_power(self):
        r=exponent_review();self.assertEqual(r['kernel_integral_coefficient'],'4')
        self.assertEqual(r['kernel_norm_radius_power'],'-3/4')

    def test_eighth_weight_and_failed_seventh_weight(self):
        c=exponent_review()['checks'];self.assertTrue(c['eighth_weight_just_sufficient']);self.assertTrue(c['seventh_weight_insufficient'])

    def test_all_absorption_and_offset_powers(self):
        r=exponent_review();self.assertEqual([x['young_decay'] for x in r['rows']],['4','4/3','14/5'])
        self.assertTrue(r['checks']['worst_decay_is_four_thirds'])

    def test_critical_negative_Sobolev_order_not_silently_included(self):
        self.assertTrue(exponent_review()['checks']['critical_negative_Sobolev_order_not_used'])

    def test_zero_strain_Gr_onwall_has_no_division(self):
        r=fixed_interval_bound(Q(3,4),0,5,2)
        self.assertEqual(r['prefactor'],'15/2');self.assertEqual(r['exponential_argument'],'0')

    def test_large_strain_kept_symbolic_not_rounded_or_overflowed(self):
        self.assertEqual(fixed_interval_bound(Q(3,4),10**100,5,128)['exponential_argument'],str(Q(3,2)*10**100))

    def test_invalid_fixed_intervals_or_radius_rejected(self):
        for args in ((1,1,1,2),(-1,1,1,2),(Q(1,2),-1,1,2),(Q(1,2),1,-1,2),(Q(1,2),1,1,0)):
            with self.assertRaises(ValueError):fixed_interval_bound(*args)

    def test_current_plateau_threshold_both_constraints(self):
        self.assertEqual(core_threshold(Q(1,16),1,1),Q(1,32))
        self.assertEqual(core_threshold(Q(1,16),Q(1,1000),17),Q(1,136000))

    def test_nonpositive_core_constants_rejected(self):
        for args in ((0,1,1),(1,0,1),(1,1,-1)):
            with self.assertRaises(ValueError):core_threshold(*args)

    def test_fixed_T_then_radius_then_core_limits(self):
        self.assertTrue(all(self.checks('fixed_interval_and_protected_core_limits').values()))
        self.assertFalse(FLAGS['uniform_comparison_constants_through_terminal_time_assumed'])

    def test_748_inherited_hashes_and_ten_source_bindings(self):
        self.assertEqual(len(self.record['inherited_provenance']),748)
        self.assertTrue(all(self.record['inherited_provenance'].values()));self.assertEqual(len(self.record['source_hashes']),10)

    def test_frozen_record_and_JSON_roundtrip(self):
        r=json.loads((HERE/'evidence.json').read_text());validate(r);compare_records(self.record,r)
        validate(json.loads(json.dumps(self.record)))

    def test_unconditional_promotions_rejected(self):
        for n,v in FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][n]=True
                with self.assertRaises(ValueError):validate(r)

    def test_missing_gate_or_provenance_or_source_rejected(self):
        for section in ('gates','inherited_provenance','source_hashes'):
            r=deepcopy(self.record);r[section].pop(next(iter(r[section])))
            with self.assertRaises(ValueError):validate(r)

    def test_precise_pressure_ledger_cannot_be_mutated_behind_passed_gate(self):
        r=deepcopy(self.record);r['reviews']['current_hypothesis_interface']['pressure']='individual pressure is unforced'
        with self.assertRaises(ValueError):validate(r)

    def test_generated_record_does_not_alias_cached_reference_data(self):
        r=generate()
        control=r['reviews']['general_pressure_and_infinity_controls']['checks']
        control[next(iter(control))]=False
        with self.assertRaises(ValueError):validate(r)
        self.assertTrue(all(general_review()['checks'].values()))
        r=generate();d=r['inherited_provenance'];key=next(iter(d))
        d.pop(key);d['invented:binding']=True
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('invented:binding',parent_review()[1])

    def test_foundation_obligations_cannot_be_erased(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop()
        with self.assertRaises(ValueError):validate(r)

    def test_CLI_preserves_existing_output(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('preserve')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__=='__main__':unittest.main()
