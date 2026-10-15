#!/usr/bin/env python3
"""Pressure recovery, all flux powers, exact difference identities and scope."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from comparison_bounds import (cutoff_weights,flux_review,gronwall_review,holder_review,kernel_exponents,
                               kernel_review,missing_hypothesis_review,pressure_recovery_review,
                               young_absorption_review,young_remainder)
from pressure_checks import (commutator_review,fourier_pressure_review,negative_sobolev_review,weighted_energy_review)
from audit import FLAGS,HERE,compare_records,expected_protocol,inherited_provenance,main,validate


class ExactExponentTests(unittest.TestCase):
    def test_all_Holder_and_interpolation_exponents(self):
        self.assertTrue(all(holder_review()['checks'].values()))

    def test_eighth_power_supplies_the_commutator_dual_weight(self):
        self.assertTrue(cutoff_weights(8)['valid'])
        self.assertEqual(cutoff_weights(8)['margins']['commutator_pressure_dual_weight'],0)
        self.assertFalse(cutoff_weights(7)['valid'])
        self.assertFalse(cutoff_weights(4)['valid'])

    def test_larger_cutoff_powers_retain_all_weights(self):
        for power in (9,12,32):
            self.assertTrue(cutoff_weights(power)['valid'])

    def test_kernel_near_far_and_norm_radius_exponents(self):
        row=kernel_exponents()
        self.assertEqual(row['near_radial_power'],-F(2,3))
        self.assertEqual(row['far_radial_power'],-2)
        self.assertEqual(row['near_integral_coefficient'],3)
        self.assertEqual(row['far_integral_coefficient'],1)
        self.assertEqual(row['norm_radius_power'],-F(3,4))

    def test_kernel_membership_requires_both_endpoint_restrictions(self):
        for p in (1,F(3,2),2):
            with self.assertRaises(ValueError):
                kernel_exponents(p)

    def test_independent_improper_kernel_integrals(self):
        self.assertTrue(all(kernel_review()['checks'].values()))

    def test_all_three_Young_remainder_powers(self):
        self.assertEqual(young_remainder(F(3,2),1),4)
        self.assertEqual(young_remainder(F(1,2),1),F(4,3))
        self.assertEqual(young_remainder(F(3,4),F(7,4)),F(14,5))

    def test_exact_Young_maximizers_have_radius_independent_constants(self):
        self.assertTrue(all(young_absorption_review()['checks'].values()))

    def test_a_quadratic_or_superquadratic_flux_is_not_absorbed_by_this_lemma(self):
        for a,b in ((2,1),(3,1),(F(1,2),0)):
            with self.assertRaises(ValueError):
                young_remainder(a,b)

    def test_every_offset_and_flux_has_at_least_four_thirds_radius_decay(self):
        row=flux_review()
        self.assertEqual(row['uniform_absorbed_remainder_power'],F(4,3))
        self.assertTrue(all(row['checks'].values()))

    def test_missing_cancellation_integrals_are_positive_and_diverge(self):
        row=missing_hypothesis_review()
        self.assertEqual(row['uncanceled_near_singularity_integrals'],[9,999])
        self.assertTrue(all(row['checks'].values()))

    def test_raw_L1_pressure_shortcut_fails(self):
        row=missing_hypothesis_review()
        shells=row['pressure_L1_shells']
        self.assertGreater(shells[-1]['pressure_L1_shell_mass'],3*shells[0]['pressure_L1_shell_mass'])


class PressureRecoveryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.radial=negative_sobolev_review()
        cls.fourier=fourier_pressure_review()
        cls.commutator=commutator_review()

    def test_negative_Sobolev_bounds_use_only_L1_tensors(self):
        self.assertTrue(all(self.radial['checks'].values()))
        self.assertTrue(all(pressure_recovery_review()['checks'].values()))

    def test_time_averaged_momentum_has_the_correct_integration_by_parts_sign(self):
        minus_a_times_wprime=-2*(F(1,3)-F(1,4))
        a_prime_times_w=F(1,3)-F(1,2)
        self.assertEqual(minus_a_times_wprime,a_prime_times_w)
        self.assertNotEqual(minus_a_times_wprime,-a_prime_times_w)

    def test_full_old_new_and_new_square_tensor(self):
        self.assertTrue(self.fourier['checks']['full_old_new_and_new_square_tensor'])

    def test_pressure_Poisson_multiplier_sign_and_deliberate_reversal(self):
        self.assertTrue(self.fourier['checks']['pressure_multiplier_has_correct_sign'])
        self.assertTrue(self.fourier['checks']['reversed_pressure_sign_fails'])

    def test_conservative_and_advective_difference_equations_agree(self):
        self.assertTrue(self.fourier['checks']['conservative_and_advective_forms_agree'])

    def test_instantaneous_difference_momentum_and_solenoidal_time_jet(self):
        self.assertTrue(self.fourier['checks']['exact_difference_momentum_instantaneous_jet'])
        self.assertTrue(self.fourier['checks']['pressure_makes_time_derivative_solenoidal'])

    def test_commutator_matches_independent_single_mode_symbol_convolution(self):
        self.assertTrue(self.commutator['checks']['physical_product_and_independent_symbol_convolution_agree'])

    def test_omitted_commutator_is_detected(self):
        self.assertTrue(self.commutator['checks']['dropping_the_commutator_fails'])

    def test_zero_output_pressure_mode_is_retained_explicitly(self):
        self.assertTrue(self.commutator['zero_output_modes_included'])

    def test_finite_energy_is_essential_for_the_pressure_harmonic_argument(self):
        bad=self.radial['sharp_infinite_energy_counterexample']
        self.assertEqual(bad['momentum_residual'],[0,0,0])
        self.assertEqual(bad['initial_velocity'],'0')
        self.assertTrue(self.radial['checks']['smoothness_without_finite_energy_has_a_counterexample'])


class WeightedComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.energy=weighted_energy_review()

    def test_complete_weighted_energy_identity(self):
        self.assertTrue(self.energy['checks']['complete_compact_test_energy_identity'])

    def test_pressure_transport_laplacian_and_strain_terms_are_nonzero(self):
        self.assertEqual(set(self.energy['retained_terms']),{'strain','laplacian_flux','transport_flux','pressure_flux'})
        self.assertTrue(self.energy['checks']['all_flux_and_strain_terms_are_nonzero'])

    def test_a_time_dependent_constant_pressure_difference_has_zero_flux(self):
        self.assertTrue(self.energy['checks']['pressure_constant_has_zero_flux'])

    def test_replacing_v_by_w_requires_actual_support_separation(self):
        self.assertTrue(self.energy['checks']['premature_transport_replacement_fails'])

    def test_dissipation_is_retained(self):
        self.assertTrue(self.energy['checks']['dissipative_term_is_nonzero'])

    def test_fixed_T_radius_limit_and_conditional_exclusion(self):
        self.assertTrue(all(gronwall_review()['checks'].values()))
        self.assertTrue(FLAGS['no_global_derivative_or_dissipation_bound_on_competitor_assumed'])

    def test_periodic_diagnostics_are_not_promoted_to_whole_space_infinity_proof(self):
        self.assertTrue(self.energy['periodic_weight_identity_only'])
        self.assertTrue(FLAGS['periodic_diagnostics_are_not_whole_space_infinity_estimates'])


class EvidenceTests(unittest.TestCase):
    def test_protocol_and_all_inherited_source_bytes(self):
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()),expected_protocol())
        inherited=inherited_provenance()
        self.assertEqual(len(inherited),394)
        self.assertTrue(all(inherited.values()))

    def test_general_comparison_and_conditional_exclusion_do_not_promote_the_foundation(self):
        record=json.loads((HERE/'evidence.json').read_text())
        validate(record)
        self.assertTrue(record['flags']['general_smooth_compact_solution_comparison_proof_supplied'])
        self.assertTrue(record['flags']['conditional_global_smooth_bounded_energy_solution_exclusion_supplied'])
        self.assertFalse(record['flags']['inherited_foundation_and_global_auxiliary_geometry_independently_verified'])
        self.assertFalse(record['flags']['blowup_verified'])

    def test_unconditional_or_unforced_blowup_promotion_is_rejected(self):
        for flag in ('blowup_verified','unforced_blowup_verified','unconditional_admissible_construction_and_exclusion_verified',
                     'inherited_foundation_and_global_auxiliary_geometry_independently_verified'):
            row=json.loads((HERE/'evidence.json').read_text())
            row['flags'][flag]=True
            with self.assertRaises(ValueError):
                validate(row)

    def test_exact_radius_exponent_or_scope_tampering_is_rejected(self):
        row=json.loads((HERE/'evidence.json').read_text())
        bad=copy.deepcopy(row)
        bad['reviews']['flux_absorption_and_radius_decay']['reviews']['exponents']['uniform_absorbed_remainder_power']='1'
        with self.assertRaises(ValueError):
            compare_records(row,bad)
        bad=copy.deepcopy(row)
        bad['conditional_proof']['independent_profile_or_atlas_validation_not_supplied']=False
        with self.assertRaises(ValueError):
            validate(bad)

    def test_frozen_replay_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'replay.json'
            main(['--verify-record',str(HERE/'evidence.json'),'--output',str(output)])
            self.assertTrue(output.exists())
            with self.assertRaises(SystemExit):
                main(['--output',str(output)])


if __name__=='__main__':
    unittest.main()
