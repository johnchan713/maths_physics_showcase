#!/usr/bin/env python3
"""Full localization, endpoint, extension, energy and conditional scope checks."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from endpoint_bounds import (borel_log2_scale,borel_scale_review,domain_review,
                             endpoint_interface_review,jet_requirements)
from localization_checks import localization_review,support_axis_review
from force_checks import energy_review,heat_endpoint_review,terminal_gluing_review
from audit import FLAGS,HERE,compare_records,expected_protocol,inherited_provenance,main,validate


class InterfaceAndScaleTests(unittest.TestCase):
    def test_fixed_domain_margin_is_positive(self):
        row=domain_review()
        self.assertEqual(row['support_q_over_qstar_bound'],F(3,8))
        self.assertEqual(row['remaining_boundary_margin'],F(1,8))
        self.assertTrue(all(row['checks'].values()))

    def test_force_needs_three_potential_derivatives_and_Cauchy_needs_one_more(self):
        for m in (0,1,12,100):
            row=jet_requirements(m)
            self.assertEqual(row['potential_order'],m+3)
            self.assertEqual(row['uniform_time_Cauchy_potential_order'],m+4)
            self.assertEqual(row['uniform_time_Cauchy_pressure_order'],m+2)

    def test_endpoint_uses_existing_uniform_physical_estimates(self):
        row=endpoint_interface_review()
        self.assertTrue(row['no_new_boundary_label_selection'])
        self.assertTrue(all(row['checks'].values()))
        self.assertFalse(FLAGS['inherited_foundation_and_global_auxiliary_geometry_independently_verified'])

    def test_extreme_terminal_constants_fit_the_diagonal_scale_recipe(self):
        row=borel_scale_review()
        self.assertTrue(all(row['checks'].values()))
        self.assertGreater(row['rows'][-1]['dyadic_scale_power'],100000)
        self.assertFalse(row['actual_terminal_jets_evaluated'])

    def test_each_fixed_mixed_derivative_is_eventually_controlled(self):
        for a,m in ((0,0),(1,3),(20,4),(6,70)):
            first=max(1,2*(a+m))
            for j in (first,first+1,first+100):
                self.assertLessEqual(a+m,j//2)
                self.assertLess(m-j,0)

    def test_invalid_extension_constraints_or_derivative_orders_are_rejected(self):
        for args in ((0,[(0,0,0)]),(1,[(0,1,0)]),(2,[(0,0,-1)]),(2,[(0,0,float('inf'))])):
            with self.assertRaises(ValueError):
                borel_log2_scale(*args)
        for value in (-1,True,1.5):
            with self.assertRaises(ValueError):
                jet_requirements(value)


class PhysicalLocalizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.row=localization_review()
        cls.axis=support_axis_review()

    def test_complete_twelve_term_localization_identity(self):
        self.assertEqual(len(self.row['retained_term_sizes']),12)
        self.assertTrue(self.row['checks']['full_localization_identity_with_every_term'])

    def test_cylindrical_and_cartesian_residuals_agree_at_positive_viscosities(self):
        self.assertEqual(self.row['positive_viscosities'],[.4,1.,1.7])
        self.assertTrue(self.row['checks']['independent_cartesian_residual_agrees'])

    def test_curl_cutoff_matches_two_velocity_derivatives(self):
        self.assertTrue(self.row['checks']['all_velocity_jets_through_order_two_match_product_rule'])

    def test_complete_field_is_divergence_free(self):
        self.assertTrue(self.row['checks']['exact_localized_divergence'])

    def test_naive_velocity_cutoff_fails(self):
        self.assertTrue(self.row['checks']['omitted_curl_cutoff_breaks_divergence'])

    def test_theta_dependent_cutoff_fails_for_direct_swirl(self):
        self.assertTrue(self.row['checks']['axisymmetry_is_required_for_direct_swirl'])

    def test_each_retained_term_has_a_nonzero_omission_control(self):
        self.assertTrue(self.row['checks']['each_retained_term_is_nonzero'])
        for name in ('pressure_cutoff','gradient_viscosity','curl_viscosity','curl_square','cutoff_time','curl_time'):
            self.assertGreater(self.row['retained_term_sizes'][name],1e-7)

    def test_discarded_transport_term_is_identically_zero_by_orthogonality(self):
        self.assertTrue(self.row['checks']['d_dot_grad_c_is_exactly_zero'])

    def test_nonzero_Cartesian_axis_is_regular_and_solenoidal(self):
        self.assertTrue(self.axis['checks']['nonzero_smooth_axis_velocity'])
        self.assertTrue(self.axis['checks']['Cartesian_axis_has_zero_divergence'])

    def test_spatial_support_initialization_and_inner_plateaus(self):
        self.assertTrue(self.axis['checks']['spatial_and_initial_zero_extensions_are_exact'])
        self.assertTrue(self.axis['checks']['inner_cutoff_plateau_preserves_actual_field'])


class EndpointAndForceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.heat=heat_endpoint_review()
        cls.terminal=terminal_gluing_review()

    def test_gamma_integral_gives_bounded_endpoint_jets(self):
        self.assertTrue(self.heat['checks']['every_fixed_heat_derivative_has_uniform_nonnegative_argument_bound'])

    def test_heat_equation_uses_tau_with_the_correct_time_sign(self):
        self.assertTrue(self.heat['checks']['exact_heat_identity_includes_tau_equals_one_minus_t'])
        self.assertTrue(self.heat['checks']['reversed_time_sign_fails'])

    def test_exterior_pressure_majorant_is_integrable(self):
        self.assertTrue(self.heat['checks']['pressure_majorant_is_integrable_for_positive_radius'])

    def test_all_available_mixed_terminal_jets_match(self):
        self.assertTrue(self.terminal['checks']['all_available_mixed_terminal_jets_match'])

    def test_time_cutoff_product_rule_is_checked_independently(self):
        self.assertTrue(self.terminal['checks']['independent_scalar_time_product_rule_agrees'])

    def test_zero_future_force_would_jump(self):
        self.assertTrue(self.terminal['checks']['extending_the_force_by_zero_would_create_a_jump'])

    def test_origin_flatness_alone_does_not_supply_endpoint_limits(self):
        self.assertTrue(self.terminal['checks']['flatness_at_origin_alone_is_not_endpoint_regularity'])

    def test_extension_preserves_compact_support_and_flat_origin(self):
        self.assertTrue(self.terminal['checks']['spatial_support_time_support_and_flat_origin_are_preserved'])


class EnergyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.row=energy_review()

    def test_compact_manufactured_field_is_exactly_solenoidal(self):
        self.assertTrue(self.row['checks']['manufactured_field_is_solenoidal'])

    def test_energy_identity_keeps_dissipation(self):
        self.assertTrue(self.row['checks']['complete_forced_energy_identity'])
        self.assertTrue(self.row['checks']['omission_of_dissipation_is_detected'])

    def test_pressure_and_nonlinear_fluxes_integrate_to_zero(self):
        self.assertTrue(self.row['checks']['pressure_and_nonlinear_integrals_cancel'])
        self.assertTrue(self.row['checks']['nonlinear_energy_density_is_not_trivially_zero'])

    def test_quadrature_refines_failed_coarse_data_at_unchanged_tolerance(self):
        self.assertEqual(self.row['quadrature_nodes_per_axis'],96)
        self.assertGreater(self.row['retained_failed_64_node_integration_by_parts_error'],2e-9)
        self.assertLess(self.row['integration_by_parts_error'],2e-9)

    def test_whole_space_integrated_energy_bound(self):
        self.assertTrue(self.row['checks']['integrated_energy_bound_holds'])


class EvidenceTests(unittest.TestCase):
    def test_protocol_and_all_inherited_source_bytes(self):
        self.assertEqual(json.loads((HERE/'protocol.json').read_text()),expected_protocol())
        inherited=inherited_provenance()
        self.assertEqual(len(inherited),382)
        self.assertTrue(all(inherited.values()))

    def test_conditional_force_is_distinct_from_unconditional_solution_or_comparison(self):
        record=json.loads((HERE/'evidence.json').read_text())
        validate(record)
        self.assertTrue(record['flags']['conditional_all_order_compact_smooth_force_extension_supplied'])
        self.assertFalse(record['flags']['arbitrary_global_smooth_bounded_energy_solution_comparison_verified'])
        self.assertFalse(record['flags']['blowup_verified'])

    def test_premature_comparison_foundation_or_blowup_promotion_is_rejected(self):
        for flag in ('blowup_verified','unconditional_admissible_force_and_solution_verified',
                     'arbitrary_global_smooth_bounded_energy_solution_comparison_verified',
                     'inherited_foundation_and_global_auxiliary_geometry_independently_verified'):
            record=json.loads((HERE/'evidence.json').read_text())
            record['flags'][flag]=True
            with self.assertRaises(ValueError):
                validate(record)

    def test_scope_or_exact_domain_budget_tampering_is_rejected(self):
        record=json.loads((HERE/'evidence.json').read_text())
        bad=copy.deepcopy(record)
        bad['reviews']['compact_localization_domain']['support_q_over_qstar_bound']='1/2'
        with self.assertRaises(ValueError):
            compare_records(record,bad)
        bad=copy.deepcopy(record)
        bad['conditional_proof']['arbitrary_global_solution_comparison_not_assumed']=False
        with self.assertRaises(ValueError):
            validate(bad)

    def test_frozen_replay_and_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'replay.json'
            main(['--verify-record',str(HERE/'evidence.json'),'--output',str(output)])
            self.assertTrue(output.exists())
            with self.assertRaises(SystemExit):
                main(['--output',str(output)])


if __name__=='__main__':
    unittest.main()
