#!/usr/bin/env python3
"""Moving-normal Gaussian adaptation, explicit norms and scope regressions."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from budgets import P, ForcedBudget, fixture, product
from controls import review, raw_bound


class ForcedInverseBudgetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=review();cls.record=audit.generate()

    def control(self, part, fragment):
        rows={k:v for k,v in self.result["parts"][part]["checks"].items() if fragment in k}
        self.assertTrue(rows);self.assertTrue(all(rows.values()))

    def test_tangent_parameter_derivative_requires_normal_connection(self):
        for name in ("manufactured_velocity_is_tangent","parameter_derivative_is_not_tangent","differentiated_constraint_keeps_normal_connection"):
            self.control("moving_constraint",name)

    def test_natural_ambient_extension_can_lose_Gaussian_weight(self):
        self.control("moving_constraint","ungauged_normal_damping_cannot_preserve_arbitrary_Gaussian")

    def test_adapted_extension_preserves_tangent_generator_and_profile_rate(self):
        self.control("moving_constraint","adapted_extension_leaves_tangent_generator_unchanged")
        self.control("moving_constraint","adapted_normal_component_has_exact_profile_logarithmic_rate")

    def test_basis_free_lift_and_arbitrary_ambient_initial_data(self):
        for name in ("full_normal_lift_has_unit_pairing","normal_lift_connection_is_tangent","arbitrary_ambient_initial_vector_decomposition_exact"):
            self.control("moving_constraint",name)

    def test_omission_of_moving_normal_is_detected(self):
        self.control("moving_constraint","omitting_moving_normal_breaks_lift_tangency")
        self.control("explicit_budgets","moving_normal_not_erased_in_pressure_or_ambient_bound")

    def test_full_unprojected_normal_source_pressure_is_needed(self):
        self.control("moving_constraint","normal_source_projects_to_zero")
        self.control("moving_constraint","normal_source_pressure_does_not_vanish")

    def test_exact_complex_constraint_pressure_and_conjugate_controls_replayed(self):
        self.control("moving_constraint","inherited:")

    def test_manufactured_solution_satisfies_full_damped_equation_and_zero_entry(self):
        self.control("full_box_norm","source_computed_from_full_damped_equation")
        self.control("full_box_norm","whole_entrance_parameter_slice")

    def test_all_mixed_solution_norms_bounded_over_whole_box(self):
        self.control("full_box_norm","all_mixed_full_box_solution_norm_is_bounded")

    def test_box_coefficients_use_exact_rationals_without_floor_rounding(self):
        from controls import inherited
        J=inherited().maps().JetClass();x=J.variable(0,Q(0),n=2,order=5)
        self.assertEqual(raw_bound([x*x/Q(3)],1),Q(2,3))
        self.assertEqual(raw_bound([x*x/Q(3)],2),Q(2,3))

    def test_zero_source_values_at_endpoints_do_not_give_zero_exit(self):
        self.control("full_box_norm","zero_endpoint_source_values_do_not_give_zero_exit")

    def test_independent_rational_recurrence_at_multiple_orders_and_scales(self):
        self.control("explicit_budgets","independent_rational_recurrence")

    def test_normal_reciprocal_recurrence_retains_higher_jet_products(self):
        b=fixture()
        expected=2*(b.square(2)*b.reciprocal(0)+2*b.square(1)*b.reciprocal(1))
        self.assertEqual(b.reciprocal(2).record(),expected.record())

    def test_profile_logarithmic_rate_cannot_be_dropped_from_generator(self):
        b=fixture();removed=fixture();removed.profile_rate=lambda n:P()
        self.assertGreater(b.generator(1).evaluate(1),removed.generator(1).evaluate(1))
        self.assertGreater(b.velocity(1).evaluate(1),removed.velocity(1).evaluate(1))

    def test_tangent_and_length_costs_remain_in_ambient_norm(self):
        b=fixture();small=fixture();small.tangent=P(1);small.length=P(Q(1,2))
        self.assertGreater(b.ambient().evaluate(1),small.ambient().evaluate(1))
        self.assertGreater(b.velocity(0).evaluate(1),small.velocity(0).evaluate(1))

    def test_conjugate_input_harmonics_have_same_norms(self):
        self.control("explicit_budgets","conjugate_majorants_match")

    def test_whole_path_zero_source_has_zero_response(self):
        self.control("explicit_budgets","whole_path_zero_source_has_zero_response_and_pressure")

    def test_each_order_uses_finite_source_jets_and_lift_connection_jet(self):
        self.control("explicit_budgets","order_three_only_finitely_many_full_source_jets")
        self.control("explicit_budgets","lift_derivative_requires_at_least_first_normal_jet")

    def test_cutoff_pressure_and_potential_derivatives_retained(self):
        self.control("explicit_budgets","same_cutoff_pressure_and_potential_keep_derivative_products")
        with self.assertRaises(ValueError):fixture().cut("unknown",lambda n:P(1),0)

    def test_raw_binomial_product_not_factorial_convolution(self):
        a=lambda n:P(n+1);b=lambda n:P(n+2)
        self.assertEqual(product(a,b,2).evaluate(1),4+2*2*3+3*2)

    def test_bad_harmonics_floors_orders_and_bound_types_rejected(self):
        for m in (0,Q(1),True):
            with self.assertRaises(ValueError):fixture(harmonic=m)
        b=fixture()
        with self.assertRaises(ValueError):b.velocity(-1)
        b.source=lambda n:1
        with self.assertRaises(TypeError):b.projected_source(0)
        for floor in (0,-1):
            with self.assertRaises(ValueError):ForcedBudget(None,None,None,None,None,P(1),P(1),floor)
        with self.assertRaises(ValueError):ForcedBudget(None,None,None,None,None,P(),P(1),1)
        with self.assertRaises(TypeError):ForcedBudget(None,None,None,None,None,1,P(1),1)

    def test_five_frozen_dependencies_and_evidence_reproduce(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),5)
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_source_core_full_residual_and_blowup_not_promoted(self):
        for key,value in self.record["flags"].items():
            if not value:
                changed=deepcopy(self.record);changed["flags"][key]=True
                with self.assertRaises(ValueError):audit.validate(changed)

    def test_omitted_remaining_gate_or_historical_hash_rejected(self):
        for field in ("source_hashes","remaining_obligations"):
            changed=deepcopy(self.record)
            if isinstance(changed[field],dict):changed[field].pop(next(iter(changed[field])))
            else:changed[field].pop(0)
            with self.assertRaises(ValueError):audit.validate(changed)

    def test_existing_output_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"keep.json";p.write_text("sentinel")
            with self.assertRaises(SystemExit):audit.main(["--output",str(p)])
            self.assertEqual(p.read_text(),"sentinel")


if __name__=="__main__":unittest.main()
