#!/usr/bin/env python3
"""Canonical-root recurrences, independent jets and scoped evidence tests."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from budgets import P,RootBudget,fixture
from controls import Jet,exact_fixture,matmul,matrix_norm,review


class CanonicalRootBudgetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate()

    def controls(self,part,fragment):
        rows={k:v for k,v in self.result["parts"][part]["checks"].items() if fragment in k}
        self.assertTrue(rows);self.assertTrue(all(rows.values()))

    def test_noncommuting_matrix_inverse_jets_and_b_bounds(self):
        self.controls("matrix_roots_signed","inverse_identity")
        self.controls("matrix_roots_signed","noncommuting_inverse_majorant")
        self.controls("matrix_roots_signed","weighted_target_to_b_majorant")

    def test_root_and_reciprocal_jets_match_independent_fractional_expansions(self):
        for fragment in ("root_square","reciprocal_product","root_majorant","reciprocal_majorant"):
            self.controls("matrix_roots_signed",fragment)

    def test_actual_root_recurrence_contains_all_compositions(self):
        b=RootBudget(lambda n:P(),lambda n:P(2 if n else 1),lambda n:P(1),
                     lambda n:P(1),lambda n:P(),1,1,1,1,1)
        self.assertEqual(b.reciprocal(1).evaluate(1),1)
        self.assertEqual(b.reciprocal(2).evaluate(1),Q(5,2))
        self.assertEqual(b.reciprocal(3).evaluate(1),Q(13,2))

    def test_differentiated_inverse_cannot_be_replaced_by_constant_inverse(self):
        self.controls("matrix_roots_signed","omitted_inverse_derivative_detected")

    def test_differentiated_original_reciprocal_cannot_be_omitted(self):
        self.controls("matrix_roots_signed","omitted_reciprocal_root_derivative_detected")

    def test_same_root_signed_cross_is_exact_and_square_is_retained(self):
        self.controls("matrix_roots_signed","signed_original_root_cross_exact")
        self.controls("matrix_roots_signed","signed_square_must_remain")

    def test_factorial_conversion_keeps_mixed_derivatives(self):
        self.controls("factorial_curl_flat","ordinary_derivative_conversion")
        self.controls("edge_ledger","three_component_Euclidean_padding")
        j=Jet({(2,1):3});self.assertEqual(j.norm(3),3)
        self.assertEqual(j.diff(0).diff(0).diff(1).c,{(0,0):Q(6)})

    def test_first_order_curl_uses_all_coordinate_derivatives(self):
        self.controls("factorial_curl_flat","all_coordinate_derivative_sum")
        self.controls("factorial_curl_flat","full_variable_curl_coefficient_product")

    def test_variable_curl_coefficients_are_differentiated(self):
        self.controls("factorial_curl_flat","derivative_of_curl_coefficient_cannot_be_omitted")

    def test_flat_factor_exponential_jet_bound(self):
        self.controls("factorial_curl_flat","flat_exponential_Bell_majorant")

    def test_root_and_curl_edge_exponents_cover_every_partition(self):
        for fragment in ("reciprocal_edge_degree","root_edge_degree","curl_edge_degree","flat_product_edge_degree"):
            self.controls("edge_ledger",fragment)

    def test_width_and_actual_haar_radius_constants_retained(self):
        self.controls("edge_ledger","fixed_actual_radius_mass_padding")
        self.controls("edge_ledger","positive_edge_unweighting")
        self.controls("edge_ledger","same_flat_constants_retained")

    def test_normalized_outer_root_not_declared_bounded(self):
        self.controls("edge_ledger","outer_normalized_root_need_not_be_bounded")
        self.assertTrue(self.record["flags"]["normalized_outer_target_allowed_polynomial_edge_pole"])

    def test_unweighted_native_bounds_feed_atlas_exactly_once(self):
        self.controls("aggregate_pressure_chain","native_to_atlas_without_double_unweighting")
        self.controls("aggregate_pressure_chain","all_positive_scale_coefficients")

    def test_initial_scalar_chain_accepts_derived_native_bounds(self):
        self.controls("aggregate_pressure_chain","scalar_initial_chain_accepts_derived_root_bounds")
        self.controls("aggregate_pressure_chain","scalar_underlying_jet_cost")

    def test_whole_moving_pressure_moment_keeps_extra_underlying_jet(self):
        self.controls("aggregate_pressure_chain","whole_moving_pressure_moment_positive")
        self.controls("aggregate_pressure_chain","whole_moment_underlying_jet_cost")
        self.controls("aggregate_pressure_chain","whole_moment_curl_coefficient_cost")

    def test_signed_polynomial_budget_remains_positive(self):
        self.controls("aggregate_pressure_chain","signed_budget_is_positive")

    def test_invalid_floor_width_and_nonexact_input_rejected(self):
        data=(lambda n:P(1),)*5
        for constants in ((1,0,1,1,1),(1,1,Q(1,2),1,1),(1,1,1,0,1),(1,1,1,1,0)):
            with self.assertRaises(ValueError):RootBudget(*data,*constants)
        with self.assertRaises(TypeError):RootBudget(*data,1,1,1,1,.5)
        with self.assertRaises(ValueError):fixture().root(-1)

    def test_invalid_bound_data_not_silently_promoted(self):
        b=RootBudget(lambda n:1,*(lambda n:P(1),)*4,1,1,1,1,1)
        with self.assertRaises(TypeError):b.K(1)

    def test_positive_polynomials_evaluate_consistently(self):
        b=fixture()
        for n in range(4):
            self.assertEqual(b.K(n).evaluate(4),2*sum(b.H(j).evaluate(4)*b.K(n-j).evaluate(4) for j in range(1,n+1)) if n else 2)

    def test_all_frozen_source_bindings_and_record_reproduce(self):
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),5)
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_no_ode_source_residual_or_blowup_promotion(self):
        for key,value in self.record["flags"].items():
            if not value:
                changed=deepcopy(self.record);changed["flags"][key]=True
                with self.assertRaises(ValueError):audit.validate(changed)

    def test_tampered_recurrence_or_historical_binding_rejected(self):
        for field in ("source_hashes","dependency_evidence_sha256"):
            changed=deepcopy(self.record);changed[field].pop(next(iter(changed[field])))
            with self.assertRaises(ValueError):audit.validate(changed)

    def test_evidence_output_cannot_overwrite_existing_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"keep.json";p.write_text("sentinel")
            with self.assertRaises(SystemExit):audit.main(["--output",str(p)])
            self.assertEqual(p.read_text(),"sentinel")


if __name__=="__main__":unittest.main()
