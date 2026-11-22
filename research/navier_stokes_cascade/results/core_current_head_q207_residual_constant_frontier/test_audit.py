import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q207


class Q207Tests(unittest.TestCase):
    def test_provider_statuses_and_same_object(self):
        self.assertTrue(all(q207.statuses().values()))
        self.assertTrue(all(q207.same_object_links().values()))

    def test_rates_and_exact_stage_ceiling(self):
        r = q207.rates()
        self.assertEqual(r["delta"], Q(49997, 100000))
        self.assertEqual(r["B0"], Q(74999, 50000))
        self.assertEqual(r["Lambda"], Q(151, 100))
        sample = q207.stage_for_power()
        self.assertEqual(sample["stage"], 602)
        self.assertEqual(sample["exponent"], Q(5999, 2500000))
        self.assertEqual(sample["previous_exponent"], Q(-401, 4000000))
        self.assertTrue(sample["meets_target"])
        self.assertTrue(sample["minimal"])
        self.assertEqual(q207.stage_for_power(1, 0)["stage"], 1402)
        self.assertEqual(q207.stage_for_power(0, 1)["stage"], 1002)

    def test_constant_construction(self):
        r = q207.constant_construction()
        self.assertTrue(r["graph_finite"])
        self.assertTrue(r["graph_acyclic"])
        self.assertTrue(r["finite_monotone_composition_of_tame_and_archived_inverse_estimates"])
        self.assertFalse(r["literal_polynomial_in_all_inverse_norms_claimed"])
        self.assertTrue(r["constant_finite_if_primitive_bundle_holds"])
        self.assertFalse(r["constant_uniform_in_stage"])
        self.assertFalse(r["large_derivative_array_evaluated"])
        self.assertFalse(r["numerical_materialization_required"])

    def test_complete_residual_ledger(self):
        r = q207.residual_ledger()
        self.assertEqual(len(r["finite_state_parts"]), 4)
        self.assertEqual(len(r["six_diagonal_groups"]), 6)
        self.assertIn("unit_viscosity", r["regular_contains"])
        self.assertIn("old_new_crosses", r["regular_contains"])
        self.assertIn("square", r["regular_contains"])
        self.assertIn("source_moving_normal_pressure", r["complete_memory"])
        self.assertTrue(r["scale_cutoffs_allocated_only_in_diagonal"])
        self.assertTrue(r["potential_direct_swirl_pressure_cut_together"])
        self.assertTrue(r["cutoff_scale_chosen_after_finite_constant"])
        self.assertFalse(r["stage_uniform_constant_used"])

    def test_primitive_frontier(self):
        r = q207.primitive_frontier()
        self.assertIn("f_(j,m)/P", r["log_renormalized_inequality"])
        self.assertIn("S^-d", r["log_renormalized_inequality"])
        self.assertEqual(r["pure_power_conversion"], "for eta>0 the exponent is B_j-eta")
        self.assertTrue(r["same_exponent_pure_power_requires_d_zero"])
        self.assertTrue(r["Q196_initial_source_and_endpoint_already_derived_conditionally"])
        self.assertEqual(r["smooth_mixed_jet_zero_extension_both_endpoints"], "required")
        self.assertEqual(len(r["source_terms"]), 8)
        self.assertTrue(r["named_providers_supply_log_renormalized_bounds_conditionally"])
        self.assertFalse(r["evidence_hashes_prove_continuum_supremum"])
        self.assertFalse(r["finite_symbolic_request_proves_continuum_supremum"])
        self.assertFalse(r["new_exponent_or_algebraic_obstruction_found"])
        self.assertFalse(r["independent_primitive_continuum_proof_complete"])

    def test_concrete_stage_602_request(self):
        r = q207.concrete_request()
        self.assertEqual(r["stage"], 602)
        self.assertEqual(r["A0_integrations"], 30100002)
        self.assertEqual(r["A1_integrations"], 1)
        self.assertEqual(r["A0_strict_margin"], Q(1, 20000000))
        self.assertEqual(r["A1_strict_margin"], Q(48193, 20000000))
        self.assertEqual(r["first_real_derivative_ceiling"], 36479350300)
        self.assertEqual(r["sufficient_common_pole"], 218876101818)
        self.assertEqual(r["largest_pressure_moment_archive_pole"], 2626513221816)
        self.assertFalse(r["actual_derivatives_or_operator_constants_evaluated"])

    def test_review_scope(self):
        r = q207.review()
        self.assertTrue(all(r["checks"].values()))
        self.assertTrue(r["scope"]["fixed_stage_residual_constant_reconstructed"])
        self.assertTrue(r["scope"]["complete_pressure_viscosity_cutoffs_and_interactions_retained"])
        self.assertFalse(r["scope"]["stage_uniform_constant_required"])
        self.assertFalse(r["scope"]["primitive_same_object_continuum_bundle_proved_independently"])
        self.assertFalse(r["scope"]["unconditional_complete_residual_force_exclusion_or_blowup"])

    def test_invalid_identity_and_inputs(self):
        for call in (q207.statuses, q207.same_object_links, q207.rates,
                     q207.stage_for_power, q207.constant_construction,
                     q207.residual_ledger, q207.primitive_frontier,
                     q207.concrete_request, q207.review, q207.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError):
            q207.stage_for_power(m=-1)
        with self.assertRaises(ValueError):
            q207.stage_for_power(eta=0)

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 457)
        self.assertGreaterEqual(record["closure"]["edges"], 1661)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6142)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
