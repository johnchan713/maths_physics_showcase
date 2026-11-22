import json
import tempfile
import unittest
from pathlib import Path
import audit
import q211


class Q211Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q211.statuses().values()))

    def test_product_and_reciprocal_degrees(self):
        self.assertEqual(q211.product_degrees([0, 2, 4], [1, 3, 5]), [1, 3, 5])
        self.assertEqual(q211.reciprocal_degrees([0, 6, 8, 10]), [0, 6, 12, 18])

    def test_sourced_degree_recurrence(self):
        row = q211.affine_closure(100, 17)
        self.assertTrue(row["closes"])
        self.assertEqual(row["solution_degrees"], [7 * i + 18 for i in range(101)])
        self.assertFalse(row["actual_optimized_offset_evaluated"])

    def test_completed_route(self):
        row = q211.completed_route()
        self.assertTrue(row["Q210_next_step_already_present_conditionally"])
        self.assertFalse(row["these_records_are_unconditional_continuum_proofs"])

    def test_corrected_frontier(self):
        row = q211.corrected_frontier()
        self.assertFalse(row["source_grammar_can_create_superpolynomial_S_growth_from_polynomial_inputs"])
        self.assertIn("D_R^prim", row["primitive_degree_obligation"])
        self.assertFalse(row["unweighted_jets_moments_and_support_imply_primitive_degree_obligation"])
        self.assertTrue(row["Q199_countermodel_remains_applicable_to_that_invalid_inference"])
        self.assertFalse(row["independent_raw_map_certificate_complete"])

    def test_review_scope(self):
        row = q211.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["finite_degree_algebra_and_sourced_recurrence_proved"])
        self.assertFalse(row["scope"]["actual_primitive_degree_array_certified_independently"])
        self.assertFalse(row["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup"])

    def test_invalid_inputs(self):
        for call in (q211.statuses, q211.product_degrees, q211.reciprocal_degrees,
                     q211.sourced_solution_degrees, q211.affine_closure,
                     q211.completed_route, q211.corrected_frontier, q211.review, q211.record):
            with self.assertRaises((ValueError, TypeError)):
                if call in (q211.product_degrees, q211.sourced_solution_degrees):
                    call([0], [0], identity="wrong")
                elif call is q211.reciprocal_degrees:
                    call([0], identity="wrong")
                else:
                    call(identity="wrong")
        with self.assertRaises(ValueError):
            q211.product_degrees([0], [0, 1])
        with self.assertRaises(ValueError):
            q211.reciprocal_degrees([1, 2])
        with self.assertRaises(ValueError):
            q211.affine_closure(True)

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 460)
        self.assertGreaterEqual(record["closure"]["edges"], 1684)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6182)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
