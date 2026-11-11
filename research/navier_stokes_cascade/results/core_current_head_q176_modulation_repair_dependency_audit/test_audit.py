import json, tempfile, unittest
from fractions import Fraction as Q
from pathlib import Path
import audit, q176


class Q176Tests(unittest.TestCase):
    def test_statuses(self):
        self.assertTrue(all(q176.statuses().values()))

    def test_refined_dependency_order(self):
        order = q176.topological_order()
        pos = {node: i for i, node in enumerate(order)}
        self.assertEqual(set(order), set(q176.NODES))
        self.assertTrue(all(pos[a] < pos[b] for a, b in q176.EDGES))
        self.assertLess(pos["regional_value_selection"], pos["all_order_same_map"])
        self.assertLess(pos["all_order_same_map"], pos["matching_root_reaudit"])

    def test_dependency_order_guards(self):
        with self.assertRaises(ValueError):
            q176.topological_order(("a", "a"), ())
        with self.assertRaises(ValueError):
            q176.topological_order(("a",), (("a", "b"),))
        with self.assertRaises(ValueError):
            q176.topological_order(("a", "b"), (("a", "b"), ("b", "a")))

    def test_conservative_repair_gate(self):
        row = q176.repair_gate()
        self.assertEqual(row["jacobian_perturbation_upper_bound"], Q(400000000000, 8 ** 16))
        self.assertTrue(row["strict_half_gate"])
        self.assertEqual(row["inverse_bound"], 2000)

    def test_sharpened_repair_gate(self):
        row = q176.repair_gate(beta=800, q=25000)
        self.assertTrue(row["strict_half_gate"])
        self.assertLess(row["jacobian_perturbation_upper_bound"], q176.repair_gate()["jacobian_perturbation_upper_bound"])
        self.assertEqual(row["inverse_bound"], 1600)

    def test_repair_gate_guards(self):
        for args in ((0, 1, 8), (1, 0, 8), (1, 1, 0), (1.0, 1, 8)):
            with self.assertRaises(ValueError):
                q176.repair_gate(*args)

    def test_cone_gate(self):
        row = q176.cone_gate()
        self.assertEqual(row["left"], 512 * 4 ** 10)
        self.assertEqual(row["right_at_floor"], 4 ** 30)
        self.assertTrue(row["strict_cone_absorption"])
        with self.assertRaises(ValueError):
            q176.cone_gate(0)

    def test_compact_moment_inverse(self):
        rows = (Q(2), Q(-300), Q(3), Q(-7), Q(11))
        self.assertEqual(q176.compact_moment_inverse(rows, Q(1, 100)),
                         (Q(2), Q(3), Q(5), Q(7), Q(11)))
        self.assertEqual(q176.compact_moment_inverse((0, 0, 0, 0, 0), Q(1, 10)),
                         (Q(0),) * 5)
        with self.assertRaises(ValueError):
            q176.compact_moment_inverse((0,) * 4, 1)
        with self.assertRaises(ValueError):
            q176.compact_moment_inverse((0,) * 5, 0)

    def test_review_correction_and_scope(self):
        row = q176.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["correction"]["foundation_route_matching_roots_node_was_too_coarse"])
        self.assertFalse(row["correction"]["historical_frozen_evidence_modified"])
        self.assertTrue(row["established"]["conditional_same_map_all_order_modulation_and_exact_repair"])
        self.assertFalse(row["established"]["radial_field_errors_uniformly_small_at_all_derivative_orders"])
        self.assertFalse(row["scope"]["actual_selected_continuum_input_suprema_independently_recomputed_here"])
        self.assertFalse(row["scope"]["uniform_actual_pulse_or_complete_residual_proved"])
        self.assertFalse(row["scope"]["unconditional_force_exclusion_or_blowup_proved"])

    def test_invalid_identity(self):
        for fn in (q176.statuses, q176.review):
            with self.assertRaises(ValueError):
                fn(identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 421)
        self.assertGreater(record["closure"]["source_bindings"], 5637)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
