import json
import tempfile
import unittest
from pathlib import Path
import audit
import q188


class Q188Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q188.statuses().values()))

    def test_exact_hierarchy_at_endpoint_and_larger_A(self):
        for A in (2**16, 2**17, 10**9):
            row = q188.exact_hierarchy(A)
            self.assertTrue(row["all_exact_reductions_hold"])
            self.assertTrue(all(row["checks"].values()))
            self.assertFalse(row["scales_materialized"])

    def test_every_acceptance_gate_is_retained(self):
        gates = q188.exact_hierarchy()["accepted_gates"]
        self.assertEqual(
            set(gates),
            {"state", "coefficient_radius", "contraction", "incoming_memory", "modulation_transfer"},
        )
        self.assertIn("G^64", gates["state"])
        self.assertIn("A^128", gates["coefficient_radius"])

    def test_full_pressure_memory(self):
        row = q188.pressure_memory_contract()
        self.assertEqual(row["forward_rows"], ["M", "I", "J", "S", "Cp"])
        self.assertEqual(row["pressure"], "Pi=Pi_axis+Cp")
        self.assertTrue(row["all_order_provider_keeps_five_rows"])
        self.assertTrue(row["regional_provider_uses_exact_global_ledger"])
        self.assertTrue(row["four_row_repair_is_insufficient"])

    def test_q187_provider_inventory_correction(self):
        row = q188.provider_correction()
        self.assertTrue(row["Q187_frontier_was_too_coarse_as_a_provider_inventory"])
        for key, value in row.items():
            if key != "independent_first_principles_reproof_completed_here":
                self.assertTrue(value, key)
        self.assertFalse(row["independent_first_principles_reproof_completed_here"])

    def test_precise_remaining_obstruction(self):
        row = q188.obstruction()
        self.assertIn("<=A", row["compact_envelope"])
        self.assertIn(">= A^-1", row["regional_gap"])
        self.assertFalse(row["frequency_or_repair_arithmetic_still_blocks"])
        self.assertFalse(row["tests_discharge_continuum_suprema"])

    def test_scope(self):
        row = q188.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["low_order_compact_repair_hierarchy_closed_relative_to_actual_provider_bounds"])
        self.assertTrue(row["scope"]["same_map_pressure_complete_cone_join_bound"])
        self.assertFalse(row["scope"]["all_provider_continuum_suprema_reproved_from_initial_equations"])
        self.assertFalse(row["scope"]["uniform_actual_pulse_or_complete_residual_proved"])

    def test_invalid_inputs(self):
        for call in (
            q188.statuses,
            q188.exact_hierarchy,
            q188.pressure_memory_contract,
            q188.provider_correction,
            q188.obstruction,
            q188.review,
            q188.record,
        ):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        for value in (True, 2**16 - 1, 1.0):
            with self.assertRaises(ValueError):
                q188.exact_hierarchy(value)

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 436)
        self.assertGreater(record["closure"]["edges"], 1500)
        self.assertGreater(record["closure"]["source_bindings"], 5800)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
