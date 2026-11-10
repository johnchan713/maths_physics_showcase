import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q156


class Q156Tests(unittest.TestCase):
    def test_pressure_force_and_gauge_ledger(self):
        row = q156.pressure_ledger()
        self.assertEqual(row["force_cancellations"], 1)
        self.assertIn("-div f", row["individual"])
        self.assertIn("H^-2", row["riesz_pressure_space"])
        self.assertIn("H^-3", row["averaged_gradient_space"])
        self.assertFalse(row["competitor_pressure_decay_assumed"])

    def test_commutator_exact_arithmetic(self):
        row = q156.commutator_ledger()
        self.assertEqual(row["near_radial_power"], Q(-2, 3))
        self.assertEqual(row["far_radial_power"], -2)
        self.assertEqual(row["radial_coefficient"], 4)
        self.assertEqual(row["norm_decay"], Q(3, 4))
        self.assertTrue(row["cutoff_sufficient"])
        self.assertFalse(q156.commutator_ledger(cutoff_power=7)["cutoff_sufficient"])

    def test_refined_radius_ledger(self):
        row = q156.radius_ledger()
        self.assertEqual(
            [item["young_decay"] for item in row["terms"]],
            [4, Q(14, 5)],
        )
        self.assertEqual(
            [item["sobolev_offset_decay"] for item in row["terms"]],
            [Q(5, 2), Q(5, 2)],
        )
        self.assertEqual(row["obsolete_unsplit_decay"], Q(4, 3))
        self.assertEqual(row["worst_decay"], 2)

    def test_comparison_has_no_extra_competitor_bounds(self):
        row = q156.comparison_ledger()
        self.assertIn("v=u", row["conclusion"])
        self.assertFalse(row["competitor_derivative_bound_assumed"])
        self.assertFalse(row["competitor_boundary_flux_bound_assumed"])
        self.assertFalse(row["competitor_dissipation_bound_assumed"])
        self.assertFalse(row["terminal_uniform_constants_assumed"])

    def test_core_gate_blocks_exclusion_claim(self):
        row = q156.core_gate()
        self.assertIn("e0,A,tau0>0", row["required"])
        self.assertFalse(row["Q155_frozen_interface_binds_quantitative_path"])
        self.assertFalse(row["global_smooth_competitor_exclusion_claimed"])
        self.assertFalse(row["blowup_claimed"])

    def test_provider_contracts(self):
        row = q156.provider_contracts()
        self.assertTrue(row["Q155_same_field_force_energy"])
        self.assertTrue(row["refined_comparison_lemma"])
        self.assertTrue(row["historical_Q10_frozen"])
        self.assertFalse(row["quantitative_same_physical_core_path"])

    def test_scope(self):
        row = q156.record()
        self.assertTrue(row["preterminal_comparison_certified_conditionally"])
        self.assertFalse(row["global_competitor_exclusion_or_blowup_certified"])
        self.assertFalse(row["historical_frozen_records_modified"])

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            q156.pressure_ledger("wrong")
        with self.assertRaises(ValueError):
            q156.commutator_ledger(1)
        with self.assertRaises(ValueError):
            q156.young_decay(2, 1)
        with self.assertRaises(ValueError):
            q156.radius_ledger("wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 390)
        self.assertGreater(record["closure"]["source_bindings"], 5300)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory() as td:
            target = Path(td) / "fresh.json"
            data = audit.generate()
            target.write_text(json.dumps(data, indent=2, sort_keys=True, default=str) + "\n")
            self.assertEqual(json.loads(target.read_text()), data)


if __name__ == "__main__":
    unittest.main()
