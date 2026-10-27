from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import audit
from identities import (
    backward_pressure, core_density, forward_pressure, i4_background_radial_derivative,
    pressure_defect, review, split_identity, swirl_squared,
)


class SameFieldPressureReconciliationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_core_density_is_swirl_density(self):
        self.assertEqual(swirl_squared(7, 2, 3, 5) / 14, core_density(2, 3, 5))

    def test_forward_backward_identity_at_every_split(self):
        increments = (2, 3, 5, 7)
        axis = -sum(increments)
        for split in range(5):
            forward, backward, defect = split_identity(axis, increments, split)
            self.assertEqual(forward, backward)
            self.assertEqual(defect, 0)

    def test_axis_or_fifth_row_defect_is_detected(self):
        increments = (2, 3, 5)
        self.assertEqual(pressure_defect(-9, increments), 1)
        self.assertNotEqual(forward_pressure(-9, increments[:1]), backward_pressure(increments[1:]))

    def test_full_pressure_does_not_collapse(self):
        self.assertGreater(i4_background_radial_derivative(1, 3, 2), 0)

    def test_invalid_domains_are_rejected(self):
        with self.assertRaises(ValueError): core_density(1, 1, 0)
        with self.assertRaises(ValueError): swirl_squared(-1, 1, 1, 1)
        with self.assertRaises(ValueError): split_identity(-1, (1,), 2)
        with self.assertRaises(ValueError): i4_background_radial_derivative(1, 1, 0)

    def test_exact_controls_and_scope(self):
        result = review()
        self.assertTrue(all(result["checks"].values()))
        self.assertTrue(result["scope"]["leading_q_zero_pressure_identified_on_same_field"])
        self.assertFalse(result["scope"]["full_q_positive_pressure_identified_with_core_pressure"])
        self.assertFalse(result["scope"]["fixed_jet_closed_corner_residual_bound_proved_here"])

    def test_dependencies_and_sources_are_frozen(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 8)
        self.assertTrue(self.record["inherited_source_checks"])
        self.assertTrue(all(self.record["inherited_source_checks"].values()))

    def test_frozen_record_and_scope_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["full_q_positive_pressure_identified_with_core_pressure"] = True
        with self.assertRaises(ValueError): audit.validate(bad)

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit): audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
