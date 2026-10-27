from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from corner import (
    boundary_loss, corner_transfer_constant, first_stage, jacobian_floor,
    normalized_identity, pointwise_spike, review, stage_exponents,
)


class CornerFlatnessObstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_normalized_chart_identity_and_jacobian(self):
        self.assertEqual(normalized_identity(Q(3, 4), Q(1, 4)), 1)
        self.assertEqual(jacobian_floor(Q(1, 100)), Q(49, 50))

    def test_q_to_cartesian_corner_transfer(self):
        for power in (0, 1, 8, 32):
            self.assertEqual(corner_transfer_constant(power), 4**power)

    def test_closed_eta_boundary_loses_an_order(self):
        q = Q(1, 1024)
        for order in (0, 1, 5):
            self.assertEqual(boundary_loss(q, 1, order), q**order)

    def test_pointwise_flatness_is_not_compact_X_uniform(self):
        for n in (4, 8, 16):
            q = 2.0**(-n)
            self.assertGreater(pointwise_spike(q, q), 1/3)
        self.assertLess(pointwise_spike(1, 1/4), pointwise_spike(1/4, 1/4))

    def test_fixed_comparison_stage_exists(self):
        for m, target in ((0, 0), (0, 10), (3, 5)):
            stage = first_stage(m, target)
            exponents = stage_exponents(stage, m+2)
            self.assertGreaterEqual(exponents["fixed_stage_residual"], target+1)
            self.assertGreaterEqual(exponents["cutoff_tail"], target+1)

    def test_invalid_inputs_rejected(self):
        with self.assertRaises(ValueError): normalized_identity(-1, 2)
        with self.assertRaises(ValueError): jacobian_floor(1/10)
        with self.assertRaises(ValueError): corner_transfer_constant(-1)
        with self.assertRaises(ValueError): boundary_loss(0, 0, 1)
        with self.assertRaises(ValueError): pointwise_spike(-1, 1)
        with self.assertRaises(ValueError): first_stage(-1, 0)

    def test_controls_dependencies_and_scope(self):
        result = review()
        self.assertTrue(all(result["checks"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 7)
        self.assertTrue(all(self.record["inherited_source_checks"].values()))
        self.assertFalse(result["scope"]["closed_chart_constants_proved_finite_for_actual_field"])

    def test_frozen_record_and_scope_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["closed_chart_constants_proved_finite_for_actual_field"] = True
        with self.assertRaises(ValueError): audit.validate(bad)

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit): audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
