from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import audit
from corner_force import (
    cartesian_transfer_constant, comparison_stage, force_jet_requirements,
    logarithmic_half_gain_bound, review, support_distances,
)


class ScaleCompatibleCornerForceReplayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_comparison_stage_satisfies_and_predecessor_fails(self):
        h, delta = F(1, 200), F(1, 2) - F(3, 100000)
        row = comparison_stage(h, delta, 13, 17, 9)
        self.assertGreaterEqual(row["fixed_exponent"], 10)
        self.assertGreaterEqual(row["tail_exponent"], 10)
        j = row["stage"] - 1
        fixed = h * (F(3, 2) - F(2, 100000) + j * delta) - 13
        tail = h * delta * (j + 1) / 2 - 17
        self.assertTrue(fixed < 10 or tail < 10)

    def test_half_gain_bound_keeps_logarithm(self):
        self.assertEqual(logarithmic_half_gain_bound(F(7, 100), 5), F(2500, 7))
        with self.assertRaises(ValueError):
            logarithmic_half_gain_bound(0, 1)

    def test_cartesian_transfer_factor(self):
        self.assertEqual(cartesian_transfer_constant(0), 1)
        self.assertEqual(cartesian_transfer_constant(7), 4**7)
        with self.assertRaises(ValueError):
            cartesian_transfer_constant(-1)

    def test_normalized_support_does_not_create_physical_floor(self):
        values = support_distances([F(1, 2**n) for n in (1, 10, 100)])
        self.assertGreater(values[0], values[1])
        self.assertGreater(values[1], values[2])
        self.assertLess(values[2], F(1, 10**20))

    def test_force_jet_ledger(self):
        self.assertEqual(force_jet_requirements(4), {
            "potential": 7,
            "direct_swirl": 6,
            "pressure": 5,
            "next_time_derivative_reserved": True,
        })
        with self.assertRaises(ValueError):
            force_jet_requirements(-1)

    def test_review_keeps_core_obstruction(self):
        result = review()
        self.assertTrue(all(result["checks"].values()))
        self.assertTrue(result["analytic_conclusions"]["all_fixed_residual_jets_are_cartesian_corner_flat_relative_to_the_certified_provider_chain"])
        self.assertFalse(result["scope"]["same_field_positive_core_lower_bound_independently_proved_here"])
        self.assertFalse(result["scope"]["unconditional_blowup_or_competitor_exclusion_proved"])

    def test_dependencies_and_sources_are_frozen(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 8)
        self.assertTrue(all(self.record["inherited_source_checks"].values()))
        self.assertTrue(all(self.record["gates"].values()))

    def test_frozen_record_reproduces_and_tampering_fails(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["unconditional_blowup_or_competitor_exclusion_proved"] = True
        with self.assertRaises(ValueError):
            audit.validate(bad)

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
