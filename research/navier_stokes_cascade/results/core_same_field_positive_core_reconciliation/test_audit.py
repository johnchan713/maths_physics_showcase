from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import audit
from core_path import (
    exact_core_threshold, leading_coefficient_square, localization_threshold,
    normalized_lower_bound, path_coordinates, review,
)


class SameFieldPositiveCoreReconciliationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_z_zero_path_coordinates(self):
        row = path_coordinates(F(1, 1000), F(1, 101))
        self.assertEqual(row["q"], F(1, 1000))
        self.assertEqual(row["eta"], 0)
        self.assertEqual(row["X"], F(1, 101))
        self.assertEqual(row["r_squared"], F(2, 101000))

    def test_leading_coefficient_has_strict_positive_square_bound(self):
        bound = leading_coefficient_square(F(1, 101), F(1, 10), F(33, 125), 10**6)
        self.assertGreater(bound, 0)
        with self.assertRaises(ValueError):
            leading_coefficient_square(0, 1, 1, 1)

    def test_exact_core_threshold_spends_half_the_leading_value(self):
        e0, constant = F(1, 10**10), F(10**5)
        qstar = exact_core_threshold(e0, constant, 100)
        base = e0 / (2 * constant)
        self.assertEqual(qstar, base**100)
        self.assertEqual(normalized_lower_bound(e0, constant, base), e0 / 2)

    def test_stricter_threshold_leaves_positive_slack(self):
        e0, constant = F(1, 100), F(7)
        base = e0 / (4 * constant)
        self.assertEqual(normalized_lower_bound(e0, constant, base), 3 * e0 / 4)

    def test_localization_plateau_uses_same_path(self):
        threshold = localization_threshold(F(1, 64), F(1, 16), F(1, 8))
        self.assertEqual(threshold, F(1, 128))
        with self.assertRaises(ValueError):
            localization_threshold(0, 1, 1)

    def test_review_corrects_parent_obstruction_without_overclaim(self):
        result = review()
        self.assertTrue(all(result["checks"].values()))
        self.assertTrue(result["analytic_conclusions"]["same_complete_diagonal_has_positive_core_lower_bound"])
        self.assertFalse(result["scope"]["annular_vanishing_incorrectly_applied_to_background_coefficients"])
        self.assertFalse(result["scope"]["unconditional_provider_chain_or_blowup_claimed"])

    def test_dependencies_and_sources_are_frozen(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 8)
        self.assertTrue(all(self.record["inherited_source_checks"].values()))
        self.assertTrue(all(self.record["gates"].values()))

    def test_frozen_record_reproduces_and_tampering_fails(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["unconditional_provider_chain_or_blowup_claimed"] = True
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
