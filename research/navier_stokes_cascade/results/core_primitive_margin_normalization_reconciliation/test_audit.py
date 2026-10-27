from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import audit
from margins import (
    NORMALIZED_SUPPORT_MARGIN, annular_contraction, background_floors,
    edge_weight_exponent, physical_support_margin, repair_contraction, review,
)


class PrimitiveMarginReconciliationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_background_certificate_has_five_positive_floors(self):
        floors = background_floors(F(1, 5), F(1, 7), F(1, 11), 4, F(1, 17))
        self.assertEqual(len(floors), 5)
        self.assertTrue(all(x > 0 for x in floors.values()))

    def test_support_margin_is_normalized_not_physical(self):
        self.assertEqual(NORMALIZED_SUPPORT_MARGIN, F(3, 8012))
        self.assertEqual(physical_support_margin(F(1, 8)), F(3, 8 * 8012))
        self.assertLess(physical_support_margin(F(1, 2 ** 100)), F(1, 10 ** 20))

    def test_annular_absolute_radius_may_shrink_but_slack_stays_positive(self):
        row = annular_contraction(F(1, 10 ** 16))
        smaller = annular_contraction(F(1, 10 ** 17))
        self.assertEqual(smaller["absolute_radius"], row["absolute_radius"] / 10)
        self.assertGreater(row["slack"], 0)

    def test_repair_uses_dimensionless_slack(self):
        row = repair_contraction()
        self.assertEqual(row["contraction"], F(1, 2 ** 219))
        self.assertEqual(row["normalized_uniqueness_margin"], 400)

    def test_edge_weight_rejects_positive_weight_floor(self):
        values = [edge_weight_exponent(F(1, 100), F(1, n), 1) for n in (100, 200, 400, 800)]
        self.assertEqual(values, sorted(values))

    def test_invalid_inputs_rejected(self):
        with self.assertRaises(ValueError): background_floors(0, 1, 1, 1, F(1, 2))
        with self.assertRaises(ValueError): physical_support_margin(0)
        with self.assertRaises(ValueError): annular_contraction(0)
        with self.assertRaises(ValueError): edge_weight_exponent(1, 0, 1)

    def test_controls_dependencies_and_scope(self):
        result = review()
        self.assertTrue(all(result["checks"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 8)
        self.assertTrue(all(self.record["inherited_source_checks"].values()))
        self.assertFalse(result["scope"]["unconditional_complete_residual_force_comparison_or_blowup"])

    def test_frozen_record_reproduces_and_tampering_fails(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["unconditional_complete_residual_force_comparison_or_blowup"] = True
        with self.assertRaises(ValueError): audit.validate(bad)

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit): audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
