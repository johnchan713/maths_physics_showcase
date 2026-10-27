from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from localization_force import (
    TERM_LEDGER, endpoint, extension_request, force_jet_request, historical, review,
)


class LogReservedLocalizationForceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_exact_twelve_term_identity(self):
        self.assertEqual(tuple(name for name, _ in TERM_LEDGER), historical.TERMS)
        self.assertEqual(len(dict(TERM_LEDGER)), 12)
        self.assertEqual(dict(TERM_LEDGER)["pressure_cutoff"], "p0 grad c")

    def test_default_force_jet_orders(self):
        row = force_jet_request()
        self.assertEqual(row["force_orders"], {
            "potential": 3, "direct_swirl": 2, "pressure": 1,
        })
        self.assertEqual(row["Cauchy_orders"], {
            "potential": 4, "direct_swirl": 3, "pressure": 2,
        })
        self.assertEqual(row["residual_Cauchy_jet"], 1)
        self.assertTrue(all(row["checks"].values()))

    def test_arbitrary_fixed_force_requests_close(self):
        for a, j, power in ((1, 3, 5), (20, 4, 11), (6, 70, 2), (100, 100, 17)):
            row = force_jet_request(a, j, power)
            total = a + j
            self.assertEqual(row["Cauchy_orders"]["potential"], total + 4)
            self.assertEqual(row["Cauchy_orders"]["direct_swirl"], total + 3)
            self.assertEqual(row["Cauchy_orders"]["pressure"], total + 2)
            self.assertTrue(all(row["checks"].values()))

    def test_parent_request_is_used_without_new_field(self):
        self.assertEqual(force_jet_request(1, 3, 5)["comparison_stage"],
                         endpoint.request(1, 3, 5)["comparison_stage"])

    def test_future_force_diagonal(self):
        for a, m in ((0, 0), (1, 3), (20, 4), (6, 70), (100, 100)):
            row = extension_request(a, m)
            self.assertTrue(all(row["checks"].values()))
            self.assertGreater(row["first_controlled_index"], m)

    def test_exact_differentiated_factor(self):
        self.assertEqual(historical.derivative_factor(5, 2, 1), Q(1, 12))
        self.assertEqual(historical.derivative_factor(1, 3, 1), 0)

    def test_invalid_inputs_rejected(self):
        for value in (-1, True, 1.5):
            with self.assertRaises(ValueError):
                force_jet_request(value, 0, 0)
            with self.assertRaises(ValueError):
                extension_request(0, value)

    def test_review_scope(self):
        result = review()
        self.assertTrue(all(result["checks"].values()))
        scope = result["scope"]
        self.assertTrue(scope["same_field_axisymmetric_localization_replayed_conditionally"])
        self.assertTrue(scope["one_smooth_compact_force_extension_constructed_conditionally"])
        self.assertFalse(scope["preterminal_energy_replayed_here"])
        self.assertFalse(scope["pressure_recovery_comparison_or_blowup_proved"])

    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["preterminal_energy_replayed_here"] = True
        with self.assertRaises(ValueError):
            audit.validate(bad)

    def test_existing_output_not_overwritten(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
