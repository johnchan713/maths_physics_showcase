from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from localization_force import (
    TERM_LEDGER, endpoint, extension_request, first_controlled,
    force_jet_request, historical, review,
)


class CurrentLogReservedLocalizationForceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()
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

    def test_current_parent_request_is_used_without_new_field(self):
        for a, j, power in ((1, 3, 5), (5, 4, 11), (3, 7, 2)):
            row = force_jet_request(a, j, power)
            parent = endpoint.endpoint_request(a, j, power)
            total = a + j
            self.assertEqual(row["comparison_stage"], parent["comparison_stage"])
            self.assertEqual(row["force_orders"]["potential"], total + 3)
            self.assertEqual(row["Cauchy_orders"]["pressure"], total + 2)
            self.assertTrue(all(row["checks"].values()))

    def test_high_fixed_order_ledger_is_algebraic(self):
        total = 200
        current = endpoint.endpoint.requirements(total)
        self.assertEqual(current["future_force_Cauchy_potential_bound"], 204)
        self.assertEqual(current["future_force_Cauchy_direct_swirl_bound"], 203)
        self.assertEqual(current["future_force_Cauchy_pressure_bound"], 202)

    def test_both_crosses_square_viscosity_and_pressure_retained(self):
        terms = dict(TERM_LEDGER)
        for name in ("cutoff_laplacian", "gradient_viscosity", "curl_viscosity",
                     "old_curl_transport", "curl_old_transport", "curl_square"):
            self.assertIn(name, terms)
        self.assertEqual(self.result["residual"]["viscosity"], 1)
        self.assertIn("-div f_minus", self.result["residual"]["pressure_Poisson"])

    def test_only_apparent_thirteenth_term_vanishes(self):
        self.assertEqual(
            self.result["localized_field"]["vanishing_term"],
            "(d.grad c)u0=0",
        )
        self.assertNotIn("d_grad_c", dict(TERM_LEDGER))

    def test_future_force_diagonal(self):
        for a, m in ((0, 0), (1, 3), (20, 4), (6, 70), (100, 100)):
            row = extension_request(a, m)
            self.assertTrue(all(row["checks"].values()))
            self.assertGreater(row["first_controlled_index"], m)
        self.assertEqual(first_controlled(0, 0), 1)

    def test_exact_differentiated_factor(self):
        self.assertEqual(historical.derivative_factor(5, 2, 1), Q(1, 12))
        self.assertEqual(historical.derivative_factor(1, 3, 1), 0)

    def test_invalid_inputs_are_rejected(self):
        for value in (-1, True, 1.5):
            with self.assertRaises(ValueError):
                force_jet_request(value, 0, 0)
            with self.assertRaises(ValueError):
                extension_request(0, value)

    def test_review_scope_and_controls(self):
        self.assertEqual(len(self.result["checks"]), 10)
        self.assertTrue(all(self.result["checks"].values()))
        scope = self.result["scope"]
        self.assertTrue(scope["current_same_field_axisymmetric_localization_replayed_conditionally"])
        self.assertTrue(scope["one_smooth_compact_force_extension_constructed_conditionally"])
        self.assertFalse(scope["preterminal_energy_replayed_here"])
        self.assertFalse(scope["pressure_recovery_comparison_or_blowup_proved"])

    def test_dependency_provenance_and_tampering(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 7)
        self.assertTrue(all(self.record["inherited_source_checks"].values()))
        self.assertTrue(all(self.record["gates"].values()))
        changed = deepcopy(self.record)
        changed["flags"]["preterminal_energy_or_comparison_replayed_here"] = True
        with self.assertRaises(ValueError):
            audit.validate(changed)

    def test_frozen_record_and_no_overwrite(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
