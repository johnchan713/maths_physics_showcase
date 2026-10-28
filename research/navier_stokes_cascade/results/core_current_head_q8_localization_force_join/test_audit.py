from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from localization_force import (
    TERM_LEDGER, extension_request, first_controlled, force_jet_request,
    historical, review,
)


class CurrentHeadQ8LocalizationForceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()
        cls.record = audit.generate()

    def test_exact_twelve_term_identity(self):
        self.assertEqual(tuple(name for name, _ in TERM_LEDGER), historical.TERMS)
        self.assertEqual(len(dict(TERM_LEDGER)), 12)
        self.assertEqual(dict(TERM_LEDGER)["pressure_cutoff"], "p0 grad c")

    def test_sharp_direct_force_ledger(self):
        row = force_jet_request()
        self.assertEqual(row["sharp_direct_force_orders"], {
            "potential": 2, "direct_swirl": 1, "pressure": 0,
        })
        self.assertEqual(row["sharp_Cauchy_orders"], {
            "potential": 3, "direct_swirl": 2, "pressure": 1,
        })
        self.assertTrue(all(row["checks"].values()))

    def test_Q7_supplies_one_extra_field_derivative(self):
        for a, j, power in ((0, 0, 0), (1, 3, 5), (5, 4, 11), (3, 7, 2)):
            row = force_jet_request(a, j, power)
            total = a + j
            self.assertEqual(row["Q7_supplied_orders"], {
                "potential": total + 4,
                "direct_swirl": total + 3,
                "pressure": total + 2,
            })
            self.assertEqual(row["Q7_order_reserve"], {
                "potential": 1, "direct_swirl": 1, "pressure": 1,
            })
            self.assertEqual(row["residual_Cauchy_jet"], total + 1)

    def test_higher_order_changes_stage_not_endpoint_premise(self):
        low = force_jet_request(0, 0, 0)
        high = force_jet_request(5, 5, 10)
        self.assertGreater(high["comparison_stage"], low["comparison_stage"])
        self.assertEqual(
            low["endpoint_blocking_inequality"],
            high["endpoint_blocking_inequality"],
        )

    def test_unit_viscosity_crosses_square_and_pressure_retained(self):
        terms = dict(TERM_LEDGER)
        for name in (
            "cutoff_laplacian", "gradient_viscosity", "curl_viscosity",
            "old_curl_transport", "curl_old_transport", "curl_square",
        ):
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

    def test_scope_keeps_endpoint_premise_explicit(self):
        scope = self.result["scope"]
        self.assertTrue(scope["exact_Q1_Q7_field_localized_conditionally"])
        self.assertFalse(scope["Q7_endpoint_uniform_all_order_bound_proved_here"])
        self.assertFalse(scope["terminal_velocity_or_pressure_regular_at_origin_proved"])
        self.assertFalse(scope["preterminal_energy_replayed_here"])

    def test_dependency_closure_and_controls(self):
        self.assertEqual(len(self.result["checks"]), 12)
        self.assertTrue(all(self.result["checks"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 6)
        self.assertTrue(all(self.record["gates"].values()))

    def test_scope_tampering_is_rejected(self):
        changed = deepcopy(self.record)
        changed["review"]["scope"][
            "terminal_velocity_or_pressure_regular_at_origin_proved"
        ] = True
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
