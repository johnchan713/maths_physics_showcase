#!/usr/bin/env python3
"""Tests for explicit same-cone tilt selection and pulse-tail arithmetic."""

from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import audit
from tail_bounds import (
    ceil_fraction, ceil_log2, ceil_root, physical_exponents,
    tail_certificate, tilt_certificate,
)


class ExplicitTiltPulseTailTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_exact_tilt_slack_identity(self):
        for n in range(1, 20):
            epsilon = F(1, 2**n)
            row = tilt_certificate(epsilon)
            self.assertEqual(row["u_squared"], 2 / epsilon)
            self.assertEqual(row["kappa"], F(2, 2 + epsilon))
            self.assertEqual(
                row["actual_deficit_floor"],
                epsilon * (1 + epsilon) / (2 + epsilon),
            )
            self.assertEqual(
                row["unused_slack"], epsilon * epsilon / (2 * (2 + epsilon))
            )

    def test_unspecified_admissible_tilt_need_not_pass_coarse_bound(self):
        epsilon = F(1, 2**40)
        ract = F(1, 4)
        old_u2 = F(1)
        self.assertGreater(old_u2 / (1 + old_u2), ract)
        self.assertGreater((1 + 1 / old_u2) * (1 - epsilon), 1)

    def test_exact_ceiling_helpers(self):
        self.assertEqual(ceil_fraction(F(7, 3)), 3)
        self.assertEqual(ceil_log2(F(1)), 0)
        self.assertEqual(ceil_log2(F(9, 8)), 1)
        self.assertEqual(ceil_log2(F(17)), 5)
        self.assertEqual(ceil_root(F(64), 6), 2)
        self.assertEqual(ceil_root(F(65), 6), 3)

    def fixture(self):
        return tail_certificate(
            F(1, 200), F(1000), F(1, 9), F(1, 8), F(13), F(17),
            F(19), F(6), F(1, 23), F(1, 64),
        )

    def test_tail_certificate_passes_complete_selected_gate(self):
        row = self.fixture()
        self.assertTrue(all(row["checks"].values()))
        self.assertGreaterEqual(row["H_lower"], row["gamma"] / 16)

    def test_every_later_band_preserves_gate(self):
        row = self.fixture()
        args = (F(1, 200), F(1000), F(1, 9), F(1, 8), F(13), F(17),
                F(19), F(6), F(1, 23), F(1, 64))
        base = row["ell_tail"]
        for extra in (0, 1, 2, 7, 31):
            ell = base + extra
            h, cost, rb, omega, lpsi, cmesh, cpert, dbox, dchart, epsilon = args
            tilt = tilt_certificate(epsilon)
            gamma = omega * omega * tilt["chosen_d0"] / (9 * tilt["kappa"])
            exponent = (2 * h * (ell - 2)).numerator // (2 * h * (ell - 2)).denominator
            error = cost / F(2**exponent)
            H = gamma / 8 - cmesh / F(ell**6)
            self.assertLessEqual(error, rb / 8)
            self.assertLessEqual(error, omega * omega * tilt["chosen_d0"] / (16 * lpsi))
            self.assertLessEqual(error, H / (8 * cpert))
            self.assertGreaterEqual(F(ell**6), dbox / dchart)

    def test_physical_exponent_ledger(self):
        row = physical_exponents(F(1, 200), 3, 2)
        self.assertEqual(row["sx"], F(201, 400))
        self.assertEqual(row["st"], F(403, 400))
        self.assertEqual(row["pressure"] - row["velocity"], row["sx"])
        self.assertEqual(row["velocity"] - row["potential"], row["sx"])

    def test_scope_not_promoted(self):
        for key in (
            "unspecified_historical_first_integer_verified",
            "old_localized_force_automatically_reidentified_with_new_tail",
            "all_generated_sources_complete_covariance_residual_or_blowup_verified",
        ):
            self.assertIs(self.record["flags"][key], False)

    def test_dependencies_bound(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 7)

    def test_frozen_record(self):
        audit.validate(json.loads((audit.HERE / "evidence.json").read_text()))

    def test_false_scope_promotions_rejected(self):
        for key, value in audit.FLAGS.items():
            if not value:
                changed = deepcopy(self.record)
                changed["flags"][key] = True
                with self.assertRaises(ValueError):
                    audit.validate(changed)

    def test_invalid_inputs(self):
        with self.assertRaises(TypeError):
            tilt_certificate(0.1)
        with self.assertRaises(ValueError):
            tilt_certificate(F(0))
        with self.assertRaises(ValueError):
            physical_exponents(F(1, 200), -1, 0)

    def test_fresh_replay_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory(dir=audit.PROJECT) as directory:
            output = Path(directory) / "fresh.json"
            audit.main([
                "--verify-record", str(audit.HERE / "evidence.json"),
                "--output", str(output),
            ])
            audit.validate(json.loads(output.read_text()))
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(output)])


if __name__ == "__main__":
    unittest.main()
