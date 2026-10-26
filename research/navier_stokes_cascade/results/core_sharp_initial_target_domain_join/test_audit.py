from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from rate_join import EDGES, NODES, covariance_rates, final_groups, rates, review, topological_order


class SharpInitialTargetDomainJoinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()
        cls.record = audit.generate()

    def test_route_is_acyclic(self):
        self.assertEqual(set(topological_order()), set(NODES))

    def test_cycle_is_rejected(self):
        with self.assertRaises(ValueError):
            topological_order(edges=EDGES + (("strict_first_real_final_row", "weighted_initial_haar_rows"),))

    def test_sharp_Haar_rate_is_not_downgraded(self):
        row = rates()
        self.assertEqual(row["beta"], 2 - 2 * row["kappa"])
        self.assertEqual(row["beta"] - row["alpha"], row["kappa"])
        self.assertEqual(row["projected_source"], row["beta"])

    def test_linear_signed_maps_preserve_order(self):
        row = rates()
        self.assertEqual(row["pointwise_response"], row["alpha"])

    def test_complete_covariance_rates(self):
        row = covariance_rates()
        k = Q(1, 100000)
        self.assertEqual(row["old_new_cross"], 3 - 4 * k)
        self.assertEqual(row["real_pair_remainder"], 3 - 4 * k)

    def test_radial_viscosity_has_strict_kappa_margin(self):
        row = rates()
        final = final_groups()
        self.assertEqual(final["fast_groups"]["radial_viscosity"], 3 - 5 * row["kappa"])
        self.assertEqual(final["strict_margin"], row["kappa"])
        self.assertGreater(final["strict_margin"], 0)

    def test_stale_witness_does_not_block_corrected_rate(self):
        row = rates()
        self.assertEqual(row["stale_D"] + 1 - 2 * row["kappa"], row["nu"])
        self.assertEqual(row["radial_viscosity"], row["nu"] + row["kappa"])

    def test_invalid_kappa_is_rejected(self):
        for value in (0, Q(1, 4), Q(-1, 10)):
            with self.assertRaises(ValueError):
                rates(value)

    def test_controls_dependencies_and_scope(self):
        self.assertTrue(all(self.result["checks"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 8)
        self.assertTrue(all(self.record["gates"].values()))
        self.assertFalse(self.record["flags"]["actual_same_state_high_jet_enclosures_proved"])
        self.assertFalse(self.record["flags"]["later_complex_wave_rate_improved"])

    def test_scope_tampering_is_rejected(self):
        changed = deepcopy(self.record)
        changed["flags"]["complete_residual_force_or_blowup_certified"] = True
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
