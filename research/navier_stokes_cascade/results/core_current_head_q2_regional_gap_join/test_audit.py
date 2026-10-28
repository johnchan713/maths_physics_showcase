from copy import deepcopy
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import tempfile
import unittest

import audit
from regional_gaps import DEPENDENCIES, REGIONS, load_closure, perturbation_radius, review


class CurrentHeadQ2RegionalGapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_provider_graph_and_all_hashes(self):
        records, graph, dependencies, sources = load_closure()
        row = review()
        self.assertGreater(len(records), 239)
        self.assertEqual(len(graph), row["graph"]["nodes"])
        self.assertTrue(row["checks"]["provider_graph_is_acyclic"])
        self.assertTrue(all(dependencies.values()))
        self.assertTrue(all(sources.values()))

    def test_all_direct_Q2_provider_statuses(self):
        row = review()
        self.assertEqual(len(DEPENDENCIES), 13)
        self.assertEqual(len(row["provider_status_checks"]), 12)
        self.assertTrue(all(row["provider_status_checks"].values()))

    def test_ten_region_cover_has_no_endpoint_gap(self):
        self.assertEqual(len(REGIONS), 10)
        self.assertEqual(REGIONS[0][1], "Xminus")
        self.assertEqual(REGIONS[-1][2], "Xrep=Xstar")
        for left, right in zip(REGIONS, REGIONS[1:]):
            self.assertEqual(left[2], right[1])

    def test_corrected_outer_ramp_floor(self):
        value = Decimal(review()["interface_controls"]["corrected_outer_ramp_floor"])
        self.assertGreater(value, Decimal("0.586"))
        self.assertGreater(value, Decimal("0.5"))

    def test_same_tilt_radius_is_exact_and_positive(self):
        row = perturbation_radius(
            Fraction(1, 5), 3, Fraction(1, 7), 4,
            Fraction(1, 11), Fraction(1, 13), Fraction(1, 17), 2,
        )
        self.assertEqual(row["radius"], Fraction(33, 387618781940))
        self.assertTrue(all(value > 0 for value in row["rows"].values()))

    def test_nonpositive_or_unordered_gap_is_rejected(self):
        with self.assertRaises(ValueError):
            perturbation_radius(0, 3, Fraction(1, 7), 4, Fraction(1, 11), Fraction(1, 13), Fraction(1, 17), 2)
        with self.assertRaises(ValueError):
            perturbation_radius(3, 2, Fraction(1, 7), 4, Fraction(1, 11), Fraction(1, 13), Fraction(1, 17), 2)

    def test_flat_edge_obstruction_and_valid_replacement(self):
        edge = review()["interface_controls"]["edge"]
        self.assertEqual(edge["absolute_open_annulus_stress_infimum"], "0")
        self.assertFalse(edge["positive_absolute_stress_floor"])
        self.assertEqual(edge["valid_lower_bound"], "|T0|>=c*zeta")
        self.assertTrue(edge["normalized_directions_extend_smoothly"])

    def test_uncertified_historical_tilt_route_is_excluded(self):
        controls = review()["interface_controls"]
        self.assertIn("explicit_tilt_pulse_tail", controls["tilt_route"])
        self.assertIn("leading_cone_margin", controls["excluded_route"])
        self.assertFalse(review()["scope"]["historical_uncertified_existing_tilt_gate_used"])

    def test_scope_advances_only_Q2(self):
        scope = review()["scope"]
        self.assertTrue(scope["Q1_exact_object_is_a_direct_dependency"])
        self.assertTrue(scope["Q2_regional_continuum_and_factored_edge_chain_bound_on_that_object"])
        self.assertFalse(scope["Q3_full_frame_primitive_domain_discharged"])
        self.assertFalse(scope["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"])

    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["Q3_full_frame_primitive_domain_discharged"] = True
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
