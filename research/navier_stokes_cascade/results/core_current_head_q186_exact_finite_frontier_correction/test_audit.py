import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q186


class Q186Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q186.statuses().values()))

    def test_exact_core_frontier_is_stale(self):
        row = q186.exact_core_contract()
        self.assertTrue(row["all_index_core_and_pressure_audited"])
        self.assertTrue(row["strict_contraction"])
        self.assertTrue(row["complete_later_pressure_memory_retained"])

    def test_finite_root_and_map_frontier_is_stale(self):
        row = q186.finite_map_contract()
        self.assertTrue(row["Q2_Q5_chain_bound"])
        self.assertTrue(row["reference_angular_heat_maps_traced"])
        self.assertTrue(row["axial_map_and_strength_traced"])
        self.assertTrue(row["annular_and_repair_maps_traced"])
        self.assertTrue(row["annular_inverse_below_610"])
        self.assertTrue(row["repair_inverse_below_800"])

    def test_Q2_numeric_radius_is_fixture(self):
        row = q186.fixture_contract()
        self.assertEqual(row["radius"], Q(33, 387618781940))
        self.assertTrue(row["all_rows_positive"])
        self.assertTrue(row["manufactured_rational_fixture"])
        self.assertFalse(row["actual_selected_profile_constants"])
        self.assertFalse(row["continuum_infima_certified_by_fixture"])

    def test_actual_gate_is_not_fixture(self):
        row = q186.actual_radius_contract()
        self.assertEqual(
            row["blocking_inequality"],
            "Cval*q^(2h) <= rho_gap(actual selected-profile constants)",
        )
        self.assertTrue(row["regional_provider_derives_value_margins"])
        self.assertTrue(row["regional_provider_replays_analytic_lemmas"])
        self.assertFalse(row["actual_Cval_and_huge_constants_numerically_evaluated"])
        self.assertFalse(row["fixture_value_may_be_substituted_for_actual_rho_gap"])

    def test_review_scope(self):
        row = q186.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["current_head_frontier_corrected"])
        self.assertFalse(row["scope"]["actual_regional_infima_independently_rederived_here"])
        self.assertFalse(row["scope"]["fixture_promoted_to_actual_profile_certificate"])
        self.assertFalse(row["scope"]["complete_residual_or_blowup_proved"])

    def test_invalid_identity(self):
        for call in (
            q186.statuses,
            q186.exact_core_contract,
            q186.finite_map_contract,
            q186.fixture_contract,
            q186.actual_radius_contract,
            q186.review,
            q186.record,
        ):
            with self.assertRaises(ValueError):
                call(identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 431)
        self.assertGreater(record["closure"]["edges"], 1500)
        self.assertGreater(record["closure"]["source_bindings"], 5790)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
