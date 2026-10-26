from copy import deepcopy
from fractions import Fraction as F
import json
import tempfile
import unittest
from pathlib import Path

import audit
from gap_join import (
    EDGES, NODES, accepted_tail_margins, background_floors, chart_tail,
    reference_normal_floor, topological_order,
)


class SameBandGapJoinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_route_is_acyclic(self):
        self.assertEqual(set(topological_order()), set(NODES))

    def test_cycle_is_rejected(self):
        with self.assertRaises(ValueError):
            topological_order(edges=EDGES + (("same_band_contract", "leading_regional_gaps"),))

    def test_background_floors(self):
        floors = background_floors(F(1, 5), 3, F(1, 7), F(1, 11), 4,
                                   F(1, 13), F(1, 17))
        self.assertTrue(all(value > 0 for value in floors.values()))
        with self.assertRaises(ValueError):
            background_floors(0, 3, F(1, 7), F(1, 11), 4, F(1, 13), F(1, 17))

    def test_one_integer_selector_rows(self):
        margins = accepted_tail_margins(
            400000, F(1, 200), 100, F(1, 6057072), 1000,
            F(1, 2 ** 2000), F(1, 2 ** 1000),
        )
        self.assertTrue(all(value > 0 for value in margins.values()))

    def test_reference_normal_floor(self):
        floor = reference_normal_floor(F(1, 32), F(1, 101))
        self.assertEqual(floor["amplitude_fraction"], F(31, 32))
        self.assertTrue(floor["sharper_than_coarse"])
        with self.assertRaises(ValueError):
            reference_normal_floor(F(1, 31), F(1, 101))

    def test_common_chart_tail(self):
        tail = chart_tail(F(1, 8), 20)
        self.assertTrue(tail["inside_qstar"])
        self.assertTrue(tail["at_or_after_selected_band"])

    def test_scope_and_correction(self):
        flags = self.record["flags"]
        self.assertTrue(flags["same_rule_defined_background_and_completed_leading_profile_used"])
        self.assertTrue(flags["reference_relative_full_phase_normal_floor_joined_on_that_tail"])
        self.assertFalse(flags["optional_native_Ctan_or_Cphi_margin_required_for_original_full_frame_route"])
        self.assertFalse(flags["unspecified_historical_first_band_or_old_localized_force_certified"])
        self.assertFalse(flags["canonical_homogeneous_full_path_ODE_independently_closed_here"])
        self.assertFalse(flags["unconditional_pulse_residual_force_or_blowup_proved"])

    def test_dependencies(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 8)
        self.assertTrue(self.record["inherited_source_checks"])
        self.assertTrue(all(self.record["inherited_source_checks"].values()))

    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        changed = deepcopy(frozen)
        changed["flags"]["canonical_homogeneous_full_path_ODE_independently_closed_here"] = True
        with self.assertRaises(ValueError):
            audit.validate(changed)

    def test_no_overwrite(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
