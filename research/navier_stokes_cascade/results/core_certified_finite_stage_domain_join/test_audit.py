from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from certified_recurrence import (
    diagonal_obstruction, physical_join, pole_support_join,
    primitive_contracts, review, state_join,
)


class CertifiedFiniteStageDomainJoinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()
        cls.record = audit.generate()

    def test_primitive_contracts_are_bound(self):
        rows = primitive_contracts()
        self.assertEqual(len(rows), 6)
        self.assertTrue(all(row["joined_on_current_certified_tail"] for row in rows.values()))

    def test_state_interface_and_primitive_ceilings(self):
        zero, one = state_join(0), state_join(1)
        self.assertEqual((zero["first_real_interface_ceiling"], zero["fully_expanded_primitive_ceiling"]), (500038, 900047))
        self.assertEqual((one["first_real_interface_ceiling"], one["fully_expanded_primitive_ceiling"]), (1000067, 1400076))
        self.assertEqual(zero["provider_offset"], 400009)
        self.assertTrue(zero["finite_and_acyclic"] and one["finite_and_acyclic"])

    def test_physical_requests_are_expanded_and_strict(self):
        zero, one = physical_join(0), physical_join(1)
        self.assertEqual((zero["first_real_interface_ceiling"], zero["fully_expanded_primitive_ceiling"]), (60500042, 60900051))
        self.assertEqual((one["first_real_interface_ceiling"], one["fully_expanded_primitive_ceiling"]), (60900077, 61300086))
        self.assertEqual(min(zero["strict_margins"].values()), Q(1, 20000000))
        self.assertFalse(zero["uniform_in_stage_order_or_scale"])

    def test_common_pole_uses_primitive_ceiling(self):
        row = pole_support_join(2)
        self.assertEqual(row["requests"][0], ("physical", 61900109))
        self.assertEqual(row["archived_distinct_poles"], [Q(371400672)])
        self.assertEqual(len(row["operator_archive"]), 1)
        self.assertEqual(row["harmonics_at_most"], 16)
        self.assertTrue(row["same_closed_supports_preserved"])

    def test_diagonal_obstruction_is_exact(self):
        row = diagonal_obstruction()
        self.assertEqual(row["first_fixed_stage"], 1002)
        self.assertGreaterEqual(row["fixed_stage_power_margin"], 0)
        self.assertLess(row["preceding_power_margin"], 0)
        self.assertIn("log C_{j(q),m}", row["q_dependent_diagonal_required_inequality"])
        self.assertFalse(row["stage_constant_growth_bound_available"])
        self.assertFalse(row["one_field_all_targets_diagonal_proved"])

    def test_nonretained_kappa_is_rejected(self):
        for value in (0, Q(-1, 10), Q(1, 200000), Q(1, 10)):
            with self.assertRaises(ValueError):
                state_join(0, kappa=value)

    def test_controls_dependencies_and_scope(self):
        self.assertTrue(all(self.result["checks"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 10)
        self.assertTrue(all(self.record["gates"].values()))
        self.assertTrue(self.record["flags"]["every_fixed_stage_jet_and_power_has_a_finite_analytical_compact_graph"])
        self.assertFalse(self.record["flags"]["q_dependent_one_field_diagonal_proved"])
        self.assertFalse(self.record["flags"]["complete_forcing_energy_or_blowup_proved"])

    def test_scope_tampering_is_rejected(self):
        changed = deepcopy(self.record)
        changed["flags"]["q_dependent_one_field_diagonal_proved"] = True
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
