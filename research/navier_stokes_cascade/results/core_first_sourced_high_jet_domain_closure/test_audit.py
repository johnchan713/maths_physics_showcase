from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from sourced_join import normalized_join, physical_join, primitive_join, residual_obstruction, review


class FirstSourcedHighJetDomainClosureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()
        cls.record = audit.generate()

    def test_six_primitive_contracts_are_joined(self):
        row = primitive_join()
        self.assertEqual(len(row), 6)
        self.assertTrue(all(item["joined_on_current_certified_tail"] for item in row.values()))

    def test_interface_and_primitive_ceilings_are_distinct(self):
        row = normalized_join()
        self.assertEqual(row["first_real_interface_ceiling"], 300021)
        self.assertEqual(row["fully_expanded_primitive_ceiling"], 700030)
        self.assertGreater(row["fully_expanded_primitive_ceiling"], row["first_real_interface_ceiling"])

    def test_normalized_sources_and_strict_rates(self):
        row = normalized_join()
        self.assertEqual(row["A0_source_order"], 100005)
        self.assertEqual(row["A1_source_order"], 7)
        self.assertEqual(row["A0_tail_rate"], Q(3, 2))
        self.assertEqual(row["A1_tail_rate"], Q(3, 2))
        self.assertEqual(row["A0_margin"], row["kappa"])
        self.assertEqual(row["A1_margin"], row["kappa"])
        self.assertTrue(row["finite_and_acyclic"])

    def test_physical_request_is_finite_and_strict(self):
        row = physical_join()
        self.assertEqual(row["A0_source_order"], 60200007)
        self.assertEqual(row["A1_source_order"], 60100009)
        self.assertEqual(row["first_real_interface_ceiling"], 60400023)
        self.assertEqual(row["fully_expanded_primitive_ceiling"], 60800032)
        self.assertEqual(row["A0_margin"], Q(1, 20000000))
        self.assertEqual(row["A1_margin"], Q(1, 20000000))
        self.assertTrue(row["finite_graph"])
        self.assertFalse(row["uniform_in_all_orders_stages_or_scales"])

    def test_regular_residual_obstruction_is_retained(self):
        row = residual_obstruction()
        self.assertEqual(row["regular_order"], Q(3, 2) - Q(2, 100000))
        self.assertEqual(row["regular_physical_exponent"], Q(-15025001, 10000000))
        self.assertFalse(row["current_regular_decay_certified"])
        self.assertFalse(row["complete_residual_decay_proved"])

    def test_nonretained_kappa_is_rejected(self):
        for value in (0, Q(-1, 10), Q(1, 200000), Q(1, 10)):
            with self.assertRaises(ValueError):
                normalized_join(kappa=value)

    def test_controls_dependencies_and_scope(self):
        self.assertTrue(all(self.result["checks"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 10)
        self.assertTrue(all(self.record["gates"].values()))
        self.assertTrue(self.record["flags"]["actual_first_sourced_A0_and_A1_fixed_order_inequalities_proved_on_certified_tail"])
        self.assertFalse(self.record["flags"]["later_stage_compact_graph_recurrence_proved"])
        self.assertFalse(self.record["flags"]["complete_regular_or_physical_residual_decay_proved"])

    def test_scope_tampering_is_rejected(self):
        changed = deepcopy(self.record)
        changed["flags"]["complete_regular_or_physical_residual_decay_proved"] = True
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
