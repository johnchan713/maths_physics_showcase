from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from inverse_join import (
    EDGES, NODES, fast_request_ledger, first_real_request_ledger,
    pressure_product_review, review, slow_scale_powers,
    strict_margin_review, symbolic_two_edge, topological_order,
)


class CompactSlowDomainJoinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()
        cls.record = audit.generate()

    def test_route_is_acyclic(self):
        self.assertEqual(set(topological_order()), set(NODES))

    def test_cycle_is_rejected(self):
        with self.assertRaises(ValueError):
            topological_order(edges=EDGES + (("compact_slow_domain_contract", "fixed_profile_projectors"),))

    def test_symbolic_kernel_is_finite_without_exp_materialization(self):
        row = symbolic_two_edge(12, Q(1, 100), 8, 1)
        self.assertEqual(len(row), 3)
        self.assertTrue(all(isinstance(c, Q) and isinstance(e, Q) and c > 0 and e >= 0 for c, e in row))

    def test_invalid_kernel_domains_rejected(self):
        for args in ((-1, 1, 8, 1), (1, 0, 8, 1), (1, 1, 2, 2), (1, 1.0, 8, 1)):
            with self.assertRaises(ValueError):
                symbolic_two_edge(*args)

    def test_fast_request_channels_remain_distinct(self):
        row = fast_request_ledger(3)
        self.assertEqual(row["temporal_tangential_source"], 7)
        self.assertEqual(row["radial_potential_source"], 8)
        self.assertEqual(row["whole_axial_pressure_moment"], 12)

    def test_slow_physical_scale_powers(self):
        row = slow_scale_powers()
        self.assertEqual(row["u_from_P"], Q(101, 200))
        self.assertEqual(row["u_from_Jz"], Q(-99, 200))
        self.assertEqual(row["gamma"], Q(-199, 200))
        self.assertEqual(row["psi"], row["u_from_Jz"])

    def test_first_real_full_request_graph(self):
        row = first_real_request_ledger()
        self.assertEqual((row["p1"], row["p0"]), (99999, 199996))
        self.assertEqual(row["A1_raw_source"], 200001)
        self.assertEqual(row["A0_raw_source"], 399995)
        self.assertEqual(row["max_signed_target_input"], 400009)
        self.assertEqual(row["common_pole"], 2400066)

    def test_zero_margin_obstruction_is_exact(self):
        row = strict_margin_review()
        self.assertEqual(row["source_margin"], 0)
        self.assertEqual(row["viscosity_margin"], 0)
        self.assertFalse(row["strict_without_extra_gain"])
        self.assertEqual(row["required_inequality"], "delta_source + delta_viscosity > 0")

    def test_whole_moving_pressure_product_is_required(self):
        row = pressure_product_review()
        self.assertNotEqual(row["complete"], row["dropped_moving_compensator"])

    def test_all_controls_and_dependency_sources(self):
        self.assertTrue(all(self.result["checks"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 10)
        self.assertTrue(all(self.record["gates"].values()))

    def test_scope_tampering_is_rejected(self):
        changed = deepcopy(self.record)
        changed["flags"]["strict_final_row_gain_from_positive_majorants_proved"] = True
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
