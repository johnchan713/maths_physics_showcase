from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from inverse import (
    EDGES, NODES, fast_request_ledger, first_real_request_ledger,
    pressure_product_review, pure_power_margin_review, review,
    route_classification, slow_scale_powers, symbolic_two_edge,
    topological_order,
)


class CompactSlowInverseReconstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()
        cls.record = audit.generate()

    def test_route_is_acyclic_and_current_tail_precedes_inverse(self):
        order = topological_order()
        self.assertEqual(set(order), set(NODES))
        self.assertLess(order.index("current_fixed_profile_projectors"), order.index("native_haar_primitives"))
        self.assertLess(order.index("compact_slow_domain_contract"), order.index("log_renormalized_source_route"))

    def test_cycle_is_rejected(self):
        with self.assertRaises(ValueError):
            topological_order(edges=EDGES + (("log_renormalized_source_route", "current_fixed_profile_projectors"),))

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

    def test_pure_power_zero_margin_is_exact_but_route_specific(self):
        margin = pure_power_margin_review()
        route = route_classification()
        self.assertEqual(margin["source_margin"], 0)
        self.assertEqual(margin["viscosity_margin"], 0)
        self.assertFalse(margin["strict_without_extra_gain"])
        self.assertTrue(route["old_positive_majorant_pure_power_route_is_blocked"])
        self.assertFalse(route["zero_margin_is_the_current_global_frontier"])
        self.assertTrue(route["later_log_renormalized_source_route_is_the_current_downstream_route"])
        self.assertFalse(route["log_renormalized_source_induction_is_reproved_here"])

    def test_whole_moving_pressure_product_is_required(self):
        row = pressure_product_review()
        self.assertNotEqual(row["complete"], row["dropped_moving_compensator"])

    def test_all_controls_and_dependency_sources(self):
        self.assertTrue(all(self.result["checks"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 13)
        self.assertTrue(all(self.record["gates"].values()))

    def test_route_scope_tampering_is_rejected(self):
        changed = deepcopy(self.record)
        changed["flags"]["old_pure_power_zero_margin_is_current_global_frontier"] = True
        with self.assertRaises(ValueError):
            audit.validate(changed)

    def test_claim_tampering_is_rejected(self):
        changed = deepcopy(self.record)
        changed["flags"]["later_log_source_induction_or_residual_replay_reproved_here"] = True
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
