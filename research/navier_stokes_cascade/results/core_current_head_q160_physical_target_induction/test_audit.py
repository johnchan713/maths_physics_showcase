import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q160


class Q160Tests(unittest.TestCase):
    def test_physical_seed_matches_bridge_seed(self):
        r = q160.rates()
        self.assertEqual(r["B0"], Q(3, 2) - Q(2, 100000))
        self.assertEqual(r["C0"], 2 - Q(2, 100000))
        self.assertTrue(r["C_equals_B_plus_half"])
        self.assertTrue(q160.seed()["seed_closed"])
        self.assertFalse(q160.seed()["canonical_and_physical_entrances_identified"])

    def test_target_gain_is_stronger_than_required_cycle_gain(self):
        for j in (0, 1, 17, 10**4):
            r = q160.rates(j)
            self.assertEqual(r["target_slack"], Q(2, 100000))
            self.assertGreater(r["proved_target_order"], r["required_C_next"])

    def test_every_target_row_closes_on_whole_half_line(self):
        for C in (2 - Q(2, 100000), 3, 1000):
            r = q160.target_step(C)
            self.assertTrue(r["all_rows_close_bridge_target"])
            self.assertTrue(r["all_rows_close_physical_next_target"])
            self.assertTrue(r["whole_half_line"])

    def test_two_signed_crosses_are_the_sharp_rows(self):
        rows = q160.target_step()["rows"]
        self.assertEqual(rows["signed_old_error"]["margin_over_bridge_target"], 0)
        self.assertEqual(rows["signed_full_curl"]["margin_over_bridge_target"], 0)

    def test_weighted_closure_uses_products_and_zero_total_primitive(self):
        row = q160.weighted_closure()
        self.assertIn("sqrt(zeta)", row["covariance"])
        self.assertIn("zero R^e total", row["primitive"])
        self.assertIn("partial_Z sigma_e", row["moving_compensator"])
        self.assertFalse(row["uniform_in_stage_claimed"])

    def test_complete_moments_and_pressure_memory(self):
        row = q160.moments()
        self.assertIn("c_rho Pmean", row["axial"])
        self.assertIn("all six covariance", row["current_memory"])
        self.assertIn("Pi1/Pi2", row["current_memory"])
        self.assertTrue(row["zero_total_is_necessary"])
        self.assertTrue(row["preserved_at_every_finite_cycle"])

    def test_each_operation_preserves_moment_hypotheses(self):
        row = q160.moments()
        self.assertIn("nonzero angular", row["sourced_and_signed"])
        self.assertIn("exact zero Haar", row["fast"])
        self.assertIn("two zero velocity-moment", row["slow"])

    def test_provider_binding_preserves_historical_false_flag(self):
        row = q160.provider_contracts()
        self.assertTrue(row["same_physical_family_bound"])
        self.assertTrue(row["historical_source_flag_preserved_false"])
        self.assertTrue(row["new_Q160_theorem_supersedes_that_open_gate"])
        self.assertTrue(row["Q153_may_be_replayed_but_is_not_promoted_here"])

    def test_frontier_advances_only_P3(self):
        rows = q160.frontier()
        self.assertEqual([row["id"] for row in rows], [f"P{i}" for i in range(1, 8)])
        self.assertEqual(rows[2]["status"], "closed-by-Q160-induction")
        self.assertEqual(rows[3]["status"], "next-replay-required")
        self.assertEqual(rows[4]["status"], "downstream-conditional")

    def test_record_does_not_promote_downstream_claims(self):
        row = q160.record()
        self.assertTrue(row["P3_universal_finite_cycle_gate_discharged"])
        self.assertFalse(row["Q153_downstream_chain_replayed_without_P3_premise"])
        self.assertFalse(row["unconditional_complete_residual_force_exclusion_or_blowup_certified"])
        self.assertFalse(row["historical_frozen_records_modified"])

    def test_invalid_inputs_rejected(self):
        for j in (-1, Q(1, 2), "1"):
            with self.assertRaises(ValueError):
                q160.rates(j)
        with self.assertRaises(ValueError):
            q160.target_step(1)
        for fn in (q160.rates, q160.seed, q160.weighted_closure,
                   q160.moments, q160.provider_contracts, q160.frontier):
            with self.assertRaises(ValueError):
                if fn is q160.rates:
                    fn(identity="wrong")
                else:
                    fn("wrong")

    def test_frozen_evidence_and_recursive_closure(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 399)
        self.assertGreater(record["closure"]["edges"], 1290)
        self.assertGreater(record["closure"]["source_bindings"], 5420)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory() as td:
            target = Path(td) / "fresh.json"
            record = audit.generate()
            target.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(target.read_text()), record)


if __name__ == "__main__":
    unittest.main()
