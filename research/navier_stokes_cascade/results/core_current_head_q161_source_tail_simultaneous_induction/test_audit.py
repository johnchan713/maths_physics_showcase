import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q161


class Q161Tests(unittest.TestCase):
    def test_corrected_physical_rates(self):
        r = q161.rates()
        self.assertEqual(r["B0"], Q(3, 2) - Q(2, 100000))
        self.assertEqual(r["delta"], Q(1, 2) - Q(3, 100000))
        self.assertEqual(r["C_j"], r["B_j"] + Q(1, 2))
        self.assertEqual(r["target_slack"], Q(2, 100000))

    def test_coupled_step_closes_source_and_target(self):
        for j in (0, 1, 19, 1002):
            step = q161.coupled_step(j)
            self.assertTrue(step["all_regular_rows_reach_B_next"])
            self.assertTrue(step["closes_next_state"])
            self.assertEqual(step["target_slack"], Q(2, 100000))

    def test_full_linear_rows_are_source_gain_limiters(self):
        rows = q161.coupled_step()["source_rows"]
        self.assertEqual(rows["sourced_full_linear"]["margin_over_B_next"], 0)
        self.assertEqual(rows["signed_full_linear"]["margin_over_B_next"], 0)

    def test_source_quotient_bell_majorants_are_finite(self):
        for n in (0, 1, 4, 12):
            row = q161.source_quotient(n)
            self.assertTrue(row["all_finite"])
            self.assertEqual(len(row["quotient_majorant"]), n + 1)
            self.assertIn("source-normal", row["pressure_memory"])
            self.assertIn("no cutoff", row["endpoint"])

    def test_stage_zero_and_one_tail_requests_are_finite_and_strict(self):
        r0 = q161.tail_request()
        r1 = q161.tail_request(1)
        self.assertEqual(r0["raw_sources"], {"A0": 60200007, "A1": 60100009})
        self.assertEqual(r0["first_real_ceiling"], 60500042)
        self.assertEqual(r1["first_real_ceiling"], 60900077)
        self.assertTrue(r0["finite_request"] and r1["finite_request"])
        self.assertTrue(r0["strictly_positive"] and r1["strictly_positive"])

    def test_arbitrary_fixed_tail_request_stays_finite(self):
        row = q161.tail_request(7, 3, 11)
        self.assertTrue(row["finite_request"])
        self.assertTrue(row["strictly_positive"])
        self.assertGreater(row["first_real_ceiling"], 0)

    def test_residual_allocation_is_complete_and_double_count_free(self):
        row = q161.residual_allocation()
        self.assertTrue(row["double_count_free"])
        for required in ("pressure gradients", "viscosity", "both nonlinear crosses", "new square"):
            self.assertIn(required, row["regular"])
        self.assertIn("A0/A1", row["compact"])
        self.assertIn("never forward-inverted", row["Gaussian"])
        self.assertIn("only in the diagonal", row["scale_cutoffs"])

    def test_provider_contracts_close_P4_and_P5_only(self):
        row = q161.provider_contracts()
        self.assertTrue(row["P3_and_P4_closed_by_simultaneous_induction"])
        self.assertTrue(row["P5_weighted_inverse_and_fixed_stage_tail_contract_closed"])
        self.assertTrue(row["historical_conditional_records_preserved"])
        self.assertTrue(row["complete_potential_diagonal_provider_available_but_P6_not_replayed_here"])

    def test_frontier_advances_to_P6(self):
        rows = q161.frontier()
        self.assertEqual([row["id"] for row in rows], [f"P{i}" for i in range(1, 8)])
        self.assertEqual(rows[2]["status"], "closed-in-coupled-Q160-Q161-induction")
        self.assertEqual(rows[3]["status"], "closed-by-Q161")
        self.assertIn("fixed-stage", rows[4]["status"])
        self.assertEqual(rows[5]["status"], "next-replay-required")

    def test_scope_rejects_downstream_promotion(self):
        row = q161.record()
        self.assertTrue(row["P4_generated_source_endpoint_contract_discharged"])
        self.assertTrue(row["P5_fixed_stage_weighted_inverse_and_tail_contract_discharged"])
        self.assertFalse(row["P6_complete_diagonal_endpoint_force_chain_replayed"])
        self.assertFalse(row["unconditional_complete_residual_force_exclusion_or_blowup_certified"])
        self.assertFalse(row["historical_frozen_records_modified"])

    def test_invalid_inputs_rejected(self):
        for j in (-1, Q(1, 2), True):
            with self.assertRaises(ValueError):
                q161.rates(j)
        with self.assertRaises(ValueError):
            q161.source_quotient(-1)
        for fn in (q161.rates, q161.coupled_step, q161.source_quotient,
                   q161.tail_request, q161.residual_allocation,
                   q161.provider_contracts, q161.frontier):
            with self.assertRaises(ValueError):
                if fn in (q161.rates, q161.coupled_step):
                    fn(identity="wrong")
                elif fn is q161.source_quotient:
                    fn(identity="wrong")
                elif fn is q161.tail_request:
                    fn(identity="wrong")
                else:
                    fn("wrong")

    def test_frozen_evidence_and_recursive_closure(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 400)
        self.assertGreater(record["closure"]["edges"], 1300)
        self.assertGreater(record["closure"]["source_bindings"], 5435)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory() as td:
            target = Path(td) / "fresh.json"
            record = audit.generate()
            target.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(target.read_text()), record)


if __name__ == "__main__":
    unittest.main()
