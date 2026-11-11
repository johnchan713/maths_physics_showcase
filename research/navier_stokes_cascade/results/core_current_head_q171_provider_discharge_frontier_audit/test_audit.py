import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q171


class Q171Tests(unittest.TestCase):
    def test_successor_chain_is_exact(self):
        row = q171.provider_contracts()
        self.assertTrue(row["all_identity_links_hold"])
        self.assertTrue(all(row["identity_chain"].values()))

    def test_historical_comparison_flags_are_preserved(self):
        row = q171.provider_contracts()
        self.assertFalse(row["historical_comparison_flag_all_dependencies_discharged"])
        self.assertTrue(row["historical_comparison_flag_core_and_finite_maps_still_dependencies"])
        self.assertTrue(row["Q164_later_closes_current_head_core_pressure_frontier"])
        self.assertTrue(row["Q165_later_closes_current_head_finite_map_frontier"])

    def test_successor_reconciliation_does_not_promote_providers(self):
        row = q171.successor_discharge()
        self.assertFalse(row["old_comparison_core_finite_map_flag_is_current_frontier"])
        self.assertFalse(row["later_audits_make_every_named_provider_unconditional"])
        self.assertEqual(row["protected_core_and_limit_order"], "later-conditional-audit-Q170")

    def test_exact_selector_sample(self):
        row = q171.selector_contract()
        self.assertEqual(row["h"], Q(1, 200))
        self.assertEqual(row["m_base"], Q(1, 4))
        self.assertEqual(row["ratios"], (Q(8), Q(8)))
        self.assertEqual(row["raised_power"], 100)
        self.assertEqual(row["ell"], 302)
        self.assertTrue(row["exact_selector_holds"])
        self.assertTrue(row["previous_integer_fails_when_ell_is_nontrivial"])

    def test_selector_rejects_missing_antecedents(self):
        bad = (
            {"margins": (0,), "constants": (1,)},
            {"margins": (-1,), "constants": (1,)},
            {"margins": (1,), "constants": (0,)},
            {"margins": (1,), "constants": (-1,)},
            {"margins": (1, 2), "constants": (1,)},
            {"margins": (), "constants": ()},
            {"margins": (1,), "constants": (float("inf"),)},
            {"margins": (float("nan"),), "constants": (1,)},
        )
        for kwargs in bad:
            with self.assertRaises(ValueError):
                q171.selector_contract(**kwargs)

    def test_selector_rejects_changed_h(self):
        with self.assertRaises(ValueError):
            q171.selector_contract(h=Q(1, 201))

    def test_first_open_premise_is_precise(self):
        row = q171.first_open_premise()
        self.assertIn("m_base=min", row["required_margin"])
        self.assertIn("log2(C_i/r_i)", row["blocking_selector"])
        self.assertTrue(row["Q2_provider_chain_bound_but_not_fully_reproved_in_join"])
        self.assertTrue(row["Q3_provider_chain_bound_but_not_fully_reproved_in_join"])
        self.assertFalse(row["smaller_q_alone_can_prove_missing_positivity_or_finiteness"])
        self.assertTrue(row["positive_finite_inputs_imply_a_finite_selector"])
        self.assertFalse(row["provider_claims_declared_false"])

    def test_frontier_moves_to_independent_continuum_audit(self):
        self.assertEqual(q171.frontier()[-1], {
            "id": "U9", "status": "next-independent-Q2-Q3-continuum-provider-audit",
        })

    def test_scope_rejects_unconditional_promotion(self):
        row = q171.record()
        self.assertFalse(row["historical_frozen_records_modified"])
        self.assertFalse(row[
            "unconditional_actual_pulse_residual_force_exclusion_or_blowup_certified"
        ])
        self.assertIn("do not independently reprove", row["obstruction"])

    def test_invalid_identity_rejected(self):
        for fn in (q171.provider_contracts, q171.successor_discharge,
                   q171.selector_contract, q171.first_open_premise,
                   q171.frontier):
            with self.assertRaises(ValueError):
                fn(identity="wrong")

    def test_frozen_evidence_and_recursive_closure(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 410)
        self.assertGreater(record["closure"]["edges"], 1370)
        self.assertGreater(record["closure"]["source_bindings"], 5546)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as td:
            target = Path(td) / "fresh.json"
            record = audit.generate()
            target.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(target.read_text()), record)


if __name__ == "__main__":
    unittest.main()
