import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q172

class Q172Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q172.statuses().values()))

    def test_q2_region_registry(self):
        self.assertEqual(len(q172.REGIONS), 10)
        self.assertEqual(q172.REGIONS[0][0], "activation")
        self.assertEqual(q172.REGIONS[-1][0], "intermediate_power")

    def test_orientation_identity(self):
        row = q172.orientation_identity(2, -1, 7, 3, 2)
        self.assertEqual(row["minus_T_dot_g"], row["factorized"])
        self.assertTrue(row["orientation_positive"])

    def test_orientation_rejects_bad_scale(self):
        for args in ((0, 1, 2, 3, 1), (1, 1, 2, 3, 0)):
            with self.assertRaises(ValueError): q172.orientation_identity(*args)

    def test_infinite_compact_union_counterexample(self):
        row = q172.compactness_counterexample()
        self.assertTrue(row["each_band_compact_and_finite"])
        self.assertTrue(row["prefix_supremum_increases"])
        self.assertFalse(row["compactness_alone_gives_uniform_infinite_tail_bound"])

    def test_q23_exact_tail_split(self):
        row = q172.q23_uniform_repair()
        self.assertTrue(row["tail_error_gate"])
        self.assertTrue(row["length_gate"])
        self.assertGreaterEqual(row["ell0"], 16)

    def test_q23_rejects_nonpositive_inputs(self):
        for kwargs in ({"h": 0}, {"r0": 0}, {"lreq": 0}):
            with self.assertRaises(ValueError): q172.q23_uniform_repair(**kwargs)

    def test_complete_pressure_terms(self):
        moving = q172.constrained_pressure((1,1,0),(1,-1,0),(1,-1,0),(0,0,0),3,1)
        source = q172.constrained_pressure((0,0,1),(0,0,0),(1,0,0),(0,0,2),5,-1)
        self.assertNotEqual(moving["moving_over_i"], 0)
        self.assertNotEqual(source["source_over_i"], 0)

    def test_pressure_rejects_nontangent_data(self):
        with self.assertRaises(ValueError):
            q172.constrained_pressure((1,0,0),(0,0,0),(1,0,0),(0,0,0))

    def test_scope_and_frontier(self):
        row = q172.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["Q2_provider_premise_source_audited"])
        self.assertFalse(row["scope"]["Q3_compactness_only_argument_valid"])
        self.assertTrue(row["scope"]["Q3_upper_propagator_repaired_by_Q23_on_same_field"])
        self.assertFalse(row["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"])

    def test_invalid_identity(self):
        for fn in (q172.statuses, q172.review):
            with self.assertRaises(ValueError): fn(identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreater(record["closure"]["nodes"], 410)
        self.assertGreater(record["closure"]["edges"], 1370)
        self.assertGreater(record["closure"]["source_bindings"], 5550)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as td:
            target = Path(td) / "fresh.json"
            target.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(target.read_text()), audit.generate())

if __name__ == "__main__": unittest.main()
