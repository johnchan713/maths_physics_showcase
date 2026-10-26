from copy import deepcopy
from fractions import Fraction as F
import json
import tempfile
import unittest
from pathlib import Path

import audit
from ode_join import (
    EDGES, NODES, ambient_lift, cutoff_review, projected_ode,
    topological_order, varying_growth_review,
)


class CanonicalFullPathODEJoinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_route_is_acyclic(self):
        self.assertEqual(set(topological_order()), set(NODES))

    def test_cycle_is_rejected(self):
        with self.assertRaises(ValueError):
            topological_order(edges=EDGES + (("full_path_contract", "same_background_tail"),))

    def test_projected_equation_and_moving_constraint(self):
        row = projected_ode(
            (2, 1, 1), (1, -1, 0),
            ((0, -2, 1), (3, 0, -1), (1, 2, 0)),
            (3, -2, 5), (1, -2, 0), 7, -2, F(3, 5),
        )
        self.assertEqual(row["constraint_derivative"], 0)
        self.assertTrue(all(value == 0 for value in row["original_residual"]))

    def test_pure_normal_source_pressure_obstruction(self):
        row = projected_ode(
            (0, 0, 1), (0, 0, 0), ((0, 0, 0),) * 3,
            (0, 0, 2), (1, 0, 0), 5, 1, 1,
        )
        self.assertTrue(all(value == 0 for value in row["projected_forcing"]))
        self.assertEqual(row["pressure_over_i"], F(2, 5))

    def test_ambient_lift_normal_identity(self):
        row = ambient_lift(
            (2, 1, 1), (1, -1, 0),
            ((0, -2, 1), (3, 0, -1), (1, 2, 0)), F(3, 5), F(4, 7),
        )
        self.assertEqual(row["normal_row"], row["expected_normal_row"])
        self.assertTrue(row["tangent_agreement"])

    def test_cutoff_derivative_terms_are_nonzero(self):
        row = cutoff_review()
        self.assertNotEqual(row["potential_derivative_commutator"], 0)
        self.assertNotEqual(row["pressure_derivative_commutator"], 0)

    def test_varying_growth_higher_jets(self):
        row = varying_growth_review()
        self.assertTrue(row["inverse_identity"])
        self.assertTrue(row["third_or_higher_nonzero"])
        self.assertTrue(row["frozen_rate_misses_jets"])

    def test_scope_and_correction(self):
        flags = self.record["flags"]
        self.assertTrue(flags["finite_interval_homogeneous_primary_existence_uniqueness_joined"])
        self.assertTrue(flags["adjacent_sourced_inverse_keeps_unprojected_normal_source_pressure"])
        self.assertFalse(flags["raw_outgoing_response_declared_zero_or_compactly_supported"])
        self.assertFalse(flags["canonical_and_inherited_primary_families_identified"])
        self.assertFalse(flags["fixed_compensator_projector_or_weighted_inverse_contracts_closed_here"])
        self.assertFalse(flags["unconditional_pulse_residual_force_or_blowup_proved"])

    def test_dependencies(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 8)
        self.assertTrue(self.record["inherited_source_checks"])
        self.assertTrue(all(self.record["inherited_source_checks"].values()))

    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        changed = deepcopy(frozen)
        changed["flags"]["raw_outgoing_response_declared_zero_or_compactly_supported"] = True
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
