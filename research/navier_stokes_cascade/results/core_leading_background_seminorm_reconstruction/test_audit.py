from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
import audit
from seminorms import derivative_slots, expected_slot_bound, finite_selector, locally_finite_support_count, maximum_surviving_slots, normalized_stage_slacks, review, word_majorant

class LeadingBackgroundSeminormReconstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.record = audit.generate()
    def test_surviving_derivative_slots(self):
        for k in range(13): self.assertEqual(maximum_surviving_slots(k), expected_slot_bound(k))
        with self.assertRaises(ValueError): derivative_slots((1, 1))
    def test_ordered_simplex_majorants_are_finite(self):
        rows = [word_majorant(k, F(7, 3), F(1, 5)) for k in range(16)]
        self.assertTrue(all(row > 0 for row in rows))
    def test_selector_uses_complete_fixed_row_list(self):
        with self.assertRaises(ValueError): finite_selector(F(1, 200), 3, (1, 2, 3), 4, 5)
        self.assertGreater(finite_selector(F(1, 200), 3, (4, 7, 9, 11), 13, 5), 0)
    def test_selector_has_strict_normalized_margin(self):
        b = finite_selector(F(1, 200), 3, (4, 7, 9, 11), 13, 5)
        row = normalized_stage_slacks(F(1, 200), 3, 11, b)
        self.assertGreaterEqual(row["strict_after_unit_padding"], 0)
    def test_local_finiteness_requires_separation(self):
        self.assertEqual(locally_finite_support_count(F(1, 20), (1, 2, 4, 8, 16, 32)), 5)
        with self.assertRaises(ValueError): locally_finite_support_count(F(1, 20), (1, 3, 4))
    def test_P1_is_split_at_exact_leading_leaf(self):
        result = review(); self.assertTrue(result["scope"]["coefficient_and_borel_part_of_P1_reconstructed"]); self.assertFalse(result["leading_certificate"]["independently_reproved_here"])
    def test_pressure_and_commutator_dependencies_are_bound(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 5); self.assertTrue(self.record["inherited_source_checks"]); self.assertTrue(all(self.record["inherited_source_checks"].values()))
    def test_scope_rejects_promotion(self):
        scope = review()["scope"]; self.assertFalse(scope["scale_selected_before_coefficient_seminorm"]); self.assertFalse(scope["numerical_infinite_scale_list_required"]); self.assertFalse(scope["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"])
    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text()); audit.validate(frozen); bad = deepcopy(frozen); bad["review"]["leading_certificate"]["independently_reproved_here"] = True
        with self.assertRaises(ValueError): audit.validate(bad)
    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"; path.write_text("keep")
            with self.assertRaises(SystemExit): audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")
if __name__ == "__main__": unittest.main()
