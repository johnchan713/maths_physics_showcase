import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q194

class Q194Tests(unittest.TestCase):
    def test_provider_statuses(self): self.assertTrue(all(q194.statuses().values()))

    def test_base_scalings(self):
        row=q194.base_scalings(); self.assertTrue(row["velocity_matches"]); self.assertTrue(row["pressure_matches"]); self.assertTrue(row["potential_matches"])
        self.assertEqual(row["pressure"],-2*row["sx"]); self.assertEqual(row["potential"],0)

    def test_complete_curl_ledger(self):
        row=q194.curl_ledger(); self.assertEqual(row["component_count"],3); self.assertEqual(row["term_count"],7)
        self.assertTrue(row["axial_velocity_retained"]); self.assertTrue(row["cylindrical_connections_retained"]); self.assertTrue(row["cut_before_curl"])

    def test_every_small_allocation(self):
        for ns in range(8):
            for nt in range(6):
                for ps in range(ns+1):
                    for conn in range(ns-ps+1):
                        cs=ns-ps-conn
                        for pt in range(nt+1):
                            row=q194.allocation_budget(ns,nt,ps,cs,conn,pt,nt-pt)
                            self.assertTrue(row["fits"]); self.assertGreaterEqual(row["margin"],0)

    def test_field_losses(self):
        sx=q194.base_scalings()["sx"]
        for a in range(7):
            for b in range(5):
                v=q194.field_loss("velocity",a,b); p=q194.field_loss("pressure",a,b); A=q194.field_loss("potential",a,b)
                self.assertEqual(v["Q_loss"]-A["Q_loss"],sx); self.assertEqual(p["Q_loss"]-v["Q_loss"],sx)

    def test_support_contract(self):
        row=q194.support_contract(); self.assertTrue(row["cutoffs_on_potential_and_pressure_before_curl"]); self.assertTrue(row["derivatives_and_curl_do_not_enlarge_support"])
        self.assertTrue(row["primary_cutoff_supplies_smooth_pulse_endpoint_extension"]); self.assertTrue(row["generated_source_still_needs_own_endpoint_extension"])
        self.assertTrue(row["T_zero_jets_are_one_sided_not_zero_extended"])

    def test_scope_and_frontier(self):
        row=q194.review(); self.assertTrue(all(row["checks"].values())); self.assertTrue(row["scope"]["complete_primary_fixed_order_field_bounds_closed_conditionally"])
        self.assertFalse(row["scope"]["complete_wave_covariance_means_rederived"]); self.assertFalse(row["scope"]["complete_residual_improvement_or_blowup"])
        self.assertFalse(row["primary_assembly"]["positive_pointwise_lower_bound_for_complete_velocity"])

    def test_invalid_inputs(self):
        for call in (q194.statuses,q194.curl_ledger,q194.support_contract,q194.primary_assembly,q194.frontier,q194.review,q194.record):
            with self.assertRaises(ValueError): call(identity="wrong")
        with self.assertRaises(ValueError): q194.base_scalings(0)
        with self.assertRaises(ValueError): q194.field_loss("other",0,0)
        with self.assertRaises(ValueError): q194.allocation_budget(1,0,0,0,0,0,0)

    def test_frozen_evidence(self):
        record=audit.generate(); self.assertEqual(record,json.loads((audit.HERE/"evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"],444); self.assertGreater(record["closure"]["edges"],1580); self.assertGreater(record["closure"]["source_bindings"],6000)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path=Path(folder)/"fresh.json"; path.write_text(json.dumps(audit.generate(),indent=2,sort_keys=True)+"\n"); self.assertEqual(json.loads(path.read_text()),audit.generate())

if __name__=="__main__": unittest.main()
