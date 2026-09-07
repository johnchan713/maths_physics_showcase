"""Scientific gate tests: analytic decay, corrupted evidence, and trade-offs."""
import copy
import importlib.util
import math
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location(
    "robust_search", Path(__file__).resolve().parents[1] / "scripts" / "robust_search.py")
SEARCH = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SEARCH)


def shear_evidence():
    rows = []
    for resolution in ("coarse", "fine"):
        for i in range(9):
            t = i * 0.1 / 8
            amplitude = math.exp(-0.2 * t)
            values = {"resolution": resolution, "grid": 8 if resolution == "coarse" else 32,
                      "cutoff": 2, "sampling_grid": 32, "snapshot": i, "time": t,
                      "energy": 0.25 * amplitude**2, "enstrophy": amplitude**2,
                      "palinstrophy": 4 * amplitude**2, "h_half": amplitude,
                      "l3_sample": amplitude, "vorticity_sample": 2 * amplitude,
                      "vorticity_fourier_upper_bound": 2 * amplitude, "k_rms": 2,
                      "stretching": 0, "viscous_destruction": 0.4 * amplitude**2,
                      "net_enstrophy_rate": -0.4 * amplitude**2,
                      "production_to_dissipation": 0, "cutoff_fraction": 0}
            rows.append({k: str(v) for k, v in values.items()})
    return rows


def traces():
    return [{"resolution": grid, "stage": stage, "search_track": "amplification",
             "valid": "true", "peak_cutoff_fraction": "0.001",
             "steps": "100",
             "refinement_eligible": "false", "peak_vorticity_ratio": "1"}
            for grid, stage in (("coarse", "selected"), ("fine", "fine-validation"))]


class ScientificGates(unittest.TestCase):
    def evidence(self, rows=None):
        return SEARCH.summarize_evidence(rows or shear_evidence(), 0.1, 8, 32, 0.05, 0.25)

    def test_exact_decay_is_resolved_but_not_promising(self):
        result = SEARCH.assess(traces(), self.evidence())
        self.assertEqual(result["numerical_failures"], [])
        self.assertIn("h_final", result["physics_failures"])
        self.assertGreater(len(result["physics_failures"]), 2)
        self.assertAlmostEqual(result["conservative"]["h_final"], math.exp(-0.02))

    def test_missing_clock_corruption_and_nan_fail_closed(self):
        mutations = (("time", "0.33"), ("stretching", "99"), ("energy", "nan"),
                     ("sampling_grid", "16"), ("vorticity_sample", "200"))
        for key, value in mutations:
            rows = shear_evidence()
            rows[3][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.evidence(rows)
        with self.assertRaises(ValueError):
            self.evidence(shear_evidence()[:-1])

    def test_initial_field_mismatch_is_not_a_resolution_comparison(self):
        rows = shear_evidence()
        rows[9]["h_half"] = "1.01"
        with self.assertRaisesRegex(ValueError, "initial field mismatch"):
            self.evidence(rows)

    def test_profile_failure_does_not_block_amplification(self):
        evidence = self.evidence()
        for metrics in evidence.values():
            metrics.update(h_final=1.03, h_late=1.015, h_late_min_step=1.001,
                           scale_final=1.10, stretch_late_min=2.0)
        result = SEARCH.assess(traces(), evidence)
        self.assertFalse(result["profile_route_eligible"])
        self.assertEqual(result["physics_failures"], [])
        changed = copy.deepcopy(result)
        changed["evidence"]["fine"]["vorticity_final"] *= 1.03
        self.assertTrue(SEARCH.perturbation_failures(result, changed, "sampling"))

    def test_cutoff_margin_and_small_gain_do_not_pass(self):
        trace = traces()
        trace[0]["peak_cutoff_fraction"] = "0.0099"
        result = SEARCH.assess(trace, self.evidence())
        self.assertIn("cutoff margin", result["numerical_failures"])
        result["evidence"]["coarse"]["h_final"] = 1.004
        result["evidence"]["fine"]["h_final"] = 1.006
        self.assertFalse(SEARCH.gain_error_gate([result])["passed"])

    def test_pareto_preserves_tradeoffs(self):
        self.assertTrue(SEARCH.dominates((2, 3, -0.001), (1, 2, -0.002)))
        self.assertFalse(SEARCH.dominates((2, 3, -0.002), (1, 2, -0.001)))
        self.assertFalse(SEARCH.dominates((1.0000001, 2), (1, 2)))
        # Epsilon-relaxed weak comparisons can form a cycle among these points.
        cycle = ((0.00020, 0, 0.00010), (0.00010, 0.00020, 0), (0, 0.00010, 0.00020))
        for a, b in zip(cycle, cycle[1:] + cycle[:1]):
            self.assertFalse(SEARCH.dominates(a, b))
        base = SEARCH.assess(traces(), self.evidence())
        a = {"id": "a", "family": "wave-packets", "result": copy.deepcopy(base)}
        b = {"id": "b", "family": "vortex-tubes", "result": copy.deepcopy(base)}
        c = {"id": "c", "family": "orthogonal-bundle", "result": copy.deepcopy(base)}
        a["result"]["conservative"]["h_final"] = 1.08
        b["result"]["conservative"]["h_late"] = 1.07
        c["result"]["numerical_failures"] = ["unresolved"]
        frontier, selected = SEARCH.select_finalists([a, b, c], 2)
        self.assertEqual({r["id"] for r in frontier}, {"a", "b"})
        self.assertEqual({r["id"] for r in selected}, {"a", "b"})

    def test_unchanged_adaptive_steps_are_not_a_timestep_check(self):
        reference = SEARCH.assess(traces(), self.evidence())
        repeat = copy.deepcopy(reference)
        self.assertTrue(SEARCH.perturbation_failures(reference, repeat, "half-dt"))
        repeat["steps"] = {"coarse": 200, "fine": 200}
        self.assertEqual(SEARCH.perturbation_failures(reference, repeat, "half-dt"), [])

    def test_late_rate_clock_aggregation_and_negative_minimum(self):
        horizon, temperature, samples = 0.1, 1.0, 3
        rates = [-0.01, 0.8, 1.0]
        exponentials = [math.exp(-(rate - min(rates)) / temperature) for rate in rates]
        weights = [e / sum(exponentials) for e in exponentials]
        soft = min(rates) - temperature * math.log(sum(exponentials) / samples)
        self.assertGreater(soft, 0)
        trace, rows = traces(), []
        for row in trace:
            row.update(growth_objective="late-rate", horizon=str(horizon), late_window_start="0.5",
                       late_rate_samples=str(samples), late_rate_temperature=str(temperature),
                       late_rate_minimum=str(min(rates)), late_rate_maximum=str(max(rates)),
                       late_rate_soft_minimum=str(soft))
            for i, (rate, weight) in enumerate(zip(rates, weights)):
                rows.append(dict(resolution=row["resolution"], snapshot=str(i),
                                 time=str(horizon * (0.5 + 0.5 * i / (samples - 1))),
                                 critical_log_rate=str(rate), soft_minimum_weight=str(weight)))

        def assess_late(evidence_rows=rows, trace_rows=trace):
            result = SEARCH.assess(trace_rows, self.evidence())
            SEARCH.add_late_growth_assessment(result, evidence_rows, trace_rows, horizon, 0.5, samples, temperature)
            return result

        result = assess_late()
        self.assertEqual(result["numerical_failures"], [])
        self.assertIn("non-positive sampled late critical rate", result["physics_failures"])
        self.assertEqual(len(SEARCH.pareto_vector(result)), 7)
        for key, value in (("time", "0.03"), ("soft_minimum_weight", "0.99"),
                           ("critical_log_rate", "nan"), ("snapshot", "0.5")):
            corrupted = copy.deepcopy(rows)
            corrupted[1][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                assess_late(corrupted)
        with self.assertRaises(ValueError):
            assess_late(rows[:-1])
        bad_trace = copy.deepcopy(trace)
        bad_trace[0]["late_rate_soft_minimum"] = "4"
        with self.assertRaises(ValueError):
            assess_late(rows, bad_trace)
        repeated = copy.deepcopy(result)
        repeated["late_growth"]["evidence"]["fine"]["minimum"] += 0.1
        self.assertIn("sampling: fine late critical rate",
                      SEARCH.perturbation_failures(result, repeated, "sampling"))

    def test_late_rate_cross_resolution_gate_is_separate_from_norms(self):
        trace = traces()
        rows = []
        for row, rate in zip(trace, (0.2, 0.5)):
            row.update(growth_objective="late-rate", horizon="0.1", late_window_start="0.5",
                       late_rate_samples="3", late_rate_temperature="0.1",
                       late_rate_minimum=str(rate), late_rate_maximum=str(rate),
                       late_rate_soft_minimum=str(rate))
            for i in range(3):
                rows.append(dict(resolution=row["resolution"], snapshot=str(i), time=str(0.05 + 0.025 * i),
                                 critical_log_rate=str(rate), soft_minimum_weight=str(1 / 3)))
        result = SEARCH.assess(trace, self.evidence())
        SEARCH.add_late_growth_assessment(result, rows, trace, 0.1, 0.5, 3, 0.1)
        self.assertIn("cross-resolution late critical rate", result["numerical_failures"])


if __name__ == "__main__":
    unittest.main()
