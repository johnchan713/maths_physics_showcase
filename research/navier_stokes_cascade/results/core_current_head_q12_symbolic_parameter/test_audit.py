import json
import unittest
import audit
import parameter as p

Q = p.Q
recurrence = audit.load("q12_test_recurrence", "core_finite_stage_compact_recurrence/recurrence.py")
diagonal = audit.load("q12_test_diagonal", "core_complete_potential_diagonal_join/diagonal.py")
q5 = audit.load("q12_test_q5", "core_current_head_q5_weighted_inverse_source_join/weighted_inverse.py")
chart = audit.load("q12_test_chart", "core_corner_flatness_obstruction/corner.py")


class SymbolicParameterTests(unittest.TestCase):
    def test_exact_closed_recurrence_against_full_provider(self):
        for j in (0, 1, 2, 7, 27):
            for n in (0, 3, 91):
                self.assertEqual(p.state_ceiling(j, n),
                                 recurrence.state_ceiling(j, n)["first_real_ceiling"])

    def test_layer_coefficients_against_actual_call_graph(self):
        for i in (0, 1, 8, 10**30):
            row = recurrence.layer(i)
            self.assertEqual(row["sourced_integrations"], {"A0": 100000+49997*i, "A1": 4})
            self.assertEqual(row["generic_integrations"], {"A0": 149997+49997*i, "A1": 50001})
            self.assertEqual(row["sourced_backward_offset"], 200018+99994*i)
            self.assertEqual(row["generic_backward_offset"], 300011+99994*i)

    def test_enlarged_request_matches_recursive_provider(self):
        for h in (Q(1, 200), Q(1, 10**30)):
            for j, m, N in ((0, 0, 0), (1, 2, 3), (5, 7, 11)):
                old = recurrence.physical_request(j, m, N, h)
                new = p.compact_request(h, j, m, N)
                self.assertEqual(new["integrations_before_extra"], old["integrations"])
                self.assertEqual(new["base_strict_margins"], old["strict_margins"])
                self.assertEqual(new["raw_before_extra"], old["raw_sources"])
                # Every newly generated generic branch's maximum shifts by two
                # when BOTH integration counts increase by one.
                counts = old["integrations"]
                trace = recurrence.generic.dependency_trace(m, counts["A1"]+1, counts["A0"]+1)
                base = recurrence.generic.dependency_trace(m, counts["A1"], counts["A0"])
                self.assertEqual(trace["maxima"]["next_signed_target_input"],
                                 base["maxima"]["next_signed_target_input"]+2)
                self.assertEqual(new["enlarged_physical_primitive_ceiling"], old["first_real_ceiling"]+2)

    def test_symbolic_bounds_have_strict_target_margin(self):
        for m, N in ((0, 0), (1, 3), (5, 11), (10**30, 10**30)):
            row = p.universal_bounds(m, N)
            self.assertGreater(row["cross_lower"], N+Q(49, 25))
            self.assertGreater(row["regular_lower"], N+1)
            self.assertGreater(row["A1_rate_minus_threshold_times_h_lower"], 0)
            self.assertEqual(row["jacobian_lower"], Q(49, 50))

    def test_rational_fixtures_include_tiny_positive_parameters(self):
        for h in (Q(99, 10000), Q(1, 200), Q(1, 513), Q(1, 10**100)):
            for m, N in ((0, 0), (1, 3), (5, 11)):
                row = p.rational_plan(h, m, N)
                lower = p.universal_bounds(m, N)
                self.assertTrue(0 <= row["rounding_remainder"] < h)
                self.assertGreaterEqual(row["groups"]["old_new"], lower["cross_lower"])
                self.assertGreaterEqual(row["fixed_state_exponent"], lower["regular_lower"])
                self.assertGreater(row["complete_minimum"], N+1)
                self.assertEqual(row["request"]["integrations_before_extra"]["A1"], 1)
                self.assertGreater(row["request"]["base_strict_margins"]["A0"], 0)
                self.assertLessEqual(row["request"]["base_strict_margins"]["A0"], h*p.KAPPA)

    def test_sufficient_stage_dominates_existing_minimal_stage(self):
        for h in (Q(1, 200), Q(1, 10**100)):
            for m, N in ((0, 0), (1, 3), (5, 11)):
                row = p.rational_plan(h, m, N)
                old = diagonal.comparison(m, N, h)
                self.assertGreaterEqual(row["base_stage"], old["comparison_stage"])
                self.assertGreaterEqual(row["base_stage"], m+2)
                j = row["stage"]
                loss = diagonal.cut_loss(m+2, h)
                background = diagonal.total_loss(m+2, h)
                tail = diagonal.gain(j+1, h)/2-loss
                expected = dict(zip(p.GROUPS, (tail, tail, tail, tail-background, tail-background, 2*tail)))
                self.assertEqual(row["groups"], expected)
                regular = (h*diagonal.current.rates(j)["regular"]
                           - diagonal.residual_loss(m, h) - h*p.DELTA/2)
                self.assertEqual(row["fixed_state_exponent"], regular)
                self.assertEqual(row["complete_minimum"], min(row["fixed_state_exponent"], *expected.values()))

    def test_four_reserve_accounts_and_q5_absorption(self):
        h = Q(1, 10**30)
        row = p.rational_plan(h)
        for owner in ("A0", "A1"):
            account = q5.log_reserve(h, p.KAPPA, 37)
            self.assertEqual(row["reserves"][owner], account["retained"])
            self.assertEqual(account["constant_exponent"], Q(37**2)/(4*row["reserves"][owner]))
        self.assertEqual(set(row["reserves"]), {"A0", "A1", "regular", "cutoff"})
        self.assertEqual(row["reserves"]["cutoff"], h*p.DELTA*(row["stage"]+1)/2)

    def test_cutoff_majorant_for_all_s_above_threshold(self):
        for h in (Q(1, 200), Q(1, 10**100)):
            for j, c, P in ((1, 0, 0), (4, 1000, 300), (17, 10**30, 10**12)):
                gamma = h*p.DELTA*j/2
                required = p.cutoff_requirement(j, c, P, h)
                selected = diagonal.choose_log_scale(j, gamma, [(c, P)])
                self.assertGreater(selected, required)
                # This envelope decreases with s; checking the left endpoint
                # certifies the whole half-line, not a finite sample grid.
                self.assertLess(c + Q(P*P)/(2*gamma)-gamma*selected/2+j, 0)

    def test_cutoff_requirement_cannot_use_diagnostic_h(self):
        small = Q(1, 10**30)
        self.assertGreater(p.cutoff_requirement(1, 1, 1, small),
                           p.cutoff_requirement(1, 1, 1, Q(1, 200)))

    def test_chart_uses_explicit_parameter(self):
        for h in (Q(1, 513), Q(1, 10**100)):
            self.assertEqual(p.rational_plan(h)["jacobian_floor"], chart.jacobian_floor(h))
            self.assertGreater(chart.jacobian_floor(h), Q(49, 50))

    def test_invalid_domains_and_no_parameter_default(self):
        for h in (0, -1, Q(1, 100), 0.005, True):
            with self.assertRaises(ValueError):
                p.rational_plan(h)
        with self.assertRaises(TypeError):
            p.rational_plan()
        for n in (-1, True, Q(1, 2)):
            with self.assertRaises(ValueError):
                p.symbolic_root_plan(n, 0)

    def test_frozen_evidence_and_scopes(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE/"evidence.json").read_text()))
        self.assertTrue(all(value is False for value in record["scope"].values()))
        for row in record["symbolic_root_plans"]:
            self.assertEqual(row["h"], p.ROOT_H)
            self.assertFalse(row["actual_cutoff_condition_verified"])
            self.assertFalse(row["same_field_join_restored"])


if __name__ == "__main__":
    unittest.main()
