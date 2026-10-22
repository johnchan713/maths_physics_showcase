#!/usr/bin/env python3
"""Full signed-prefix source estimate, update order, support/premise boundary."""
from copy import deepcopy
from fractions import Fraction as Q
from math import factorial
import json
from pathlib import Path
import tempfile
import unittest
import audit
from budgets import P,full_mean,source_recipe,fixture,STATE,pole_ledger
from controls import review


class SignedPrefixSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate()

    def control(self,part,fragment):
        rows={k:v for k,v in self.result["parts"][part]["checks"].items() if fragment in k}
        self.assertTrue(rows);self.assertTrue(all(rows.values()))

    def test_full_local_conservative_and_advective_equations_match(self):
        for name in ("radial","angular","axial"):
            self.control("full_current_equations",f"full_{name}_conservative_and_advective_equations_match")

    def test_nonzero_base_and_mean_divergence_constraints(self):
        self.control("full_current_equations","full_local_base_and_mean_divergence_exact")

    def test_full_auxiliary_time_not_deleted_by_initial_Haar_argument(self):
        self.control("full_current_equations","pointwise_auxiliary_time_term_survives")

    def test_graph_laplacian_and_its_coefficient_derivative_retained(self):
        self.control("full_current_equations","full_graph_viscosity_cannot_be_replaced_by_Haar_laplacian")
        self.control("full_current_equations","graph_coefficient_derivative_survives_second_derivative")

    def test_radial_velocity_covariances_and_unaveraged_crosses_remain(self):
        self.control("full_current_equations","radial_mean_and_all_six_covariances_nonzero")
        self.control("full_current_equations","mean_cross_products_nonzero_before_any_Haar_average")

    def test_full_current_pressure_gradient_not_deleted(self):
        self.control("full_current_equations","current_pressure_gradient_not_zero_or_deleted")

    def test_every_retained_group_fits_coarse_mu_order(self):
        self.control("coarse_rates","meets_mu")

    def test_mean_phase_loss_and_no_improved_full_row_assumption(self):
        self.control("coarse_rates","mean_wave_phase_loss_is_essential")
        self.control("coarse_rates","no_improved_full_row_assumption_needed")

    def test_generic_current_row_polynomials_reproduce_exact_substitution(self):
        self.control("fresh_prefix_and_source","full_mean_positive_substitution")

    def test_pressure_time_graph_and_wave_covariance_costs_detected(self):
        for key in ("Kp","c","graph","w"):
            self.control("fresh_prefix_and_source",f"current_{key}_omission_detected")

    def test_source_recipe_reproduces_exact_substitution(self):
        self.control("fresh_prefix_and_source","source_recipe_exact_substitution")

    def test_pressure_complete_curl_and_actual_total_mean_enter_source(self):
        for key in ("Pi","R","M"):
            self.control("fresh_prefix_and_source",f"source_{key}_omission_detected")

    def test_pressure_recomputed_after_fast_before_slow(self):
        self.control("fresh_prefix_and_source","full_pressure_recomputed_after_fast")
        self.control("fresh_prefix_and_source","fresh_slow_uses_postfast_scalars")

    def test_both_maps_accumulate_old_means_and_actual_A1(self):
        self.control("fresh_prefix_and_source","total_means_accumulate_both_maps")
        self.control("fresh_prefix_and_source","actual_fast_A1_changes_final_source_budget")
        self.control("fresh_prefix_and_source","actual_fast:")

    def test_fast_source_factorials_are_separate_from_output_jet_count(self):
        b=fixture();m=0
        rawless=24*2**4*b.time_inverse(m)*(b.rows(4)["Etheta"]+(1+b.KA1(m))*b.rows(4)["Ez"])
        self.assertEqual(b.fast(0)["u"].record(),(factorial(4)*rawless).record())
        rawless_rad=b.KT1(1)*24*2**5*b.time_inverse(m)*b.rows(5)["Ez"]
        self.assertEqual(b.fast(0)["x"].record(),(factorial(5)*rawless_rad).record())
        self.assertGreater(b.fast(0)["u"].evaluate(1),rawless.evaluate(1))

    def test_complete_longitudinal_transport_negative_and_zero_outputs_replayed(self):
        self.control("fresh_prefix_and_source","complete_transport:")
        self.control("fresh_prefix_and_source","complete_linear_operator_exact_viscosity_pressure_connections")

    def test_product_poles_are_enlarged_before_weighted_inverse_constants(self):
        self.control("fresh_prefix_and_source","nonlinear_poles_enlarged_before_weighted_operator_calls")
        self.control("fresh_prefix_and_source","closed_hull_slow_response_and_final_source_poles")
        self.assertGreater(pole_ledger(2)["fresh_post_fast_rows_and_moments"],pole_ledger(2)["initial"])
        with self.assertRaises(ValueError):pole_ledger(-1)

    def test_first_source_uses_finite_initial_jets_and_no_later_state(self):
        self.control("fresh_prefix_and_source","source_zero_uses_finitely_many_initial_field_jets")
        self.control("fresh_prefix_and_source","source_does_not_request_later_complex_or_forced_states")
        for n in (1,2):
            trace=set();fixture(trace).source(n)
            self.assertEqual(max(j for name,j in trace if name=="initial"),n+15)

    def test_factorial_to_raw_conversion_is_explicit(self):
        b=fixture()
        for n in (0,1,2,3):self.assertEqual(b.raw_source(n).record(),(factorial(n)*b.source(n)["F"]).record())

    def test_new_source_callback_feeds_existing_forced_inverse_interface(self):
        from budgets import load,HERE
        inverse=load("new_source_inverse",HERE.parent/"core_first_sourced_inverse_budget/budgets.py")
        b=inverse.fixture();prefix=fixture();b.source=prefix.inverse_source(inverse)
        self.assertEqual(b.source(2).c,prefix.raw_source(2).c)
        self.assertGreater(b.velocity(1).evaluate(1),0)
        self.assertGreater(b.pressure(1).evaluate(1),0)

    def test_full_pressure_source_jet_order_is_available(self):
        b=fixture();n=2;N=n+3
        row=full_mean(n,dict(b.fixed(N),**b.initial(N),w=b.wave(N)["U"]*b.wave(N)["U"]))
        self.assertEqual(row["field_order"],n+3);self.assertEqual(row["pressure_source_order"],n+1)

    def test_exceptional_pressure_tails_not_frozen_or_removed(self):
        self.assertTrue(all(self.result["parts"]["complete_pressure_tails"]["checks"].values()))

    def test_invalid_orders_or_incomplete_or_mixed_bound_inputs_rejected(self):
        with self.assertRaises(ValueError):full_mean(-1,{})
        with self.assertRaises(ValueError):source_recipe(0,{})
        b=fixture();data=dict(b.fixed(3),**b.initial(3),w=P(1));data["x"]=1
        with self.assertRaises(TypeError):full_mean(0,data)
        data={k:Q(1) for k in data};data["x"]=-1
        with self.assertRaises(ValueError):full_mean(0,data)

    def test_six_frozen_dependencies_and_evidence_reproduce(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),6)
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_actual_lower_premises_full_residual_and_blowup_not_promoted(self):
        for key,value in self.record["flags"].items():
            if not value:
                changed=deepcopy(self.record);changed["flags"][key]=True
                with self.assertRaises(ValueError):audit.validate(changed)

    def test_omitted_remaining_obligation_or_historical_hash_rejected(self):
        for field in ("source_hashes","remaining_obligations"):
            changed=deepcopy(self.record)
            if isinstance(changed[field],dict):changed[field].pop(next(iter(changed[field])))
            else:changed[field].pop(0)
            with self.assertRaises(ValueError):audit.validate(changed)

    def test_existing_output_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"keep.json";p.write_text("sentinel")
            with self.assertRaises(SystemExit):audit.main(["--output",str(p)])
            self.assertEqual(p.read_text(),"sentinel")


if __name__=="__main__":unittest.main()
