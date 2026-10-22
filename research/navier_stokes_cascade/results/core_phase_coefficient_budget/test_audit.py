#!/usr/bin/env python3
"""Actual-phase jet substitutions, varying growth, value-floor/premise guards."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from budgets import P,PhaseBudget,fixture,normal_gate,carrier_bounds,load,HERE,inverse
from controls import review


class PhaseCoefficientTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate()

    def control(self,part,fragment):
        rows={k:v for k,v in self.result["parts"][part]["checks"].items() if fragment in k}
        self.assertTrue(rows);self.assertTrue(all(rows.values()))

    def test_all_actual_phase_normal_and_moving_normal_components(self):
        self.control("actual_phase","all_three_actual_phase_normal_components_exact")
        self.control("actual_phase","moving_normal_is_path_derivative_not_full_normal_jet")
        self.control("actual_phase","normal_is_affine_in_path_with_frozen_carriers")

    def test_full_eikonal_retains_epsilon_before_estimation(self):
        self.control("actual_phase","full_eikonal_identity_has_exact_epsilon_factor")
        self.control("actual_phase","phase_cancellation_uses_full_F_and_G")

    def test_axial_transport_radial_flow_and_axial_normal_not_deleted(self):
        for name in ("full_axial_normal_error_not_dropped","axial_transport_term_in_eikonal_is_essential","nonzero_radial_background_term_in_eikonal_is_essential"):
            self.control("actual_phase",name)

    def test_logarithmic_growth_rate_is_not_frozen_lambda0(self):
        self.control("varying_growth_profile","growth_is_lambda0_over_sqrt_shape_not_constant")
        self.control("varying_growth_profile","normalization_center_has_zero_logarithmic_profile_rate")

    def test_profile_third_derivative_survives_quadratic_damping(self):
        self.control("varying_growth_profile","third_profile_derivative_is_not_zero")
        self.control("varying_growth_profile","damping_has_zero_third_derivative_but_growth_does_not")

    def test_all_requested_inverse_square_root_profile_jets_are_bounded(self):
        self.control("varying_growth_profile","actual_varying_growth_profile_jet_bounded")

    def test_actual_formula_coefficient_bounds_over_whole_manufactured_box(self):
        for name in ("whole_box_normal_bound","whole_box_matrix_bound","whole_box_eikonal_over_epsilon_bound"):
            self.control("full_box",name)

    def test_independent_raw_multibinomial_substitution(self):
        self.control("coefficient_polynomials","normal_component_substitution")
        self.control("coefficient_polynomials","moving_normal_has_no_extra_path_factor")

    def test_full_background_G_and_radial_flow_enter_eikonal_budget(self):
        self.control("coefficient_polynomials","full_G_transport_and_nonzero_radial_flow_affect_eikonal_budget")

    def test_frozen_nonzero_rounding_uses_unit_error_cost(self):
        self.control("coefficient_polynomials","nearest_nonzero_rounding_unit_cost_retained")
        self.assertEqual(carrier_bounds(1,0,1,2,1)["theta"],3)

    def test_five_variable_raw_to_factorial_conversion_is_explicit(self):
        self.control("coefficient_polynomials","raw_to_factorial_multinomial_count")
        self.assertGreater(fixture().normalized("normal",2,5).evaluate(1),fixture().normalized("normal",2,4).evaluate(1))

    def test_full_value_floor_uses_both_whole_path_inequalities(self):
        self.control("normal_value_gate","good_value_budget_meets_both_whole_path_gates")
        self.control("normal_value_gate","insufficient_normal_value_margin_is_rejected")
        self.control("normal_value_gate","large_accumulated_rounding_error_is_not_hidden")

    def test_no_component_floor_inferred_from_nonzero_angular_frequency(self):
        self.control("normal_value_gate","nonzero_angular_carrier_alone_is_not_a_uniform_normal_floor")
        self.control("normal_value_gate","floor_uses_full_or_tangential_reference_length")

    def test_rational_outward_growth_and_tilt_floor(self):
        self.control("normal_value_gate","squared_tilt_denominator_is_a_valid_outward_bound")
        gate=normal_gate(1,2,1,100,Q(1,20000))
        self.assertEqual(gate["beta"],Q(1,100));self.assertEqual(gate["inverse_normal_square"],400)

    def test_background_input_orders_and_ambient_lift_extra_jet(self):
        self.control("coefficient_polynomials","only_background_jets_through_n_plus_one")
        trace=set();b=fixture(trace);gate=normal_gate(1,2,1,100,Q(1,20000))
        out=b.forced(lambda j:P(1),P(2),P({1:2}),gate)
        out.velocity(0)
        self.assertEqual(max(j for name,j in trace if name in ("F","G","btilde")),2)

    def test_explicit_coefficients_feed_unchanged_forced_inverse(self):
        b=fixture();gate=normal_gate(1,2,1,100,Q(1,20000))
        out=b.forced(lambda j:P(j+1),P(2),P({1:2}),gate,harmonic=2)
        self.assertEqual(out.profile_rate(3).c,b.profile_rate(3).c)
        self.assertGreater(out.velocity(2).evaluate(1),0);self.assertGreater(out.pressure(2).evaluate(1),0)

    def test_failed_or_inconsistent_value_gate_cannot_feed_inverse(self):
        b=fixture()
        for gate in (normal_gate(1,2,1,10,Q(1,20000)),normal_gate(1,1,1,100,Q(1,20000)),normal_gate(4,2,1,1000,Q(1,2000000))):
            with self.assertRaises(ValueError):b.forced(lambda j:P(1),P(2),P({1:2}),gate)

    def test_coefficients_feed_actual_signed_prefix_source_interface(self):
        module=load("coefficient_to_prefix",HERE.parent/"core_signed_prefix_source_budget/budgets.py")
        prefix=module.fixture();b=fixture();fixed,linear=prefix.fixed,prefix.linear
        convert=lambda n:{k:module.P(v.c) for k,v in b.source_coefficients(n).items()}
        prefix.fixed=lambda n:dict(fixed(n),**{k:v for k,v in convert(n).items() if k in ("b","G","V","r")})
        prefix.linear=lambda n:dict(linear(n),**{k:v for k,v in convert(n).items() if k in ("normal","K","Eik","bR","Bz")})
        self.assertGreater(prefix.raw_source(1).evaluate(1),0)

    def test_large_higher_background_jets_change_constants_not_value_floor(self):
        b=fixture();large=fixture();old=large.F
        large.F=lambda n:1000000*old(n) if n>=5 else old(n)
        self.assertEqual(b.normal(0).record(),large.normal(0).record())
        self.assertGreater(large.normal(4).evaluate(1),b.normal(4).evaluate(1))

    def test_bad_bounds_types_orders_or_shell_inputs_rejected(self):
        with self.assertRaises(ValueError):normal_gate(0,1,1,100,0)
        with self.assertRaises(ValueError):carrier_bounds(1,1,0,1,1)
        b=fixture()
        with self.assertRaises(ValueError):b.normal(-1)
        with self.assertRaises(ValueError):b.normalized("normal",0,0)
        b.F=lambda n:1
        with self.assertRaises(TypeError):b.Hder(1)
        with self.assertRaises(ValueError):PhaseBudget(None,None,None,None,2,1,1,1,1,dict(theta=1,axial=1,radial=1))

    def test_five_frozen_dependencies_and_new_evidence_reproduce(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),5)
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_actual_core_gate_propagator_support_and_residual_not_promoted(self):
        for key,value in self.record["flags"].items():
            if not value:
                changed=deepcopy(self.record);changed["flags"][key]=True
                with self.assertRaises(ValueError):audit.validate(changed)

    def test_remaining_obligation_or_historical_hash_cannot_be_deleted(self):
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
