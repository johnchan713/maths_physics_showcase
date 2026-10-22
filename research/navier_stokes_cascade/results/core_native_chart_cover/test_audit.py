#!/usr/bin/env python3
"""Existing-box domain, support-aware primitive and full-pressure scope guards."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from budgets import cover_existing,mesh_gate_compatibility,leading_global_rows,supported_reconstruct,join,fixture
from controls import review


class NativeCoverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate();cls.f=fixture()

    def control(self,part,fragment):
        rows={k:v for k,v in self.result['parts'][part]['checks'].items() if fragment in k}
        self.assertTrue(rows);self.assertTrue(all(rows.values()))

    def test_exact_physical_chart_witnesses_include_both_one_sided_poles(self):
        self.control('geometry','exact_physical_chart_equation')
        self.control('geometry','eta-1')
        self.control('geometry','eta1')
        self.assertTrue(all(r['physical_T_nonnegative'] for r in self.result['parts']['geometry']['rows']))

    def test_original_to_doubled_box_displacement_and_native_hull(self):
        self.control('geometry','within_selected_doubled_box_displacement')
        self.control('geometry','s_and_native_norm_hull_bounds')

    def test_native_log_distance_has_independent_rational_upper_witnesses(self):
        self.control('geometry','rational_native_log_distance_upper')

    def test_cartesian_norm_hull_is_not_a_growth_or_cone_domain(self):
        self.control('geometry','entire_norm_hull_includes_radii_outside_active_annulus')
        self.control('geometry','norm_hull_does_not_promote_growth_or_cone_gaps')

    def test_existing_selector_checks_are_retained(self):
        self.control('cover','all_original_selector_fixture_checks_still_pass')

    def test_same_enlargement_and_hull_work_on_all_later_bands(self):
        self.control('cover','same_all_later_cover_and_native_neighborhood')

    def test_existing_axis_clearance_failure_does_not_select_another_band(self):
        self.control('cover','existing_axis_clearance_failure_is_reported_not_reselected')

    def test_original_native_value_neighborhood_is_a_separate_gate(self):
        self.control('cover','existing_native_neighborhood_failure_is_separate_from_norm_cover')
        c=deepcopy(self.f['inputs']['cover'])
        c=cover_existing(c['existing_ell'],c['Xa'],c['Xb'],c['rA'],c['rB'],c['native_log_distance']/2)
        inputs=dict(self.f['inputs'],cover=c)
        out=join(**inputs)
        self.assertFalse(out['existing_value_neighborhood_passes'])
        self.assertFalse(out['infinite_Borel_schedule_and_actual_selector_verified'])

    def test_no_band_or_neighborhood_is_reselected(self):
        self.control('cover','supplied_band_and_neighborhood_are_never_changed')

    def test_original_mesh_gate_compatibility_is_explicit(self):
        self.control('cover','old_mesh_gate_can_imply_all_three_cover_inequalities')
        bad=mesh_gate_compatibility(1,Q(1,100),6,1)
        self.assertFalse(bad['implicit_s_barrier']);self.assertFalse(bad['existing_native_neighborhood'])

    def test_derived_enclosures_meet_all_parent_coverage_boundaries(self):
        self.control('join','derived_native_hull_exactly_meets_parent_coverage_guard')
        self.control('join','all_radius_enclosures_come_from_same_active_shell')

    def test_leading_velocity_rows_extend_past_the_old_finite_heat_interface(self):
        self.control('join','leading_F_rows_are_independent_of_finite_heat_interface_endpoint')
        self.control('join','axial_rows_keep_whole_support_primitive_inputs')

    def test_all_support_aware_raw_rows_are_dominated_by_original_hull_bound(self):
        self.control('support','all_support_aware_rows_dominated_by_old_whole_hull_bound')

    def test_radial_and_primitive_values_improve_with_original_axial_support(self):
        self.control('support','radial_and_primitive_value_norms_strictly_improve')
        self.control('support','only_axial_primitive_volume_uses_common_support')

    def test_full_hull_swirl_and_axial_bounds_are_not_restricted_to_axial_support(self):
        self.control('support','full_hull_is_retained_for_swirl_and_axial_jets')

    def test_chart_prefix_does_not_truncate_future_pressure_sources(self):
        self.control('support','chart_prefix_primitive_does_not_truncate_future_pressure_source')
        self.control('support','full_future_pressure_changes_the_first_moment_repair')

    def test_shortened_pressure_endpoint_is_rejected(self):
        inputs=deepcopy(self.f['inputs']);inputs['first_data']['x_velocity']/=2
        with self.assertRaisesRegex(ValueError,'SAME original full source support'):join(**inputs)

    def test_same_implemented_cutoff_and_complete_phase_budget_are_used(self):
        self.control('join','same_implemented_Borel_K2_is_used')
        self.control('join','full_phase_budget_remains_positive_and_conditional')

    def test_phase_bound_improves_without_shortening_pressure_memory(self):
        self.control('join','support_aware_norms_improve_manufactured_Cphi_without_shortening_pressure')

    def test_negative_T_extension_is_not_certified(self):
        self.assertFalse(self.f['inputs']['cover']['negative_T_extension_verified'])
        c=deepcopy(self.f['inputs']['cover']);c['negative_T_extension_verified']=True
        with self.assertRaises(ValueError):join(**dict(self.f['inputs'],cover=c))

    def test_geometry_cannot_be_overridden_or_given_a_different_core_anchor(self):
        inputs=deepcopy(self.f['inputs']);inputs['geometry']['Rmin']=1
        with self.assertRaises(ValueError):join(**inputs)
        inputs=deepcopy(self.f['inputs']);inputs['Xa']=Q(2)
        with self.assertRaisesRegex(ValueError,'SAME active annulus'):join(**inputs)

    def test_bad_roots_nonexact_inputs_and_unaccepted_geometry_floor_rejected(self):
        for args in ((15,1,64,1,8,Q(1,1000)),(16,1,64,-1,8,Q(1,1000)),(16,1,64,1,-8,Q(1,1000)),(16,1,64,2,8,Q(1,1000))):
            with self.assertRaises(ValueError):cover_existing(*args)
        with self.assertRaises(TypeError):cover_existing(16,1,64,1.0,8,Q(1,1000))

    def test_failed_existing_norm_cover_cannot_feed_phase_recipe(self):
        f=self.f;c=f['inputs']['cover'];delta=c['delta']
        bad=cover_existing(c['existing_ell'],delta*delta,64,delta,8,Q(1,1000))
        with self.assertRaisesRegex(ValueError,'norm coverage'):join(**dict(f['inputs'],cover=bad,Xa=delta*delta))

    def test_incomplete_hull_or_invalid_support_inputs_are_rejected(self):
        with self.assertRaises(ValueError):supported_reconstruct(Q(1,200),0,1,8,1,2,(1,1),(1,1),16)
        with self.assertRaises(ValueError):leading_global_rows(Q(1,200),2,1,800,7,128,13,101,303)

    def test_seven_frozen_records_and_declared_source_hashes_match(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),7)
        self.assertTrue(all(self.record['gates'].values()))

    def test_new_frozen_evidence_reproduces_exactly(self):
        self.assertEqual(json.loads(Path(__file__).with_name('evidence.json').read_text()),self.record)

    def test_actual_native_value_schedule_and_residual_scope_cannot_be_promoted(self):
        for key,value in self.record['flags'].items():
            if not value:
                changed=deepcopy(self.record);changed['flags'][key]=True
                with self.assertRaises(ValueError):audit.validate(changed)
        changed=deepcopy(self.record);changed['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(changed)

    def test_frozen_output_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')


if __name__=='__main__':unittest.main()
