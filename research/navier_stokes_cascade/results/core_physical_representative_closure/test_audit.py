#!/usr/bin/env python3
"""Physical representative membership, original hull and floor-only scope."""
from copy import deepcopy
from fractions import Fraction as Q
import inspect
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review,geometry
from budgets import verify_representative,norm_cover,floor_join


class PhysicalRepresentativeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()

    def control(self,fragment):
        rows=[v for k,v in self.r['checks'].items() if fragment in k]
        self.assertTrue(rows);self.assertTrue(all(rows))

    def test_both_physical_poles_and_native_shell_edges_are_covered(self):self.control('original_physical_active_clauses_and_box_pass')
    def test_original_representatives_are_not_moved(self):self.control('representative_not_moved')
    def test_nonconstant_chart_scales_are_checked_exactly(self):self.control('nonconstant_s_exact_physical_chart')
    def test_selected_boundary_box_may_have_negative_time_center(self):self.control('selected_boundary_box_can_have_negative_time_center')
    def test_negative_center_does_not_satisfy_original_representative_rule(self):self.control('negative_center_is_not_an_original_active_representative')
    def test_box_membership_alone_does_not_prove_native_active_membership(self):self.control('arbitrary_box_point_outside_native_annulus_is_rejected')
    def test_same_parent_norm_hull_applies_on_all_later_supplied_bands(self):self.control('same_original_norm_hull')
    def test_no_value_floor_extended_to_cartesian_hull(self):self.control('no_value_neighborhood_or_negative_time_certificate')
    def test_original_clearance_failure_is_reported_not_repaired(self):self.control('failed_original_axis_clearance_is_reported_without_reselection')
    def test_original_exact_root_clearance_boundary_matches(self):self.control('original_geometric_clearance_matches_exact_square_root_boundary')
    def test_same_complete_tangential_formula_is_used(self):self.control('same_complete_tangential_recipe')
    def test_unrelated_native_mu_is_not_a_floor_requirement(self):self.control('arbitrarily_narrow_unrelated_mu_can_fail_without_blocking_floor_recipe')
    def test_same_band_and_floor_scope_retained(self):self.control('original_band_and_all_scope_limits_retained')
    def test_supplied_normal_margin_success_and_failure_reported(self):
        self.control('supplied_old_small_scale_margin_failure_is_reported');self.control('supplied_large_manufactured_scale_passes_without_proving_actual_band')

    def test_no_mu_or_negative_time_extension_input_added(self):
        self.assertEqual(set(inspect.signature(norm_cover).parameters),{'ell','Xa','Xb','rA','rB'})
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])

    def test_primary_definition_and_all_active_set_clauses_retained(self):
        self.control('source_defines_intersection_representative');self.control('source_keeps_physical_native_and_original_q_domain_clauses')

    def test_relative_boundary_convention_and_global_extension_scope_separate(self):
        self.control('source_physical_relative_value_neighborhood');self.control('global_smooth_extension_remains_separate')

    def test_pdf_metadata_is_not_relabelled_as_fresh_binary_verification(self):self.control('historical_pdf_digest_metadata_matches_without_claiming_recomputed_binary_hash')
    def test_manuscript_theorems_are_not_assumed(self):self.control('manuscript_theorems_not_assumed')
    def test_invalid_or_float_representative_arguments_rejected(self):
        args=(Q(1,200),16,Q(1,2),2,(16**6,0,16**6),(1,0,1),1)
        with self.assertRaises(TypeError):verify_representative(*((0.005,)+args[1:]))
        with self.assertRaises(ValueError):verify_representative(*(args[:4]+((1,2),)+args[5:]))
        with self.assertRaises(ValueError):verify_representative(*(args[:6]+(0,)))
        r=verify_representative(Q(1,200),16,Q(1,2),2,(-16**6,0,16**6),(-1,0,1),1)
        self.assertFalse(r['gates']['positive_cylindrical_radius'])
        self.assertFalse(r['all_physical_native_clauses_pass'])

    def test_supplied_cover_tampering_or_failed_gate_cannot_feed_floor(self):
        c=norm_cover(16,1,4,1,2);bad=deepcopy(c);bad['negative_T_extension_certified']=True
        with self.assertRaises(ValueError):floor_join(bad,3,1,geometry())
        d=Q(1,16**6);bad=norm_cover(16,d*d,4,d,2)
        with self.assertRaisesRegex(ValueError,'coverage fails'):floor_join(bad,3,1,geometry())

    def test_geometry_or_growth_cannot_be_silently_overridden(self):
        c=norm_cover(16,1,4,1,2)
        for extra in ({'Dbox':4},{'Rmin':1}):
            with self.assertRaises(ValueError):floor_join(c,3,1,dict(geometry(),**extra))
        with self.assertRaises(ValueError):floor_join(c,3,4,geometry())
        with self.assertRaises(ValueError):norm_cover(15,1,4,1,2)
        with self.assertRaises(TypeError):norm_cover(16,1,4,1.0,2)

    def test_five_frozen_records_and_their_sources_match(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),5);self.assertTrue(all(self.record['gates'].values()))

    def test_frozen_new_evidence_reproduces_exactly(self):self.assertEqual(json.loads(Path(__file__).with_name('evidence.json').read_text()),self.record)
    def test_scope_or_remaining_obligations_cannot_be_promoted(self):
        for k,v in self.record['flags'].items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
        r=deepcopy(self.record);r['remaining_obligations'].pop()
        with self.assertRaises(ValueError):audit.validate(r)

    def test_frozen_output_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')


if __name__=='__main__':unittest.main()
