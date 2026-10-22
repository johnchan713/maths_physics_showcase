#!/usr/bin/env python3
"""Full primitive, derivative conventions, repaired inputs and domain guards."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from budgets import reconstruct,finalized_first,low_comparison,fixture,head
from controls import review


class RadialReconstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate();cls.f=fixture()

    def control(self,part,fragment):
        rows={k:v for k,v in self.result['parts'][part]['checks'].items() if fragment in k}
        self.assertTrue(rows);self.assertTrue(all(rows.values()))

    def test_all_four_components_raw_mixed_jets_bounded(self):
        for name in ('E','U','W','J'):self.control('reconstruction',name+'_all_raw_mixed_jets_bounded')

    def test_both_eta_endpoints_and_both_coefficient_orders_are_included(self):
        for s in ('_eta-1_','_eta1_','_n0_','_n1_'):
            self.control('reconstruction',s)

    def test_exact_divergence_derivative_avoids_extra_radial_order(self):
        self.control('reconstruction','exact_divergence_derivative')

    def test_incoming_radial_primitive_survives_a_zero_velocity_interval(self):
        self.control('memory','zero_local_U_does_not_erase_incoming_radial_memory')
        self.control('memory','leading_plateau_V_is_exactly_minus_mu')

    def test_coefficient_order_term_cannot_use_leading_plateau_formula(self):
        self.control('memory','positive_order_radial_memory_is_not_leading_plateau_formula')

    def test_plateau_pressure_source_is_not_reset(self):
        self.control('memory','leading_plateau_pressure_source_does_not_vanish')

    def test_U_first_angular_norm_alone_cannot_bound_radial_first_jet(self):
        self.control('memory','low_angular_U_norm_bounded_but_radial_first_jet_grows')

    def test_factorial_row_to_raw_second_derivative_cost_is_required(self):
        self.control('memory','factorial_row_is_not_identical_raw_second_jet')

    def test_full_pressure_first_repair_is_reused(self):
        self.control('join','same_full_pressure_Q03_Q13_enter_first_repair')
        self.control('join','removing_full_pressure_memory_changes_the_finalized_first_norm')

    def test_nonzero_axis_viscosity_boundary_is_preserved(self):
        self.control('join','nonzero_axis_viscosity_boundary_not_dropped')

    def test_first_field_value_matches_previous_K13(self):
        self.control('join','first_U_value_agrees_with_frozen_K13')

    def test_radial_coefficient_cutoff_and_fixed_repair_bump_jets_remain(self):
        self.control('join','first_mixed_rows_keep_coefficient_cutoff_derivative')
        self.control('join','repair_radial_bump_derivative_is_retained')

    def test_low_order_derivative_demands_are_precise(self):
        self.control('join','first_tuple_first_jet_uses_angular_order_two')
        self.control('join','leading_tuple_second_jet_uses_angular_order_three')

    def test_independent_first_order_closed_radial_and_primitive_formulas(self):
        self.control('join','independent_first_order_W_J_closed_form_equals_majorant')

    def test_Cphi_retains_every_parent_group(self):
        self.control('join','Cphi_has_all_seven_parent_groups')
        self.assertGreater(self.f['result']['comparison']['Cphi'],0)

    def test_derived_Cphi_does_not_silently_pass_the_original_selector(self):
        self.control('join','manufactured_derived_Cphi_does_not_automatically_certify_selector')

    def test_central_native_slice_does_not_cover_whole_band_chart(self):
        f=self.f
        with self.assertRaisesRegex(ValueError,'ENTIRE band chart'):
            low_comparison(f['h'],2,1,2,(7,11,13),(5,17,19),f['state'],3,f['geometry'])

    def test_lower_native_boundary_and_upper_chart_boundary_each_checked(self):
        f=self.f
        for Xhi,Rlo,Rhi in ((8,Q(3,4),4),(7,Q(1,2),4)):
            with self.assertRaisesRegex(ValueError,'ENTIRE band chart'):
                low_comparison(f['h'],Xhi,Rlo,Rhi,(7,11,13),(5,17,19),f['state'],3,f['geometry'])

    def test_exact_chart_coverage_boundaries_are_accepted(self):
        f=self.f
        out=low_comparison(f['h'],8,Q(1,2),4,(7,11,13),(5,17,19),f['state'],3,f['geometry'])
        self.assertEqual(out['comparison']['Cphi'],f['result']['comparison']['Cphi'])

    def test_same_inner_seed_and_complete_source_data_required(self):
        f=self.f;data=dict(f['first_data']);data['local_seed']+=1
        with self.assertRaises(ValueError):finalized_first(f['h'],f['inner'],f['cutoff'],data,Q(1,2),1)
        data=dict(f['first_data']);data.pop('nu')
        with self.assertRaises(ValueError):finalized_first(f['h'],f['inner'],f['cutoff'],data,Q(1,2),1)

    def test_Borel_cutoff_order_is_distinct_and_padded(self):
        f=self.f
        with self.assertRaisesRegex(ValueError,'Borel cutoff'):
            low_comparison(f['h'],8,Q(1,2),4,(7,11,13),(5,17,19),f['state'],Q(1,2),f['geometry'])

    def test_floats_negative_rows_or_missing_radial_rows_are_rejected(self):
        for args in ((.005,0,1,2,1,2,(1,1),(1,1)),(Q(1,200),0,1,2,1,2,(-1,1),(1,1)),(Q(1,200),0,1,2,1,2,(1,),(1,1))):
            with self.assertRaises((TypeError,ValueError)):reconstruct(*args)

    def test_axis_radius_and_unproved_order_cannot_be_used(self):
        with self.assertRaises(ValueError):reconstruct(Q(1,200),0,1,2,0,2,(1,1),(1,1))
        with self.assertRaises(ValueError):reconstruct(Q(1,200),0,3,2,1,2,(1,)*4,(1,)*4)
        with self.assertRaises(ValueError):reconstruct(Q(1,200),0,2,2,1,1,(1,)*3,(1,)*3)

    def test_seven_frozen_records_and_declared_sources_match(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),7)
        self.assertTrue(all(self.record['gates'].values()))

    def test_fresh_evidence_reproduces(self):
        self.assertEqual(json.loads(Path(__file__).with_name('evidence.json').read_text()),self.record)

    def test_scope_cannot_be_promoted_to_actual_selector_or_full_residual(self):
        for key,value in self.record['flags'].items():
            if not value:
                changed=deepcopy(self.record);changed['flags'][key]=True
                with self.assertRaises(ValueError):audit.validate(changed)

    def test_remaining_domain_and_Borel_obligations_cannot_be_removed(self):
        changed=deepcopy(self.record);changed['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(changed)
        changed=deepcopy(self.record);changed['source_hashes'].pop(next(iter(changed['source_hashes'])))
        with self.assertRaises(ValueError):audit.validate(changed)

    def test_frozen_evidence_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'keep.json';path.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(path)])
            self.assertEqual(path.read_text(),'sentinel')


if __name__=='__main__':unittest.main()
