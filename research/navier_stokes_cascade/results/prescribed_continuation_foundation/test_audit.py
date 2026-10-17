#!/usr/bin/env python3
"""Continuation definitions, exact flows, pressure memory, seminorms and honest scope."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from math import factorial
from pathlib import Path
import tempfile
import unittest
import mpmath as mp
from continuation_bounds import step_cost,primitive_bound,product_convolution,bounds_review,join_review
from continuation_checks import reference,prescribed,exponential,primitive,zero_through,jet,stats,modules
from audit import HERE,FLAGS,generate,validate,compare_records,parent_review,main


class ContinuationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_current_parent_and_unchanged_width(self):
        self.assertEqual(self.record['protocol']['parent_commit'],'bde6e3d9d212a706b54d8688941379400f60f79d')
        self.assertEqual(self.record['protocol']['fixed_width'],'t1=kappa0=C^-20')

    def test_pre_modulation_domain_not_final_actual_profile(self):
        self.assertIn('pre-modulation',self.record['protocol']['domain'])
        self.assertFalse(FLAGS['pre_modulation_profile_equals_final_modulated_profile_on_whole_region'])

    def test_common_nested_strip_and_no_derivative_feedback(self):
        self.assertTrue(all(bounds_review()['checks'].values()))
        self.assertEqual(bounds_review()['final_strip'],'delta=parent_delta/16=rho/(12800*S)')

    def test_factorial_product_bound_on_independent_positive_derivative_lists(self):
        for n in range(1,12):
            a=[Q(i+1,3) for i in range(n)];b=[Q(i+2,7) for i in range(n)]
            c=product_convolution(a,b)
            self.assertLessEqual(sum(v/factorial(i) for i,v in enumerate(c)),
                sum(v/factorial(i) for i,v in enumerate(a))*sum(v/factorial(i) for i,v in enumerate(b)))

    def test_Volterra_bound_keeps_original_trace_and_long_interval(self):
        self.assertEqual(primitive_bound(7,3,1000),3010)
        self.assertEqual(primitive_bound(7,3,0),10)

    def test_narrow_step_width_cost_not_replaced_by_value_bound(self):
        self.assertEqual(step_cost([1,1,1,1],Q(1,100)),Q(515303,3))
        self.assertGreater(step_cost([1]*8,Q(1,10000)),10**24)

    def test_invalid_norm_inputs_rejected(self):
        for t in (0,-1):
            with self.assertRaises(ValueError):step_cost([1],t)
        with self.assertRaises(ValueError):primitive_bound(-1,1,1)
        with self.assertRaises(ValueError):product_convolution([],[])

    def test_exact_exponential_and_primitive_Taylor_operations(self):
        y=modules()[0].Jet.variable(0,0,n=2,order=5)
        z=exponential(y)
        for i in range(6):self.assertEqual(z.partial((i,0)),1)
        self.assertTrue(zero_through(primitive(z).derivative(0)-z,4))
        with self.assertRaises(ValueError):exponential(jet(1))

    def test_reference_fields_frozen_in_F_and_U_not_in_E(self):
        r=reference();self.assertTrue(zero_through(2*r['E'].derivative(0)-r['E'],3))
        self.assertTrue(zero_through(r['U'].derivative(0),3))
        self.assertNotEqual(r['E'].derivative(0).value(),0)

    def test_reference_physical_shear_zero_but_integrated_p1_nonzero(self):
        for center in (Q(-1),Q(0),Q(1,3),Q(1)):
            r=reference(center)
            self.assertTrue(zero_through(1-2*r['E'].derivative(0)/r['E'],3))
            self.assertGreater(r['state']['p1'].value(),0)

    def test_actual_slopes_use_integrated_reference_coordinates(self):
        a=prescribed();r=a['reference'];k=Q(1,10)
        self.assertEqual((1-2*a['E'].derivative(0)/a['E']).value(),k*r['state']['p1'].value())
        self.assertEqual(a['U'].derivative(0).value(),-k*(r['X']*r['state']['ns']).value()/2)

    def test_naive_post_taper_copy_of_first_collar_formula_detected(self):
        self.assertTrue(stats()['shortcut_detected'])

    def test_all_240_exact_prescribed_shear_coefficients(self):
        self.assertEqual(stats()['counts']['shear_coefficients'],240)
        self.assertTrue(stats()['all_flows'])

    def test_all_600_exact_five_moment_flow_coefficients(self):
        self.assertEqual(stats()['counts']['moment_coefficients'],600)
        self.assertTrue(stats()['all_flows'])

    def test_all_120_full_pressure_memory_coefficients(self):
        self.assertEqual(stats()['counts']['pressure_coefficients'],120)
        self.assertTrue(stats()['all_memory'])

    def test_independent_216_raw_second_jet_stress_components(self):
        self.assertEqual(stats()['counts']['old_coordinate_components'],216)
        self.assertLess(stats()['maximum_old_coordinate_error'],mp.mpf('1e-55'))

    def test_axis_datum_has_nonzero_parameter_derivative(self):
        a=prescribed();self.assertNotEqual(a['axis'].partial((0,1)),0)
        self.assertEqual((a['state']['Pi']-a['mom'][4]).partial((0,1)),a['axis'].partial((0,1)))

    def test_pressure_axis_datum_changes_reference_Ns_by_exact_full_term(self):
        self.assertTrue(stats()['axis_term'])

    def test_omitting_pressure_datum_changes_actual_axial_slope(self):
        a=prescribed();b=prescribed(drop_axis=True)
        self.assertNotEqual(a['U'].partial((1,0)),b['U'].partial((1,0)))

    def test_actual_pressure_second_radial_jet_not_reference_replacement(self):
        a=prescribed()
        self.assertNotEqual(a['state']['Pi'].partial((2,0)),a['reference']['state']['Pi'].partial((2,0)))

    def test_all_five_nonzero_incoming_moments_preserved(self):
        a=prescribed()
        for m,ref in zip(a['mom'],a['reference']['mom']):
            self.assertNotEqual(m.value(),0)
            for k in range(4):self.assertEqual(m.partial((0,k)),ref.partial((0,k)))

    def test_zero_kappa_allowed_without_inverse_and_keeps_reference_fields(self):
        a=prescribed(kappa=Q(0));r=a['reference']
        self.assertTrue(zero_through(a['E']-r['E'],3));self.assertTrue(zero_through(a['U']-r['U'],3))

    def test_final_shear_shape_log_slope_and_ell_are_distinct(self):
        self.assertTrue(all(join_review()['checks'].values()))
        self.assertEqual(Q(1,2)*(1-Q(4,5)),Q(1,10))
        self.assertNotEqual(Q(1,10),Q(3,5))

    def test_flat_shape_endpoint_jets_match_both_adjacent_power_slopes(self):
        # All positive jets of sigma vanish at either flat endpoint: insert
        # its exact endpoint jets, without asserting a radial analytic series.
        y=modules()[0].Jet.variable(0,0,n=2,order=5);ell=jet(Q(7,3));lf=jet(Q(-2,5))
        for value in (0,1):
            shape=y/10+(1-value)*ell+value*lf
            self.assertEqual(shape.partial((1,0)),Q(1,10))
            for r in range(2,6):self.assertEqual(shape.partial((r,0)),0)

    def test_shape_amplitude_radius_identity_has_no_moved_endpoint(self):
        # log [P*(X/XR)^(1/10)] at X=110*exp(z),
        # log XR=log110+10(log C+log P), equals -log C+z/10.
        for logC,logP,z in [(Q(10),Q(3),Q(7)),(Q(10000),Q(9),Q(2000))]:
            self.assertEqual(logP+(z-10*(logC+logP))/10,-logC+z/10)

    def test_frozen_record_replays_and_all_historical_bindings_pass(self):
        validate(self.record);f=json.loads((HERE/'evidence.json').read_text());validate(f);compare_records(self.record,f)
        self.assertEqual(len(self.record['inherited_provenance']),785)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_independent_comparison_sources_bound_directly(self):
        self.assertEqual(len(self.record['source_hashes']),16)
        self.assertIn('results/compact_jet_envelope/jets.py',self.record['source_hashes'])

    def test_false_unconditional_scope_flags_cannot_be_promoted(self):
        for name,value in FLAGS.items():
            if not value:
                r=deepcopy(self.record);r['flags'][name]=True
                with self.assertRaises(ValueError):validate(r)

    def test_formula_tampering_rejected_even_behind_passed_gates(self):
        r=deepcopy(self.record);r['reviews']['defining_fields_full_memory_and_smooth_shape_join']['full_pressure']='Pi=0'
        with self.assertRaises(ValueError):validate(r)

    def test_same_count_provenance_key_substitution_rejected(self):
        r=generate();r['inherited_provenance'].pop(next(iter(r['inherited_provenance'])));r['inherited_provenance']['invented:key']=True
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('invented:key',parent_review()[1])

    def test_global_obligations_and_gates_cannot_be_erased(self):
        for key in ('remaining_obligations','gates','source_hashes'):
            r=deepcopy(self.record)
            if key=='remaining_obligations':r[key].pop()
            else:r[key].pop(next(iter(r[key])))
            with self.assertRaises(ValueError):validate(r)

    def test_cli_refuses_existing_output(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'record.json';p.write_text('preserve')
            with self.assertRaises(SystemExit):main(['--output',str(p)])
            self.assertEqual(p.read_text(),'preserve')


if __name__=='__main__':unittest.main()
