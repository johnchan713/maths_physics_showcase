#!/usr/bin/env python3
"""Signed roots, moving phases, same-root jets, honest powers and full memory."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from math import comb,factorial
from pathlib import Path
import tempfile
import unittest
from all_order_maps import frequency_powers,coefficient_bounds,map_review,transfer_review
from map_checks import (jet,variable,zero,sine,signed_F,variance_series,root_jets,phase_jets,
    repair_jets,matvec,quadratic,B,invert_matrix,frequency_jets,diagnostic_review)
from audit import HERE,FLAGS,generate,validate,compare_records,parent_review,main


class AllOrderMapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=generate()

    def test_current_parent_and_same_accepted_frequency(self):
        self.assertEqual(self.record['protocol']['parent_commit'],'9f4ae8973fdd1384530085505010597033e206cf')
        self.assertEqual(self.record['protocol']['unchanged_frequency'],'N=1+floor(H^32)')

    def test_whole_compact_smooth_input_not_silently_supplied_by_inner_parent(self):
        self.assertIn('WHOLE',self.record['protocol']['input'])
        self.assertFalse(FLAGS['whole_annular_exterior_and_global_profile_foundation_independently_verified'])

    def test_no_global_parameter_analyticity_or_radial_smallness_promotion(self):
        self.assertFalse(FLAGS['modulated_compact_profile_globally_complex_eta_analytic'])
        self.assertFalse(FLAGS['all_radial_derivative_errors_small_or_N_uniform'])

    def test_all_order_map_and_transfer_ledgers_pass(self):
        self.assertTrue(all(map_review()['checks'].values()));self.assertTrue(all(transfer_review()['checks'].values()))

    def test_frequency_powers_distinguish_fields_from_moment_primitives(self):
        for r in (0,1,2,8,1000):
            p=frequency_powers(r);self.assertEqual(p['field'],r-1);self.assertEqual(p['moment'],max(r-2,-1))
        self.assertEqual(frequency_powers(1)['moment'],-1)

    def test_invalid_radial_orders_and_inverse_bounds_rejected(self):
        for r in (-1,True,Q(1,2)):
            with self.assertRaises(ValueError):frequency_powers(r)
        with self.assertRaises(ValueError):coefficient_bounds(0,1,[1])
        with self.assertRaises(ValueError):coefficient_bounds(1,1,[-1])

    def test_removable_variance_series_has_exact_positive_zero_limit(self):
        self.assertEqual(variance_series()[0],Q(1,2))
        self.assertEqual(variance_series()[1],0)
        self.assertTrue(all(variance_series()[i]==0 for i in (1,3,5)))

    def test_signed_variance_root_84_coefficients_match_Lagrange_inversion(self):
        for p in (Q(0),Q(1,3),Q(-1,3)):
            a,b,res,_,_=root_jets(p)
            self.assertTrue(zero(a-b));self.assertTrue(zero(res))
        self.assertEqual(diagnostic_review()['root_coefficients'],84)

    def test_signed_zero_pressure_root_has_no_cusp_or_zero_denominator(self):
        mu=variable(0);self.assertTrue(zero(signed_F(mu,jet(0))-mu))
        self.assertTrue(zero(signed_F(-mu,jet(Q(1,3)))+signed_F(mu,jet(Q(1,3)))))
        self.assertEqual(signed_F(mu,jet(Q(1,3))).partial((1,0)),1)

    def test_pressure_variation_is_retained_in_high_root_jets(self):
        a,_,_,p,s=root_jets(Q(1,3))
        self.assertFalse(zero(a-s))
        self.assertNotEqual(a.partial((2,1)),0)

    def test_moving_phase_84_coefficients_match_independent_Newton_inverse(self):
        theta,newton,res,frozen=phase_jets();self.assertTrue(zero(theta-newton));self.assertTrue(zero(res))
        self.assertEqual(diagnostic_review()['phase_coefficients'],84)

    def test_frozen_inverse_phase_misses_nonzero_mixed_derivative(self):
        theta,_,_,frozen=phase_jets()
        self.assertNotEqual(theta.partial((1,0,1)),0)
        self.assertEqual(frozen.partial((1,0,1)),0)

    def test_exact_five_by_five_inverse_has_no_sampled_determinant(self):
        inv=invert_matrix(B)
        for i in range(5):
            for j in range(5):self.assertEqual(sum(inv[i][k]*B[k][j] for k in range(5)),Q(i==j))

    def test_same_quadratic_repair_180_coefficients_match_fixed_point_series(self):
        for N in (16,32,64,128):
            c,other,target=repair_jets(N)
            self.assertTrue(all(zero(a-b) for a,b in zip(c,other)))
        self.assertEqual(diagnostic_review()['repair_coefficients'],180)

    def test_all_five_nonlinear_repair_rows_vanish_through_order_eight(self):
        c,_,target=repair_jets()
        self.assertTrue(all(zero(a+b-z) for a,b,z in zip(matvec(B,c),quadratic(c),target)))

    def test_omitted_U_square_energy_term_detected(self):
        c,_,target=repair_jets();wrong=quadratic(c,drop_u_square=True)
        self.assertFalse(zero((matvec(B,c)[3]+wrong[3]-target[3])))

    def test_omitted_pressure_square_term_detected(self):
        c,_,target=repair_jets();wrong=quadratic(c,drop_pressure_square=True)
        self.assertFalse(zero((matvec(B,c)[4]+wrong[4]-target[4])))

    def test_high_repair_derivative_bound_uses_only_lower_pairs_after_Jacobian(self):
        d=coefficient_bounds(2,3,[1]*9)
        self.assertEqual(d[1],4)
        self.assertEqual(d[2],4*(1+3*d[1]**2))
        self.assertEqual(len(d),9)

    def test_fixed_N_does_not_require_all_derivative_norms_small(self):
        N=32;k=100;catalan=Q(comb(2*(k-1),k-1),k)
        raw=factorial(k)*catalan/Q(N**k)
        self.assertLess(Q(1,N),Q(1,10));self.assertGreater(raw,10**50)
        # c+c^2=eta/N has Jacobian 1 at eta=0 at every order.

    def test_exact_chain_rule_all_135_mixed_coefficients_includes_lower_slow_terms(self):
        for N in (16,32,64):self.assertTrue(frequency_jets(N)['chain'])
        self.assertEqual(diagnostic_review()['chain_coefficients'],135)

    def test_third_field_derivative_can_grow_quadratically_with_frequency(self):
        values=[abs(frequency_jets(N)['field'].partial((3,0))) for N in (16,64)]
        self.assertGreater(values[1],10*values[0])

    def test_moment_radial_gain_keeps_second_derivative_order_one(self):
        for N in (16,32,64):
            m=frequency_jets(N)['moment'];self.assertEqual(m.partial((2,0)),1)
            self.assertEqual(m.partial((4,0)),-N**2)
            self.assertEqual(m.partial((8,0)),-N**6)

    def test_full_incoming_moment_memory_survives_every_local_flow(self):
        for N in (16,32,64):
            a=frequency_jets(N);self.assertTrue(zero(a['moment'].derivative(0)-a['density']))
            self.assertEqual(a['moment'].partial((0,1)),Q(1,N))
            self.assertEqual(a['moment'].partial((0,2)),Q(2,7*N))

    def test_fifth_row_defect_can_have_zero_value_but_nonzero_pressure_gradient(self):
        eta=variable(0,n=1,order=8);c=eta+eta*eta/7;A=Q(101,200)
        defect=4*A*eta*c-(1-eta*eta)*c.derivative(0)
        self.assertEqual(c.value(),0);self.assertEqual(defect.value(),-1)
        self.assertNotEqual(c.partial((1,)),0)

    def test_frozen_record_JSON_roundtrip_and_all_inherited_bindings(self):
        f=json.loads((HERE/'evidence.json').read_text());validate(f);compare_records(self.record,f)
        validate(json.loads(json.dumps(self.record)))
        self.assertEqual(len(self.record['inherited_provenance']),801)
        self.assertTrue(all(self.record['inherited_provenance'].values()))

    def test_direct_source_and_dependency_record_bindings(self):
        self.assertEqual(len(self.record['source_hashes']),17)
        self.assertIn('results/stress_realization_audit/loop.py',self.record['source_hashes'])

    def test_unconditional_scope_promotion_rejected(self):
        for key,value in FLAGS.items():
            if not value:
                r=deepcopy(self.record);r['flags'][key]=True
                with self.assertRaises(ValueError):validate(r)

    def test_tampered_pressure_memory_behind_passed_gates_rejected(self):
        r=deepcopy(self.record);r['reviews']['frequency_powers_full_memory_and_input_limits']['full_pressure']='Pi=0'
        with self.assertRaises(ValueError):validate(r)

    def test_generated_record_does_not_alias_cached_provenance(self):
        r=generate();r['inherited_provenance'].pop(next(iter(r['inherited_provenance'])));r['inherited_provenance']['invented:key']=True
        with self.assertRaises(ValueError):validate(r)
        self.assertNotIn('invented:key',parent_review()[1])

    def test_generated_diagnostic_review_does_not_alias_cached_canonical_checks(self):
        r=generate();key='fixed_phase_chain_includes_every_N_power'
        r['reviews']['independent_exact_higher_jet_controls']['checks'][key]=False
        with self.assertRaises(ValueError):validate(r)
        self.assertTrue(diagnostic_review()['checks'][key])

    def test_remaining_input_obligation_gate_or_binding_cannot_be_erased(self):
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
