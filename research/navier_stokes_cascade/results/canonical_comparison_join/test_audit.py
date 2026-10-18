"""Exact moment and derivative controls, omission checks and frozen scope."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
import audit
from join_bounds import (unit,kernel_unit,kernel_at,directional_kernel,contract,trace_free,
    moment_remainder,moment_review,canonical_interface,universal_bounds,universal_checks)

class CanonicalComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate()

    def test_canonical_force_parent_not_historical_comparison_field(self):
        self.assertEqual(audit.protocol()['parent_commit'],'634356dff6d7b6c0b350c714cd11403d52c216ae')
        self.assertTrue(all(canonical_interface(audit.parent_record())['checks'].values()))
    def test_same_force_pressure_cancels_once(self):
        self.assertIn('only AFTER same force cancels',audit.protocol()['pressure'])
        self.assertTrue(universal_checks().local_review()['checks']['full_forced_difference_pressure_keeps_minus_div_H'])
    def test_full_difference_energy_identity_detects_all_omissions(self):
        r=universal_checks().local_review();self.assertEqual(r['exact_stress_components'],9)
        self.assertEqual(r['exact_weighted_local_energy_identities'],3)
        self.assertEqual(r['nonzero_terms'],['dissipation','force','laplacian','pressure','strain','time','transport'])
    def test_radial_kernel_sign_and_trace(self):
        self.assertEqual(kernel_unit((0,0,1)),[[-1,0,0],[0,-1,0],[0,0,2]])
    def test_kernel_scaling_exact(self):
        for R in (1,2,7):self.assertEqual(kernel_at((0,0,R)),[[F(q,R**3) for q in row] for row in kernel_unit((0,0,1))])
    def test_radial_directional_derivative_is_minus_three_kernel(self):
        for n in ((1,0,0),(0,0,1),(F(3,5),F(4,5),0)):
            self.assertEqual(directional_kernel(n,n),[[-3*q for q in row] for row in kernel_unit(n)])
    def test_tangential_directional_derivative_has_offdiagonal_three(self):
        self.assertEqual(directional_kernel((0,0,1),(1,0,0)),[[0,0,3],[0,0,0],[3,0,0]])
    def test_derivative_trace_and_bound_controls(self):
        r=moment_review();self.assertEqual(r['derivative_controls'],16)
        self.assertTrue(r['checks']['directional_derivative_bound_24']);self.assertTrue(r['checks']['differentiated_kernel_trace_zero'])
    def test_invalid_directions_are_rejected(self):
        for n in ((0,0,0),(1,1,0),(1,0)):
            with self.assertRaises(ValueError):unit(n)
    def test_origin_and_irrational_fixture_radius_rejected(self):
        for z in ((0,0,0),(1,1,0)):
            with self.assertRaises(ValueError):kernel_at(z)
    def test_remainder_has_exact_384_constant(self):self.assertEqual(moment_remainder(2,1,1),24)
    def test_remainder_scales_linearly_in_mass(self):self.assertEqual(moment_remainder(7,1,3),3*moment_remainder(7,1,1))
    def test_remainder_zero_mass_or_support(self):self.assertEqual(moment_remainder(2,0,7),0);self.assertEqual(moment_remainder(2,1,0),0)
    def test_insufficient_radius_or_negative_mass_rejected(self):
        for args in ((1,1,1),(F(1,2),0,1),(2,-1,1),(2,1,-1)):
            with self.assertRaises(ValueError):moment_remainder(*args)
    def test_all_162_exact_point_remainders(self):
        r=moment_review();self.assertEqual(r['point_remainder_controls'],162)
        self.assertTrue(all(abs(F(a['error']))<=F(a['bound']) for a in r['rows']))
    def test_point_remainders_include_genuinely_nonzero_errors(self):
        self.assertTrue(any(F(r['error']) for r in moment_review()['rows']))
    def test_isotropic_moment_kernel_contraction_is_zero(self):
        M=[[F(5*int(i==j)) for j in range(3)] for i in range(3)]
        for n in ((1,0,0),(0,0,1),(F(2,3),F(1,3),F(2,3))):self.assertEqual(contract(kernel_unit(n),M),0)
    def test_compact_divergence_free_tensor_witness_is_anisotropic(self):
        r=moment_review();self.assertEqual(r['axis_value_without_inverse_four_pi'],'-4')
        self.assertTrue(r['checks']['nonzero_trace_free_moment']);self.assertFalse(r['actual_same_force_NS_trajectories'])
    def test_trace_free_projection_preserves_kernel_contraction(self):
        M=[[F(2*int(i==j and i<2)) for j in range(3)] for i in range(3)]
        self.assertEqual(trace_free(M),[[F(2,3),0,0],[0,F(2,3),0],[0,0,F(-4,3)]])
        self.assertEqual(contract(kernel_unit((0,0,1)),trace_free(M)),-4)
    def test_drop_mixed_moment_is_not_licensed(self):self.assertFalse(audit.FLAGS['compactness_and_divergence_freedom_assumed_to_force_moment_isotropy'])
    def test_generic_and_isotropic_flux_powers_differ(self):
        r=moment_review()['checks'];self.assertTrue(r['mixed_flux_generic_decay']);self.assertTrue(r['mixed_flux_isotropic_decay_only_if_moment_isotropic'])
    def test_retained_quadratic_Young_powers_and_full_radius_two(self):
        r=universal_bounds().flux_rows();self.assertEqual([a['young'] for a in r['rows']],['4','14/5']);self.assertEqual(r['worst'],'2')
    def test_eighth_cutoff_power_required_by_pairing(self):
        self.assertTrue(universal_bounds().old.cutoff_weights(8)['valid']);self.assertFalse(universal_bounds().old.cutoff_weights(7)['valid'])
    def test_pressure_kernel_both_endpoints_rejected(self):
        for p in (1,F(3,2)):
            with self.assertRaises(ValueError):universal_bounds().kernel(p)
    def test_pressure_recovery_strict_negative_sobolev_order(self):
        self.assertIn('H^-2',audit.protocol()['recovery']);self.assertIn('H^-3',audit.protocol()['recovery'])
        self.assertTrue(universal_bounds().pressure_review()['checks']['L1_Fourier_bound_embeds_in_H_minus_two_not_critical_H_minus_three_halves'])
    def test_dissipation_integrated_budget_algebra(self):
        # X stands for exp(2*L*T); use an independent rational placeholder.
        for C,L,T,t,X in ((F(3),F(0),F(1,2),F(1,4),F(1)),(F(2),F(3),F(2,3),F(1,3),F(17))):
            self.assertEqual(2*L*t*(2*C*T*X)+2*C*t,2*C*t*(1+2*L*T*X))
    def test_same_localization_radius_in_terminal_threshold(self):self.assertIn('rL^2/(8*Xin)',audit.protocol()['core'])
    def test_no_terminal_uniformity_or_global_competitor_derivative_assumption(self):
        self.assertFalse(audit.FLAGS['comparison_constants_claimed_uniform_through_terminal_time'])
        self.assertTrue(audit.FLAGS['no_global_competitor_gradient_dissipation_or_energy_inequality_hypothesis_added'])
    def test_parent_scope_is_preserved_after_join(self):
        self.assertFalse(audit.parent_record()['flags']['canonical_global_competitor_pressure_recovery_and_comparison_joined'])
        self.assertTrue(audit.FLAGS['canonical_global_competitor_comparison_and_exclusion_relative_to_retained_premises'])
    def test_frozen_record_matches_fresh_and_json_roundtrip(self):
        r=json.loads((audit.HERE/'evidence.json').read_text());audit.validate(r);audit.compare_records(self.record,r)
        audit.compare_records(self.record,json.loads(json.dumps(self.record)))
    def test_unconditional_scope_promotion_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_new_moment_error_bound_tamper_rejected(self):
        r=deepcopy(self.record);r['reviews']['new_compact_mixed_moment_controls']['rows'][0]['bound']='0'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_false_moment_cancellation_tamper_rejected(self):
        r=deepcopy(self.record);r['reviews']['new_compact_mixed_moment_controls']['axis_value_without_inverse_four_pi']='0'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_same_count_historical_binding_substitution_rejected(self):
        r=deepcopy(self.record);b=r['inherited_provenance']['bindings'];b['fake:key']=b.pop(next(iter(b)))
        with self.assertRaises(ValueError):audit.validate(r)
    def test_historical_bindings_include_canonical_force_and_comparison_sources(self):
        r=self.record['inherited_provenance'];self.assertGreater(len(r['bindings']),1202)
        self.assertTrue(all(r['bindings'].values()))
    def test_source_tamper_rejected(self):
        r=deepcopy(self.record);r['source_hashes'][next(iter(r['source_hashes']))]='0'*64
        with self.assertRaises(ValueError):audit.validate(r)
    def test_remaining_foundation_obligation_cannot_be_erased(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop()
        with self.assertRaises(ValueError):audit.validate(r)
    def test_cached_reviews_not_mutable_aliases(self):
        r=audit.reviews();r['new_compact_mixed_moment_controls']['rows'][0]['bound']='bad'
        self.assertNotEqual(audit.reviews()['new_compact_mixed_moment_controls']['rows'][0]['bound'],'bad')
    def test_cli_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('keep')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'keep')
if __name__=='__main__':unittest.main()
