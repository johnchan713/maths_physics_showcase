"""Exact shape/scale algebra and correction guards; not actual profile verification."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import json,tempfile,unittest
from pathlib import Path
import audit
import gate_bounds as b

class GateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record=audit.generate();cls.r=cls.record['reviews'];cls.p=b.fixture();cls.out=b.affine_gate(**cls.p);cls.rows=cls.out['rows']
    def check(self,k):self.assertTrue(self.r['checks'][k])
    def symbolic(self,k):self.assertTrue(self.r['symbolic']['checks'][k])
    def test_retained_tc_not_native_t1(self):self.symbolic('retained_cutoff_log_width_is_tc_over_sixteen')
    def test_first_width_cost_ten_sqrt_log(self):self.symbolic('corrected_first_derivative_cost_has_ten_sqrt_logC')
    def test_second_width_cost_has_g_squared_and_g(self):self.symbolic('corrected_second_cost_has_distinct_g_squared_and_g_terms')
    def test_all_four_exact_width_corrections(self):
        self.assertEqual(len(self.r['width_correction_controls']),4);self.check('correct_retained_width_chain_and_viscous_shape_coefficient')
    def test_retained_additive_bump_jets(self):self.assertEqual(b.cutoff().bump_bounds(),(320,121600,220000000))
    def test_value_jet_has_no_shape_derivative(self):self.assertEqual(self.out['jets']['value'][1:],(0,0))
    def test_first_jet_has_no_second_shape_cost(self):self.assertEqual(self.out['jets']['radial'][2],0)
    def test_first_jet_inverse_tc(self):
        p=self.p;self.assertEqual(self.out['jets']['radial'][1],16*p['core']['value']/(p['tc']*p['xkeep']))
    def test_second_jet_keeps_log_jacobian(self):
        p=self.p;d=p['tc']/16;x=p['xkeep'];self.assertEqual(self.out['jets']['radial2'][1],2*p['core']['radial']/(d*x)+p['core']['value']/(d*x*x))
    def test_second_jet_inverse_tc_squared(self):
        p=self.p;self.assertEqual(self.out['jets']['radial2'][2],p['core']['value']/((p['tc']/16)**2*p['xkeep']**2))
    def test_all_twelve_complete_interface_controls(self):self.assertEqual(len(self.r['exact_composition_controls']),12);self.check('all_twelve_exact_complete_operator_compositions')
    def test_four_affine_rows_not_unspecified_Cstar(self):self.assertEqual(set(self.rows),{'unit','value','derivative','stress'})
    def test_all_affine_coefficients_nonnegative(self):self.check('all_affine_coefficients_nonnegative')
    def test_unit_row_exactly_one(self):self.assertEqual(self.rows['unit'],(1,0,0))
    def test_value_row_has_borel_half(self):
        p=self.p;S=self.out['jets']['value'][0];self.assertEqual(self.rows['value'],(5*(2*p['c']*S+F(1,2))/p['a'],0,0))
    def test_derivative_row_slope_keeps_two_radial_measures(self):
        c=self.p['c'];self.assertEqual(self.rows['derivative'][1],10*c*(2*c+1)*self.out['jets']['radial'][1])
    def test_derivative_row_has_no_L2(self):self.assertEqual(self.rows['derivative'][2],0)
    def test_stress_row_keeps_borel_half(self):
        p=self.p;self.assertEqual(self.rows['stress'][0],F(1,2)+(p['c']**2*self.out['theta'][0]/2+p['c']*self.out['axial'][0])/p['tau'])
    def test_positive_viscous_L2_coefficient(self):self.check('strictly_positive_L2_stress_cost');self.symbolic('exact_positive_L2_stress_coefficient')
    def test_numeric_exact_L2_formula(self):
        p=self.p;self.assertEqual(self.rows['stress'][2],p['nu']*p['c']**3*(p['c']+2)*p['core']['value']/(p['tau']*(p['tc']/16)**2*p['xkeep']**2))
    def test_affine_not_quadratic_shape_dependence(self):self.symbolic('no_L1_L2_product_or_quadratic_shape_terms')
    def test_omitted_viscous_shape_cost_detected(self):self.symbolic('missing_viscous_shape_cost_detected')
    def test_complete_inhomogeneous_terms_not_counted_in_slopes(self):self.check('complete_fixed_pressure_swirl_and_axial_costs_retained')
    def test_future_signed_pressure_memory_retains_both_repair_rows(self):self.check('future_signed_pressure_memory_enters_both_normalized_repair_rows')
    def test_explicit_first_shape_source_slopes(self):
        p=self.p;X=p['c']**2;L=1-2*p['h'];u=p['leading'][1];v=p['leading'][4]
        x1=self.out['jets']['radial'][1];z1=self.out['jets']['radial2'][1]
        self.assertEqual(self.out['theta'][1],(X*(1+2*u)/L+4*p['nu']+v)*x1+2*p['nu']*X*z1)
        self.assertEqual(self.out['axial'][1],(X*(1+2*u)/L+2*p['nu']+v)*x1+2*p['nu']*X*z1)
    def test_each_inhomogeneous_source_omission_individually_detected(self):
        for i in (5,6,7):
            lead=list(self.p['leading']);lead[i]=0;q=b.affine_gate(**dict(self.p,leading=lead))
            self.assertLess(q['rows']['stress'][0],self.rows['stress'][0]);self.assertEqual(q['rows']['stress'][1:],self.rows['stress'][1:])
    def test_wider_wrong_width_underestimates_shape_costs(self):self.check('wider_mislabelled_tc_underestimates_both_shape_slopes')
    def test_exact_boundary_and_strict_outside(self):self.check('exact_closed_shape_region_and_strict_outside')
    def test_all_sixteen_independent_geometry_viscosity_compositions(self):
        for a,c,nu,r in product((F(1,2),F(1)),(F(2),F(4)),(F(1,3),F(2)),(F(1,2),F(2))):
            p=dict(self.p,a=a,c=c,nu=nu,rho_root=r);o=b.affine_gate(**p);L1,L2=F(13,2),F(123,4)
            jets=b.cutoff().propagate(p['core'],p['coefficient'],a,p['xkeep'],p['tc'],L1,L2,r)['extended']
            s=b.first().source_bounds(p['h'],nu,a,c,jets,*p['leading']);g=b.first().compact_gate(a,c,jets,s,p['tau'])
            self.assertEqual(b.cstar(o['rows'],L1,L2),g['Cval'])
    def test_shape_caps_are_the_exact_row_rearrangement(self):
        C=b.cstar(self.rows,32,1216);cap=b.shape_caps(self.rows,2*C,32)
        self.assertEqual(b.evaluate(self.rows['derivative'],cap['L1_upper'],0),2*C)
        self.assertEqual(b.evaluate(self.rows['stress'],32,cap['L2_upper_at_L1']),2*C)
    def test_lower_corner_failure_is_only_certificate_failure(self):self.check('lower_shape_corner_is_only_necessary')
    def test_sub_lower_shape_norms_never_accepted(self):
        self.assertFalse(b.shape_region(self.rows,10**100,0,4)['accepted']);self.assertFalse(b.shape_region(self.rows,10**100,1,3)['accepted'])
    def test_compressed_shape_fixed_certificate_obstruction(self):self.check('compressed_admissible_shapes_defeat_every_finite_fixed_certificate')
    def test_root_free_compression_integer_witness_for_tiny_and_large_budgets(self):
        for B in (F(1,10**50),F(1),F(10**100)):
            r=b.compression_barrier(self.rows,B);self.assertIs(type(r['compression']),int);self.assertTrue(r['exceeds_fixed_budget']);self.assertGreater(r['stress_row_lower'],B)
    def test_same_existing_tilt_not_reselected(self):self.check('same_tilt_positive_radius_and_two_band_dyadic_witness');self.assertFalse(audit.FLAGS['actual_cutoff_shape_existing_tilt_or_retained_scale_acceptance_verified'])
    def test_exact_dyadic_log_boundaries(self):
        for k in (1,2,11,1000):
            self.assertEqual(b.ceil_log2_positive(F(2)**k),k);self.assertEqual(b.ceil_log2_positive(F(2)**k+F(1,3)),k+1)
    def test_subunit_ratio_uses_zero_safe_exponent(self):
        for q in (F(1,10**100),F(1,2),F(1)):self.assertEqual(b.ceil_log2_positive(q),0)
    def test_rational_log_boundary_with_denominator(self):self.assertEqual(b.ceil_log2_positive(F(17,3)),3);self.assertEqual(b.ceil_log2_positive(F(8,3)),2)
    def test_dyadic_witness_whole_chart_plus_two(self):
        w=self.r['dyadic_witness'];self.assertEqual(w['sufficient_ell'],2+100*w['dyadic_exponent']);self.assertTrue(all(w['checks'].values()))
    def test_dyadic_witness_fractional_scale_exponent_ceiling(self):
        w=b.dyadic_witness(self.rows,32,1216,F(1,68992),F(3,1000));self.assertGreaterEqual(2*F(3,1000)*(w['sufficient_ell']-2),w['dyadic_exponent'])
    def test_dyadic_witness_already_small_ratio(self):
        C=b.cstar(self.rows,32,1216);w=b.dyadic_witness(self.rows,32,1216,2*C,F(1,200));self.assertEqual(w['sufficient_ell'],2);self.assertEqual(w['dyadic_exponent'],0)
    def test_invalid_geometry_and_negative_norms(self):
        for changes in (dict(tc=0),dict(a=0),dict(c=F(1,2)),dict(rho_root=0),dict(coefficient=-1),dict(xkeep=F(1,2)),dict(tau=0),dict(nu=0)):
            with self.assertRaises(ValueError):b.affine_gate(**dict(self.p,**changes))
    def test_invalid_h_and_leading_length(self):
        for changes in (dict(h=0),dict(h=F(1,100)),dict(leading=[]),dict(leading=[-1]*8)):
            with self.assertRaises(ValueError):b.affine_gate(**dict(self.p,**changes))
    def test_non_exact_inputs_rejected(self):
        for value in (True,.1):
            with self.assertRaises(TypeError):b.affine_gate(**dict(self.p,tc=value))
            with self.assertRaises(TypeError):b.ceil_log2_positive(value)
    def test_invalid_shape_and_log_domain(self):
        for L1,L2 in ((0,4),(1,0)):
            with self.assertRaises(ValueError):b.cstar(self.rows,L1,L2)
        for q in (0,-1):
            with self.assertRaises(ValueError):b.ceil_log2_positive(q)
    def test_invalid_budget_caps_and_witness(self):
        with self.assertRaises(ValueError):b.shape_region(self.rows,0,32,1216)
        with self.assertRaises(ValueError):b.shape_caps(self.rows,1,0)
        with self.assertRaises(ValueError):b.dyadic_witness(self.rows,32,1216,0,F(1,200))
        with self.assertRaises(ValueError):b.compression_barrier(self.rows,-1)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'4de29921a9af3fc05e5ade7badb9a0c2361919ba');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()))
    def test_frozen_record(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_bytes_preserved(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_all_unsupported_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_scope_type_mutation_rejected(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_width_correction_cannot_be_erased(self):
        r=deepcopy(self.record);r['protocol']['correction']='tc=C^-20'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_affine_slope_mutation_rejected(self):
        r=deepcopy(self.record);r['reviews']['affine_output']['rows']['stress'][2]='0'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_remaining_actual_shape_obligation_not_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_cli_fresh_replay_and_refusal_to_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
