#!/usr/bin/env python3
"""Regression and deliberate-failure checks; not a theorem verifier."""
import json
from pathlib import Path
import unittest
import mpmath as mp
from loop import (positive_series_iv,variance,variance_iv,me_minus_one,lower,upper,
                  root_bracket,shear_ratio,cone_gaps,zero_pressure_primitives)
from repair import Repair,divided_weight,precondition_target
from bounds import (repair_certificate,uniform_mu_log_bound,active_gap_bounds,
                    cone_tolerance,frequency_requirement)
from audit import FLAGS,STATUS,compare,inherited_provenance,modulation_controls,validate_protocol


class Checks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps=80;mp.iv.dps=110
        cls.repair=Repair('.01',32,8)

    def small(self,value,tolerance='1e-60'):
        self.assertLess(abs(value),mp.mpf(tolerance))

    def test_me_matches_independent_bessel(self):
        for z in ('0','1e-30','.3','2','10','64'):
            x=mp.mpf(z);reference=mp.besseli(0,x)
            self.small((1+me_minus_one(x))/reference-1)
            enclosure=positive_series_iv(x)
            self.assertLessEqual(lower(enclosure),reference+mp.mpf('1e-75')*reference)
            self.assertGreaterEqual(upper(enclosure),reference-mp.mpf('1e-75')*reference)

    def test_variance_zero_pressure_exact_limit(self):
        self.assertEqual(variance(3,0,2),18)
        enclosure=variance_iv(3,0,2)
        self.assertLessEqual(lower(enclosure),18)
        self.assertGreaterEqual(upper(enclosure),18)

    def test_variance_tiny_pressure_not_subtracted_away(self):
        for p in ('1e-20','1e-60','1e-100'):
            self.small(variance(3,p)-mp.mpf('4.5'),'1e-37')
        with mp.workdps(30):
            p=mp.mpf('1e-30')
            self.assertEqual((mp.besseli(0,6*p)/mp.besseli(0,3*p)**2-1)/p**2,0)
            self.assertGreater(variance(3,p),4)

    def test_variance_interval_matches_independent_formula(self):
        for mu,p in (('2','1'),('3','.25'),('5','-1')):
            u,z=mp.mpf(mu),mp.mpf(p)
            ref=(mp.besseli(0,2*u*z)/mp.besseli(0,u*z)**2-1)/z**2
            enclosure=variance_iv(u,z)
            self.assertLessEqual(lower(enclosure),ref+mp.mpf('1e-75'))
            self.assertGreaterEqual(upper(enclosure),ref-mp.mpf('1e-75'))

    def test_variance_even_in_pressure(self):
        self.assertEqual(variance(3,mp.mpf('.25')),variance(3,-mp.mpf('.25')))

    def test_variance_monotonic_diagnostic(self):
        for p in (0,mp.mpf('1e-30'),mp.mpf('.25'),1):
            values=[variance(mu,p) for mu in (0,1,2,4,8,16)]
            self.assertTrue(all(a<b for a,b in zip(values,values[1:])))

    def test_exact_root_bracket(self):
        for p in (0,mp.mpf('1e-30'),1,-1):
            a,b=root_bracket('1.46',p)
            self.assertLess(upper(variance_iv(a,p)),mp.mpf('1.46'))
            self.assertGreater(lower(variance_iv(b,p)),mp.mpf('1.46'))

    def test_unreached_cap_is_rejected(self):
        with self.assertRaises(ValueError):root_bracket(100,1,cap=1)

    def test_zero_target_is_zero_loop(self):
        self.assertEqual(root_bracket(0,1),(0,0))
        self.assertEqual(shear_ratio('.3',0,1,'.2'),mp.mpf('.2'))

    def test_removable_pressure_limit(self):
        reference=shear_ratio('.7',3,0,'.2')
        self.small(shear_ratio('.7',3,'1e-40','.2')-reference,'1e-38')

    def test_lift_period_and_positive_density(self):
        for theta in (mp.mpf('.2'),mp.mpf('1.1'),mp.mpf('3.3')):
            f=zero_pressure_primitives(theta)
            shifted=zero_pressure_primitives(theta+2*mp.pi)
            self.small(shifted['phi']-f['phi']-1)
            self.assertGreater(f['phase_derivative'],0)

    def test_primitive_differential_identities(self):
        theta=mp.mpf('.7');f=zero_pressure_primitives(theta,E='1.3')
        da=mp.diff(lambda x:zero_pressure_primitives(x,E='1.3')['A'],theta)/f['phase_derivative']
        db=mp.diff(lambda x:zero_pressure_primitives(x,E='1.3')['B'],theta)/f['phase_derivative']
        self.small(da+(f['a']-mp.mpf('.8'))/2)
        self.small(db-mp.mpf('1.3')*f['b']/2)

    def test_primitives_have_zero_phi_mean(self):
        for key in ('A','B'):
            val=mp.quad(lambda t:zero_pressure_primitives(t)[key]*zero_pressure_primitives(t)['phase_derivative'],
                        [0,mp.pi/2,mp.pi,3*mp.pi/2,2*mp.pi])
            self.small(val)

    def test_weighted_mean_not_unweighted(self):
        a=mp.mpf('.8')
        weighted=mp.quad(lambda t:zero_pressure_primitives(t)['a']*zero_pressure_primitives(t)['phase_derivative'],[0,mp.pi,2*mp.pi])
        raw=mp.quad(lambda t:zero_pressure_primitives(t)['a'],[0,mp.pi,2*mp.pi])/(2*mp.pi)
        self.small(weighted-a)
        self.assertGreater(abs(raw-a),mp.mpf('.1'))

    def test_exponential_denominator_cannot_be_dropped(self):
        c=modulation_controls()
        self.small(c['maximum_identity_error'])
        self.assertTrue(c['omitted_exponential_detected'])

    def test_radial_derivative_not_small(self):
        for N in (10,100,1000):
            derivative=mp.diff(lambda y:mp.sin(2*mp.pi*N*y)/N,0)
            self.small(derivative-2*mp.pi)

    def test_radial_second_derivative_grows(self):
        for N in (10,100,1000):
            derivative=mp.diff(lambda y:mp.sin(2*mp.pi*N*y)/N,mp.mpf(1)/(4*N),2)
            self.small(derivative/N+4*mp.pi**2)

    def test_eta_dependent_phase_breaks_smallness(self):
        for N in (10,100,1000):
            derivative=mp.diff(lambda eta:mp.sin(2*mp.pi*N*eta)/N,0)
            self.small(derivative-2*mp.pi)

    def test_correction_certificate_at_two_precisions(self):
        for digits in (80,110):self.assertTrue(all(repair_certificate(digits)['checks'].values()))

    def test_all_phase_gap_bound(self):
        bounds=active_gap_bounds('.8','.2',10,1,'7.8',1,64,'1e-8')
        self.assertTrue(all(lower(g)>0 for g in bounds['gaps']))
        self.assertTrue(bounds['cap_reachability_is_separate_obligation'])

    def test_excessive_loop_delta_rejected(self):
        with self.assertRaises(ValueError):active_gap_bounds('.8','.2',10,1,'7.8',1,64,'.1')

    def test_uniform_cap_is_explicit_not_a_grid(self):
        b=uniform_mu_log_bound('.8',1,1)
        self.assertGreater(b['guaranteed_variance_lower'],b['required_variance_upper'])
        self.assertGreater(b['log_mu_max'],b['Z'])

    def test_numerical_bump_mass(self):
        for j in range(5):self.small(self.repair.integral(j,lambda y:1)-1,'1e-25')

    def test_normalized_inverse_bound(self):
        B=self.repair.inverse
        self.assertLess(max(mp.fsum(abs(B[i,j]) for j in range(5)) for i in range(5)),1000)

    def test_divided_weight_stable_at_extreme_lambda(self):
        lam=mp.exp(-4*(mp.exp(64)+10))
        self.assertGreater(lam,0)
        self.assertEqual(1-mp.exp(-lam),0)
        self.small(divided_weight(lam,'1.5')-mp.mpf('1.5'))

    def test_stable_matrix_survives_naive_rank_loss(self):
        lam=mp.exp(-4*(mp.exp(64)+10))
        repair=Repair(lam,16,2)
        stable=mp.matrix([[repair.B[i,j] for j in range(2)] for i in range(2)])
        ordinary=mp.matrix([[stable[0,j] for j in range(2)],
                            [stable[0,j]-lam*stable[1,j] for j in range(2)]])
        self.assertGreater(abs(mp.det(stable)),1)
        self.assertEqual(mp.det(ordinary),0)

    def test_finite_limit_does_not_authorize_lambda_zero(self):
        self.assertEqual(divided_weight(0,2),2)
        with self.assertRaises(ValueError):Repair(0)
        with self.assertRaises(ValueError):precondition_target([0]*5,0)

    def test_discrepancy_loss_is_retained(self):
        d=precondition_target(['1e-16',0,0,0,0],'1e-5')
        self.assertEqual(d[1]/d[0],100000)

    def test_rounded_incoming_difference_cannot_be_recovered(self):
        lam=mp.exp(-4*(mp.exp(64)+10))
        lost=precondition_target([1,0,1-lam,0,0],lam)[1]
        self.assertEqual(lost,0)
        # In exact arithmetic (1-(1-lambda))/lambda is one. The target
        # must be carried separately; stabilizing the matrix cannot recover it.
        self.assertNotEqual(lost,1)

    def test_five_moments_numerically_restore(self):
        target=mp.matrix([mp.mpf('1e-16')*j for j in (1,-2,3,-1,2)])
        root=self.repair.solve(target)
        self.small(max(abs(v) for v in self.repair.value(root)-target),'1e-70')

    def test_implicit_eta_derivative_of_small_root(self):
        eps=mp.mpf('1e-16');eta=mp.mpf('.3');step=mp.mpf('1e-14')
        def target(e):return eps*mp.matrix([1+e,e*e,mp.sin(e),1/(1+e*e),e**3])
        root=self.repair.solve(target(eta));jac=self.repair.B.copy()
        for j in range(5):
            unit=mp.matrix(5,1);unit[j]=1
            column=2*self.repair.quadratic(root,unit)
            for i in range(5):jac[i,j]+=column[i]
        rhs=eps*mp.matrix([1,2*eta,mp.cos(eta),-2*eta/(1+eta*eta)**2,3*eta*eta])
        implicit=mp.lu_solve(jac,rhs)
        observed=(self.repair.solve(target(eta+step))-self.repair.solve(target(eta-step)))/(2*step)
        self.small(max(abs(v) for v in implicit-observed),'1e-40')

    def test_independent_delta_integrands(self):
        c=mp.matrix([mp.mpf('1e-6')*j for j in (1,-2,3,-1,2)])
        self.small(max(abs(v) for v in self.repair.direct_changes(c)-self.repair.ordinary_changes(c)))

    def test_disjoint_supports_remove_only_cross_term(self):
        c=mp.matrix([1]*5);q=self.repair.quadratic(c,c)
        self.assertEqual([q[j] for j in range(3)],[0,0,0])
        self.assertNotEqual(q[3],0);self.assertGreater(q[4],0)

    def test_axial_square_is_retained(self):
        q=self.repair.quadratic(mp.matrix([1,0,0,0,0]),mp.matrix([1,0,0,0,0]))
        self.assertLess(q[3],0);self.assertEqual(q[4],0)

    def test_four_moments_do_not_restore_pressure(self):
        root=self.repair.solve([0,0,0,0,mp.mpf('1e-16')])
        value=self.repair.value(root)
        self.small(max(abs(value[j]) for j in range(4)),'1e-70')
        self.assertGreater(value[4],mp.mpf('9e-17'))

    def test_large_target_refused(self):
        with self.assertRaises(ValueError):self.repair.solve([1]*5)

    def test_cone_tolerance_preserves_sampled_state(self):
        base=[mp.mpf(x) for x in ('3','.2','10','.5')]
        gaps=cone_gaps(*base);margin=min(gaps)
        epsilon=cone_tolerance(3,10,margin)
        for signs in ((1,1,1,1),(-1,-1,-1,-1),(1,-1,1,-1)):
            new=cone_gaps(*[x+s*epsilon for x,s in zip(base,signs)])
            self.assertGreaterEqual(min(new),margin/2)

    def test_zero_margin_refused(self):
        with self.assertRaises(ValueError):cone_tolerance(1,10,0)

    def test_conditional_frequency_includes_lambda_loss(self):
        kwargs=dict(epsilon='1e-10',state_constant=1,correction_constant=1,coefficient_tolerance='1e-8')
        a=frequency_requirement(discrepancy_constant=1,**kwargs)
        b=frequency_requirement(discrepancy_constant=100000,**kwargs)
        self.assertGreater(b['sufficient_N_lower'],1000*a['sufficient_N_lower'])
        self.assertFalse(a['actual_profile_inputs_verified'])

    def test_flags_keep_global_gap(self):
        for key in ('actual_compact_input_bounds_instantiated','actual_joined_profile_frequency_selected',
                    'full_admissible_stress_realized','full_PDE_corrections_verified','smooth_force_verified','blowup_verified'):
            self.assertIs(FLAGS[key],False)

    def test_changed_fixture_protocol_is_rejected(self):
        protocol=json.loads((Path(__file__).resolve().parent/'protocol.json').read_text())
        validate_protocol(protocol)
        protocol['fixture_mu_cap']='1'
        with self.assertRaises(ValueError):validate_protocol(protocol)

    def test_inherited_sources_unchanged(self):
        self.assertTrue(all(inherited_provenance().values()))

    def test_status_and_hash_tampering_rejected(self):
        with self.assertRaises(ValueError):compare({'status':'blowup'},{'status':STATUS})
        with self.assertRaises(ValueError):compare({'source_hash':'123'},{'source_hash':'124'})

    def test_false_to_true_flag_tampering_rejected(self):
        with self.assertRaises(ValueError):compare({'flag':True},{'flag':False})


if __name__=='__main__':unittest.main()
