"""Exact interface budgets and diagonal post-terminal force extension scales.

All field conclusions are conditional on inherited uniform physical jet
bounds. No actual terminal jets or full auxiliary atlas are evaluated here.
"""
from fractions import Fraction as F
import math
import mpmath as mp


def order(m):
    if type(m) is not int or m < 0:
        raise ValueError('A nonnegative integer derivative order is required')
    return m


def jet_requirements(m):
    m = order(m)
    return dict(force_jet_order=m, potential_order=m+3,
                direct_swirl_order=m+2, pressure_order=m+1,
                uniform_time_Cauchy_potential_order=m+4,
                uniform_time_Cauchy_swirl_order=m+3,
                uniform_time_Cauchy_pressure_order=m+2)


def endpoint_interface_review():
    return dict(rows=[jet_requirements(m) for m in (0, 1, 3, 12)],
                uniform_endpoint_inequality='sup_[x in E,1-delta<t<1] |D^alpha partial_t^j (A,B,p)| < infinity for all alpha,j',
                existing_input='full-potential physical jet bounds; on q>=c only finitely many stage cutoffs survive',
                Cauchy_estimate='sup_E |D^alpha partial_t^j G(t2)-D^alpha partial_t^j G(t1)| <= |t2-t1| sup |D^alpha partial_t^(j+1) G|',
                no_new_boundary_label_selection=True,
                checks=dict(force_requires_three_potential_derivatives=all(r['potential_order']-r['force_jet_order'] == 3 for r in [jet_requirements(m) for m in (0, 1, 3, 12)]),
                            Cauchy_requires_one_more_time_derivative=True,
                            existing_all_order_physical_bounds_cover_every_fixed_requirement=True,
                            global_atlas_compatibility_remains_inherited_not_independently_audited=True))


def domain_review():
    return dict(upper_q_bound='q<=2*tau+4*|z|^(1/D), 0<h<1/100',
                tau0='qstar/16', z0='(qstar/16)^D',
                support_q_over_qstar_bound=F(2+4, 16),
                remaining_boundary_margin=F(1, 2)-F(6, 16),
                cutoffs='c=chi_x(r^2,z)*chi_t(t), chi_x=chi_t=1 near (0,1)',
                checks=dict(support_stays_below_qstar_over_two=F(6, 16) < F(1, 2),
                            temporal_cutoff_zero_near_initial_time=F(1, 16) < F(1, 2),
                            exponent_for_factor_two_is_below_two=1/(1-2*F(1, 100)) < 2,
                            axisymmetric_cutoff_preserves_direct_swirl_divergence=True))


def borel_log2_scale(j, constraints, previous_power=1):
    """Choose b_j=2^n_j from finite (spatial order,time order,log constant).

    The log constant includes the cutoff derivative factors and the spatial
    norm of F_j. Each constraint satisfies a+m<=floor(j/2). Integers n_j
    keep very large scales representable without evaluating 2^n_j.
    """
    if type(j) is not int or j < 1 or type(previous_power) is not int or previous_power < 1 or not constraints:
        raise ValueError('Positive stage, previous dyadic exponent and constraints required')
    with mp.workdps(70):
        n = previous_power+1
        for a, m, logC in constraints:
            order(a), order(m)
            logC = mp.mpf(logC)
            if a+m > j//2 or not mp.isfinite(logC) or logC < 0:
                raise ValueError('Require finite padded constants and a+m<=floor(j/2)')
            # Strict unit padding avoids a rounded equality at the boundary.
            n = max(n, int(mp.ceil((logC+j*mp.log(2)+1)/((j-m)*mp.log(2)))))
        return n


def borel_scale_review():
    previous, rows = 1, []
    for j in (1, 2, 7, 20):
        constraints = [(a, m, (j+1)**4+10**(j//3))
                       for a in range(j//2+1) for m in range(j//2+1-a)]
        n = borel_log2_scale(j, constraints, previous)
        with mp.workdps(70):
            worst = max(mp.mpf(logC)+(m-j)*n*mp.log(2)+j*mp.log(2)
                        for a, m, logC in constraints)
        rows.append(dict(stage=j, dyadic_scale_power=n, finite_constraints=len(constraints),
                         all_constraint_log_margin=float(worst), scale_doubles=n >= previous+1))
        previous = n
    return dict(synthetic_constants=True, actual_terminal_jets_evaluated=False,
                term_bound='|D_x^alpha partial_sigma^m [chi0(b_j*sigma)*sigma^j/j!*F_j]| <= C_jm*b_j^(m-j)*||F_j||_C^a',
                diagonal='a+m<=floor(j/2); every fixed differentiated tail is uniformly summable',
                rows=rows,
                checks=dict(every_finite_constraint_met=all(r['all_constraint_log_margin'] < 0 for r in rows),
                            shrinking_support_scales_increase=all(r['scale_doubles'] for r in rows),
                            every_fixed_mixed_jet_is_eventually_in_the_diagonal=True,
                            no_convergence_assumption_on_raw_terminal_Taylor_series=True))


def energy_ledger():
    return dict(initial_data='u(0)=0', force_time_norm='F(t)=integral_0^t ||f(s)||_2 ds',
                energy_identity='0.5*d_t ||u||_2^2+nu*||grad u||_2^2=<f,u>',
                norm_estimate='||u(t)||_2<=F(t)',
                integrated_bound='||u(t)||_2^2+2*nu*integral_0^t ||grad u||_2^2<=F(t)^2',
                terminal_limit='monotone convergence in dissipation; only t<1 solution is used',
                checks=dict(zero_initial_data_and_compact_preterminal_fields_remove_boundary_terms=True,
                            compact_smooth_force_has_finite_F_at_one=True,
                            whole_space_energy_and_integrated_dissipation_follow_conditionally=True,
                            no_velocity_extension_or_regular_velocity_at_t_one_inferred=True,
                            no_comparison_with_an_arbitrary_global_solution_inferred=True))
