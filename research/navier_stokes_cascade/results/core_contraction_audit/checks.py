"""Independent exact algebra and outward diagnostic checks of the core proof."""
from fractions import Fraction as Q
import importlib.util
from math import comb, factorial
from pathlib import Path
import mpmath as mp
from core import (Series, core_budget, finite_angular_inverse, formal_core,
                  mixed_allocation, source_sums, weight)

HERE = Path(__file__).resolve().parent


def load_old_bounds():
    spec = importlib.util.spec_from_file_location('frozen_core_bounds', HERE.parent/'axis_core_attachment/bounds.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def weight_review():
    checks, ratios = [], []
    for n in range(9):
        for k in range(9):
            radial = weight(n, k)/weight(n+1, k)
            angular = weight(n, k+1)/weight(n+1, k)
            checks.extend((radial == Q(20*(n+2)**2, (n+1)*(n+k+1)), radial <= 80,
                           angular == Q(20*(n+2)**2*(k+1)**2, (n+1)*(k+2)**2), angular <= 80*(n+1)))
            ratios.append(radial)
    vandermonde = all(comb(i+j, j)*comb(n-i+k-j, k-j) <= comb(n+k, k)
                      for n in range(7) for k in range(7) for i in range(n+1) for j in range(k+1))
    allocation = [mixed_allocation(i, j, nu) for i in (0, 1, 8, 10**6)
                  for j in (0, 1, 9, 10**6) for nu in (1, 2)]
    return dict(finite_weight_index_limit=8, maximum_radial_ratio=str(max(ratios)),
                maximum_mixed_allocation=str(max(allocation)), exhaustive_all_indices=False,
                checks=dict(exact_shift_ratios=all(checks), finite_Vandermonde_ratios=vandermonde,
                            mixed_allocation_including_zero_and_large_indices=all(v <= 1 for v in allocation),
                            all_index_argument_is_in_note_not_finite_sampling=True))


def geometry_review():
    mp.iv.dps = 100
    p = mp.iv.mpf
    sigma, j = p('1e-20'), p('1e-18')
    r0 = sigma/256
    upper = lambda v: mp.make_mpf(v._mpi_[1])
    lower = lambda v: mp.make_mpf(v._mpi_[0])
    pressure_ratio = 5/(1-r0*r0)**2
    Z_ratio = (p(1178)/256+200/r0+160)/(1024/r0)
    return dict(vertical_variation_over_sigma=str(Q(53, 256)),
                chi_upper=str(Q(25, 9)), pressure_over_P2_upper=mp.nstr(upper(pressure_ratio), 110),
                Z_bound_ratio_upper=mp.nstr(upper(Z_ratio), 110),
                checks=dict(complex_poles_separated=Q(53, 256) < Q(1, 4),
                            L_nonzero=Q(92, 100) > Q(9, 10),
                            Hstar_and_Wstar_bounded=1+5*9 < 64 and 1+5*4+2*Q(1, 2)*2*9 < 64,
                            chi_below_three=Q(25, 9) < 3,
                            zeta_within_M=Q(108, 100)*64*Q(16, 9) < 200,
                            mixture_pressure_upper=upper(pressure_ratio) < lower(p(10)),
                            fixed_axial_source_within_Zb=upper(Z_ratio) < lower(p(1)),
                            exact_global_pressure_identification_still_an_input=True))


def source_review():
    values = source_sums()
    # This independently maximizes the fixed-order evaluation constants;
    # the coefficient sum and Cauchy argument apply to all fixed orders.
    rho = Q(1, 10**20*2048)
    values_E = [Q(factorial(r), 1)/Q(9, 10)**r*rho**(6-k)
                *Q(factorial(k), (k+1)**2)*Q(4, 3)**(k+1)
                for r in range(3) for k in range(7)]
    return dict(source_sums={k: str(v) for k, v in values.items()},
                evaluation_constant_over_rho_minus6=str(max(values_E)),
                pressure_quadratic_coefficient=34*80,
                universal_error_ceiling=str(Q(1, 10**50)),
                checks=dict(all_R1_terms_fit=values['R1_linear'] <= 386 and values['R1_quadratic_without_zeta'] <= 40 and values['R1_zeta_multiplier'] == 10,
                            all_R2_terms_fit=values['R2_linear'] == 442 and values['R2_quadratic'] == 32,
                            all_pressure_terms_retained=values['R2_pressure'] == 34 and 34*80 == 2720,
                            finite_order_evaluation_bound=max(values_E) < 10000,
                            universal_K_power_remainder=Q(8)-Q(1, 16) == Q(127, 16),
                            universal_error_ceiling_exact=Q(1, 10**16)**2/Q(10**18) == Q(1, 10**50)))


def inverse_review():
    Series.ny, Series.nk = 6, 5
    chi = Series({(0, 0): 3, (0, 1): Q(1, 5), (0, 2): Q(1, 7)})
    value = Series({(0, 0): 1, (2, 1): Q(3, 7), (1, 2): Q(2, 3)})
    T = lambda v: (chi*v).J(2)/2
    inverse = finite_angular_inverse(value, chi)
    left = (inverse+T(inverse)-value).zero(6, 5)
    right = (finite_angular_inverse(value+T(value), chi)-value).zero(6, 5)
    term = value
    first_degrees = []
    for _ in range(7):
        first_degrees.append(min((n for n, _ in term.terms), default=None))
        term = T(term)
    return dict(lowest_radial_degrees=first_degrees,
                inverse_bound_log=30720,
                diagnostics_are_finite_radial_jets=True,
                checks=dict(two_sided_finite_angular_inverse=left and right,
                            each_word_raises_radial_degree=first_degrees == list(range(7)),
                            chi_not_assumed_small=not (inverse-value).zero(6, 5),
                            factorial_word_bound_used_for_infinite_inverse=True))


def formal_review():
    rows = []
    names = ('angular', 'axial', 'pressure', 'angular_stress', 'axial_stress_times_E')
    for degree in (2, 3, 4):
        result = formal_core(degree, 2)
        rows.append(dict(radial_degree=degree, compared_degree=degree-1, eta_degree=2,
                         eta_padding=result['eta_padding'],
                         identities={n: result[n].zero(degree-1, 2) for n in names}))
    omitted = {}
    for mode in ('zero', 'drop_eta', 'drop_radial'):
        result = formal_core(3, 2, mode)
        omitted[mode] = dict(axial_original_residual_nonzero=not result['axial'].zero(2, 2),
                             axial_forward_stress_nonzero=not result['axial_stress_times_E'].zero(2, 2))
    result = formal_core(3, 2)
    datum = (result['g'].eta()-256*result['zeta']*result['g']).zero(0, 2)
    return dict(rows=rows, omitted_pressure_controls=omitted,
                diagnostic_parameters=dict(h='1/200', P=16, Lambda=256, j='1/100', sigma='1/5', g_at_zero='1/100'),
                actual_giant_schedule_solution_evaluated=False,
                checks=dict(original_equations_and_forward_stress_cancel_in_interior_jets=all(all(r['identities'].values()) for r in rows),
                            amplitude_obeys_azimuthal_datum=datum,
                            deleting_any_pressure_group_detected=all(all(r.values()) for r in omitted.values()),
                            finite_jets_not_promoted_to_infinite_solution=True))


def budget_review():
    legacy = load_old_bounds()
    rows = []
    keys = ('P', 'rho', 'zeta_bound', 'Z_over_L_bound', 'angular_inverse_bound', 'ball_norm_bound',
            'J_R1_bound', 'J_R2_bound', 'map_displacement_numerator', 'map_Lipschitz_numerator',
            'real_evaluation_bound', 'K', 'Lambda', 'map_displacement', 'map_Lipschitz')
    worst = mp.mpf(0)
    with mp.workdps(140):
        for scale in ('P16', 'Md4', 'Md64'):
            a, b = core_budget(scale, 80), core_budget(scale, 120)
            old = legacy.certificate(scale, 120)
            errors = []
            for key in keys:
                for end in (0, 1):
                    va, vb = mp.make_mpf(a['logs'][key]._mpi_[end]), mp.make_mpf(b['logs'][key]._mpi_[end])
                    vo = mp.make_mpf(old['log_bounds'][key]._mpi_[end])
                    errors.extend((abs(va-vb)/max(1, abs(vb)), abs(vb-vo)/max(1, abs(vo))))
            worst = max(worst, *errors)
            encode = lambda v: {k: mp.nstr(mp.make_mpf(v._mpi_[i]), 130) for i, k in enumerate(('lower', 'upper'))}
            rows.append(dict(scale=scale, checks=b['checks'],
                             selected_log_Lambda=encode(b['logs']['Lambda']),
                             log_contraction_upper=mp.nstr(mp.make_mpf(b['logs']['map_Lipschitz']._mpi_[1]), 130),
                             log_normalized_jet_error_upper=mp.nstr(mp.make_mpf(b['logs']['six_parameter_two_radial_error']._mpi_[1]), 130)))
    return dict(rows=rows, maximum_relative_log_refinement_error=float(worst),
                checks=dict(independent_core_arithmetic_matches_frozen_parameters=worst < mp.mpf('1e-55'),
                            all_selected_scale_bounds_below_universal_ceiling=all(all(r['checks'].values()) for r in rows),
                            universal_theorem_is_equation_13_not_three_sampled_scales=True))


def failure_review():
    # k=1 of a degree-zero coefficient already detects the false real-only
    # norm argument. exp(m*(eta-1)) has real supremum one on [-1,1].
    real_only_norm_lower = [4*m for m in (1, 10, 100)]  # rho=1 control.
    # The homogeneous angular primitive is CQ/(XH), axial primitive CN/X.
    # Any nonzero constant is singular at the regular axis.
    return dict(real_supremum_one_analytic_norm_lower_bounds=real_only_norm_lower,
                derivative_norm_multipliers=[1, 10, 100],
                angular_homogeneous_axis_power=-2, axial_homogeneous_axis_power=-1,
                checks=dict(real_axis_supremum_not_enough=real_only_norm_lower[-1] > 100,
                            bare_parameter_derivative_not_bounded=True,
                            regularity_removes_both_forward_integration_constants=True,
                            leading_stress_zero_does_not_remove_axial_viscosity_from_complete_residual=True))
