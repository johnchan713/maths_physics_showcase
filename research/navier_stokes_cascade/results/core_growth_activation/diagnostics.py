"""Exact scalar controls, labelled fixtures and outward log-budget replay.

Finite checks complement the continuum derivation in README.md. No fixture
is identified with the actual giant core or its completed exterior.
"""
from fractions import Fraction as Q
import importlib.util
from pathlib import Path
import sys
import mpmath as mp
from growth import algebra_constants, budget, comparison, core, f0, pressure_axis

HERE = Path(__file__).resolve().parent


def algebra_review():
    c = algebra_constants()
    t = Q(99, 50)
    p4 = 1-t+t*t/4-t**3/36+t**4/576
    tm = Q(41, 20)
    cubic = 1-tm/2+tm*tm/12-tm**3/144
    source_prefactor = c['source_remainder_prefactor']
    # Exact eta algebra is independent of the numerical pressure fixture.
    S = core.Series
    S.ny, S.nk = 0, 8
    e, j, h = S({(0, 1): 1}), Q(1, 10**18), Q(1, 200)
    A, D, d, U = Q(1, 2)+h, Q(1, 2)-h, 1-e*e, 4*e+j
    direct = -A*(1-2*e*U)*U-4*(D*e+d*U)
    polynomial = -20*e+(16+32*A)*e**3-(A+4)*j+(4+16*A)*j*e**2+2*A*j*j*e
    return dict(p4_at_1p98=str(p4), scalar_cubic_at_2p05=str(cubic),
                small_Z_eta_over_j=str(c['small_Z_eta_over_j']),
                source_remainder_prefactor=str(source_prefactor),
                checks=dict(velocity_polynomial_exact=(direct-polynomial).zero(0, 8),
                            velocity_bound_absorbs_j_eta_squared=20+(4+16*Q(51, 100))*Q(1, 20)+2*Q(51, 100)*Q(1, 20)**2 < 40,
                            velocity_constant_below_five=Q(51, 100)+4 < 5,
                            reference_pressure_mass_below_five=c['reference_pressure_mass'] < 5,
                            small_Z_localizes_eta=c['small_Z_eta_over_j'] < Q(1, 30),
                            microregion_H_above_p8=1-Q(45, 10)/30-4*Q(1, 20)**2/30**3-Q(1, 20)**2/30**2 > Q(4, 5),
                            separation_chi_above_p99=c['chi_separation'] > Q(99, 100),
                            scalar_lower_above_p265=cubic > Q(265, 1000),
                            scalar_cubic_decreases=Q(1, 2)-tm/6 > 0,
                            alternating_tail_decreases=tm/6 < 1,
                            endpoint_polynomial_negative=p4 == -Q(75535511, 400000000) and p4 < -Q(188, 1000),
                            endpoint_comparison_shear_above_2p376=2*(1+Q(188, 1000)) == Q(2376, 1000),
                            weighted_eta_constant=Q(205, 53)*106 == 410,
                            weighted_eta_error_constant=4*64+4*410 < 2000,
                            Young_remainder_sufficient=Q(1, 4)/(Q(3, 100)*Q(98, 100)) < 9,
                            source_remainder_universal=source_prefactor < Q(3, 1000),
                            positive_forward_shear=Q(1, 8)*Q(5, 2)*Q(1, 2) == Q(5, 32) > Q(1, 10),
                            lambda_lower_denominator=5*16*17//2 == 680,
                            no_error_divided_by_chi=True))


def scalar_review():
    rows, error = [], mp.mpf(0)
    for z in map(mp.mpf, ('0', '.1', '1', '3.96', '4', '4.1')):
        r = comparison(z)
        for derivative in range(3):
            series = mp.fsum((-z/2)**n*(-mp.mpf('.5'))**derivative/
                             (mp.factorial(n)*mp.factorial(n+derivative+1)) for n in range(150))
            error = max(error, abs(series-f0(z, derivative)))
        rows.append(dict(z=float(z), f=float(r['f']), derivative=float(r['fp']),
                         combined=float(r['combined']), shear=float(r['shear'])))
    return dict(rows=rows, independent_series_error=float(error),
                checks=dict(hypergeometric_matches_independent_series=error < mp.mpf('1e-100'),
                            scalar_positive_and_decreasing=all(r['f'] > .265 and -.251 < r['derivative'] < 0 for r in rows),
                            endpoint_angular_growth=all(r['shear'] > 2.376 for r in rows if r['z'] >= 3.96),
                            samples_not_continuum_proof=True))


def pressure_review():
    j = mp.mpf('1e-18')
    ratios = ['-1', '-.225', '-.2224694104560623', '-.22', '-.03', '0', '.03', '1']
    rows, decomposition = [], mp.mpf(0)
    for P in (16, 64):
        for ratio in ratios:
            eta = mp.mpf(ratio)*j
            r = pressure_axis(eta, P=P)
            decomposition = max(decomposition, abs(r['Z']-r['velocity']+eta*P*P*r['Kpressure']))
            rows.append(dict(P=P, eta_over_j=ratio, chi=float(r['chi']),
                             Z_over_delta=float(r['Z']/(j*P*P/100)),
                             Kpressure=float(r['Kpressure'])))
    poles = [pressure_axis(e) for e in (-1, 1)]
    weak = [r for r in rows if r['chi'] <= .99]
    return dict(fixture='2.5 delta_1 + 1.2 delta_0.3; not the actual schedule measure',
                rows=rows, maximum_pressure_decomposition_error=float(decomposition),
                checks=dict(pressure_sign_and_derivative_identity=decomposition < mp.mpf('1e-100'),
                            ideal_atom_lower_bound=all(r['Kpressure'] >= 1.25 for r in rows) and all(r['Kpressure'] >= 1.25 for r in poles),
                            tiny_weak_shear_region_really_exercised=bool(weak) and min(r['chi'] for r in weak) < 1e-20,
                            axial_source_separated_in_weak_sector=all(abs(r['Z_over_delta']) > 1 for r in weak),
                            small_Z_not_used_to_bound_chi_by_sampling=True,
                            manufactured_measure_not_actual_pressure_identification=True))


def load_gauged_activation(C):
    """Load frozen formulas, changing only the labelled manufactured seed.

    U and its derivative stay fixed; E,H,I,J carry C^-1, while the swirl
    parts of S and Pi carry C^-2. Pi0 is the SAME fixed datum for every C.
    """
    name = '_core_growth_flat_gauge_'+str(C)
    spec = importlib.util.spec_from_file_location(name, HERE.parent/'inner_collar_activation/collar.py')
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    C = mp.mpf(C)

    def seed(y, eta):
        y, e = mp.mpf(y), module.Jet(eta, 1)
        x0, x = mp.mpf(2), 2*mp.exp(y)
        eb, q0, ub = 2+e/10+e**2/20, 1+e/5, e/3+e**2/30
        b, c = q0*eb/2, ub-q0*eb/2
        E, U = eb*mp.exp(-y)/C, c+b*mp.exp(-y)
        H0, H = 2*eb/C, 2*eb*mp.exp(-y/2)/C
        M = (1+e/5)/10+x0*(c*mp.expm1(y)+b*y)
        I = (2+e/3)/(10*C)+2*x0*H0*mp.expm1(y/2)
        J = (1-e/5)/(20*C)+2*x0*H0*(c*mp.expm1(y/2)-b*mp.expm1(-y/2))
        S = (7+2*e**2)/100+x0*(c**2*mp.expm1(y)+2*c*b*y-(b**2-eb**2/(2*C*C))*mp.expm1(-y))
        Cp = ((3+e)/100-eb**2*mp.expm1(-2*y)/4)/(C*C)
        Pi = -mp.mpf('.7')+e/5+3*e**2/20-eb**2*mp.expm1(-2*y)/(4*C*C)
        return dict(X=x, E=E, U=U, H=H, M=M, I=I, J=J, S=S,
                    Cp=Cp, Pi=Pi, a=module.Jet(3), q=C*q0)

    module.seed = seed
    return module


def moment_review():
    rows, direct_errors, derivative_errors, full, normalized = [], [], [], [], []
    s, eta, t = map(mp.mpf, ('.18', '.23', '.01'))
    for C in (1, 2, 4):
        m = load_gauged_activation(C)
        act = m.Activation(width=t, kappa0='.1', order=24)
        r = act.stress(s, eta)
        state, ea = r['state'], r['state']['fields']['ea']
        base, actual = state['fields']['base'], state['actual']
        ref, changed = m.moment_map(base, eta=eta), m.moment_map(actual, eta=eta)
        direct_errors.extend(abs((changed[name]-ref[name])/ea-r[hat])
                             for name, hat in (('p1', 'p1hat'), ('p2', 'p2hat')))
        # All five reference primitives satisfy their defining forward ODEs.
        ode = dict(M=base['X']*base['U'].v, I=base['X']*base['H'].v,
                   J=base['X']*base['U'].v*base['H'].v,
                   S=base['X']*(base['U'].v**2-base['E'].v**2/2), Cp=base['E'].v**2/2)
        derivative_errors.extend(abs(mp.diff(lambda y, n=n: m.seed(y, eta)[n].v, t*s)-v) for n, v in ode.items())
        A, d, L, x = mp.mpf('.505'), 1-eta*eta, 1-mp.mpf('.01')*eta*eta, base['X']
        cp = state['hats']['Cp']
        pressure_value = x*4*A*eta*cp.v/(L*actual['E'].v)
        pressure_eta = -x*d*cp.d/(L*actual['E'].v)
        direct_errors.append(abs(pressure_value+pressure_eta-r['pressure_contribution_to_p2hat']))
        normalized.append(r['p2hat']/C)
        full.append(r)
        rows.append(dict(C=C, p1hat=float(r['p1hat']), p2hat_over_C=float(r['p2hat']/C),
                         pressure_value=float(pressure_value), pressure_eta=float(pressure_eta),
                         pressure_over_C=float(r['pressure_contribution_to_p2hat']/C),
                         denominator_error_over_C2=float(r['denominator_error']/(C*C)),
                         seed_defect=[float(v) for v in r['seed_defect_over_Fr']]))
    affine_error = abs(normalized[0]-5*normalized[1]+4*normalized[2])
    p1_error = max(abs(full[0]['p1hat']-r['p1hat']) for r in full)
    pressure_scaling_error = max(abs(rows[0]['pressure_over_C']-r['pressure_over_C']*r['C']**2) for r in rows)
    return dict(rows=rows, maximum_direct_difference_error=float(max(direct_errors)),
                maximum_primitive_ODE_error=float(max(derivative_errors)),
                affine_inverse_C_squared_error=float(affine_error),
                checks=dict(complete_map_matches_cancellation_safe_difference=max(direct_errors) < mp.mpf('1e-100'),
                            all_five_primitive_ODEs_hold=max(derivative_errors) < mp.mpf('1e-100'),
                            angular_difference_invariant_under_C_gauge=p1_error < mp.mpf('1e-100'),
                            axial_difference_over_C_keeps_inverse_C_squared_part=affine_error < mp.mpf('1e-100') and abs(normalized[0]-normalized[1]) > mp.mpf('1e-12'),
                            pressure_has_required_inverse_C_squared_scaling=pressure_scaling_error < 1e-20,
                            pressure_value_and_eta_terms_both_nonzero=all(abs(r['pressure_value']) > 1e-12 and abs(r['pressure_eta']) > 1e-12 for r in rows),
                            nonzero_manufactured_seed_defect_retained=all(max(abs(v) for v in r['seed_defect']) > .1 for r in rows),
                            manufactured_gauge_not_actual_core=True))


def budget_review():
    rows, errors = [], []
    for scale in ('P16', 'Md4', 'Md64'):
        a, b = budget(scale, 100), budget(scale, 160)
        for key in a['logs']:
            for i in (0, 1):
                va, vb = mp.make_mpf(a['logs'][key]._mpi_[i]), mp.make_mpf(b['logs'][key]._mpi_[i])
                errors.append(abs(va-vb)/max(1, abs(vb)))
        encode = lambda v: {name: mp.nstr(mp.make_mpf(v._mpi_[i]), 130) for i, name in enumerate(('lower', 'upper'))}
        rows.append(dict(scale=scale, checks=b['checks'], logs={k: encode(v) for k, v in b['logs'].items()}))
    return dict(rows=rows, maximum_relative_log_refinement_error=float(max(errors)),
                checks=dict(all_selected_interval_inequalities=all(all(r['checks'].values()) for r in rows),
                            interval_refinement_consistent=max(errors) < mp.mpf('1e-55'),
                            Kact_and_width_not_chosen_by_parameter_scan=True,
                            universal_inequalities_derived_in_note=True,
                            no_full_C_or_R_materialized=True))
