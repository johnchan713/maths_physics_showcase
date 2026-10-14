"""Independent manufactured checks of the compact mean identities.

The continuum arguments are in README.md. Quadrature and finite Fourier
fixtures test the formulas and deliberate omissions; they do not evaluate
the actual global leading profile, summed flow, or its auxiliary geometry.
"""
from fractions import Fraction
import math
import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre
from operators import (bump_moment, divisor, mean_update_budget, moment_blocks,
                       polynomial_bump, polynomial_derivative,
                       polynomial_oscillatory_integral, polynomial_value,
                       radial_mode, solve_moment_targets, temporal_inverse,
                       torus_vectors)
from mean_system import (Jet, arrays, bump, mean_balances, mean_correction,
                         mean_increment, radial_increment_remainder, scaled_data,
                         variables, waves)


def relative(a, b):
    return float(np.linalg.norm(np.asarray(a)-np.asarray(b))/(1+np.linalg.norm(a)+np.linalg.norm(b)))


def physical_mean_review():
    """Full angular momentum residual versus conservative mean equations."""
    maximum = normalized_error = remainder_error = divergence_error = 0.
    omission = dict(axial_cross_covariance=0., axial_normal_covariance=0.,
                    radial_axial_covariance=0., viscosity=0., quadratic_increment=0.,
                    radial_background_memory=0.)
    for point in ((1.12, .2, -.19, .36), (1.3, .27, .2, .47), (1.66, -.9, .31, .64)):
        coords = variables(point)
        r = coords[0]
        base, pbase = waves.background(coords)
        old, increment = mean_correction(coords), mean_increment(coords)
        new = [a+b for a, b in zip(old, increment)]
        q, _ = waves.covariance(coords, 11)
        pm = .17*bump((r-1.4)/.5)*(1+.2*coords[2].sin()+.1*coords[3].cos())
        gr, et, ez = mean_balances(coords, base, old, q, pm)
        analytic = np.array([(pm.diff(0)-gr).value, et.value, ez.value])
        wdata, wave_pg = waves.harmonic_samples(point, 11)
        mean_data = arrays([a+b for a, b in zip(base, old)])
        full_data = tuple(a+b for a, b in zip(mean_data, wdata))
        total_pg = waves.pressure_gradient(pbase+pm, point[0])+wave_pg
        direct = np.mean(waves.momentum(full_data, total_pg, point[0]), axis=0)
        direct -= waves.momentum(arrays(base), waves.pressure_gradient(pbase, point[0]), point[0])
        maximum = max(maximum, relative(analytic, direct))
        divergence_error = max(divergence_error, abs(waves.divergence(arrays(old), point[0])))
        for key, indices in (('axial_cross_covariance', (1, 2)),
                             ('axial_normal_covariance', (2, 2)),
                             ('radial_axial_covariance', (0, 2))):
            bad = [list(row) for row in q]
            i, j = indices
            bad[i][j] = bad[j][i] = Jet(0)
            bg, bt, bz = mean_balances(coords, base, old, bad, pm)
            wrong = np.array([(pm.diff(0)-bg).value, bt.value, bz.value])
            omission[key] = max(omission[key], float(np.linalg.norm(wrong-analytic)))
        gr_new = mean_balances(coords, base, new, q)[0]
        rg = radial_increment_remainder(coords, base, old, increment)
        actual_change = (gr_new-gr-2*base[1]*increment[1]/r).value
        remainder_error = max(remainder_error, abs(actual_change-rg.value))
        omission['quadratic_increment'] = max(omission['quadratic_increment'],
                                               abs((increment[1]*increment[1]/r).value))
        memory_flux = 2*base[0]*increment[0]
        memory = memory_flux.diff(0)+memory_flux/r+(base[0]*increment[2]).diff(2)
        omission['radial_background_memory'] = max(omission['radial_background_memory'], abs(memory.value))
        # A separate normalized fixture has exact normalized incompressibility.
        for eps in (.07, .19):
            rr, _, zz, tt = coords
            bnorm = [-eps*.17*rr, base[1], .34*zz]
            mn = mean_correction(coords, eps)
            zero = [[Jet(0) for _ in range(3)] for _ in range(3)]
            ng, nt, nz = mean_balances(coords, bnorm, mn, zero, pm, epsilon=eps, normalized=True)
            ndata = scaled_data(arrays([a+b for a, b in zip(bnorm, mn)]), eps)
            bdata = scaled_data(arrays(bnorm), eps)
            pg = waves.pressure_gradient(pm, point[0])*np.array([1., 1., eps])
            direct_norm = (ndata[1][:, 3]+waves.advection(ndata, ndata, point[0])
                           -eps*waves.laplacian(ndata, point[0])+pg
                           -bdata[1][:, 3]-waves.advection(bdata, bdata, point[0])
                           +eps*waves.laplacian(bdata, point[0]))
            expected = np.array([(pm.diff(0)-ng).value, nt.value, nz.value])
            normalized_error = max(normalized_error, relative(direct_norm, expected))
            omission['viscosity'] = max(omission['viscosity'],
                                        float(np.linalg.norm(eps*waves.laplacian(scaled_data(arrays(mn), eps), point[0]))))
    return dict(maximum_physical_angular_mean_relative_error=maximum,
                maximum_normalized_full_mean_relative_error=normalized_error,
                maximum_exact_radial_increment_error=remainder_error,
                maximum_curl_mean_divergence=divergence_error, omission_gaps=omission,
                checks=dict(full_six_entry_covariance_mean_identity=maximum < 2e-10,
                            normalized_time_axial_and_viscous_factors=normalized_error < 2e-10,
                            exact_increment_retains_radial_memory_and_quadratic_terms=remainder_error < 2e-10,
                            mean_curl_is_divergence_free=divergence_error < 2e-11,
                            all_selected_omissions_fail=min(omission.values()) > 1e-5))


def collocation_integral(nodes, left=.8, right=2.1):
    """Gauss integration and an independently integrated Legendre interpolant."""
    x, w = roots_legendre(nodes)
    table = np.polynomial.legendre.legvander(x, nodes)
    antiderivatives = np.empty((nodes, nodes))
    antiderivatives[:, 0] = x+1
    for j in range(1, nodes):
        antiderivatives[:, j] = (table[:, j+1]-table[:, j-1])/(2*j+1)
    projection = ((2*np.arange(nodes)+1)/2)[:, None]*table[:, :nodes].T*w[None, :]
    half = (right-left)/2
    return (left+right)/2+half*x, half*w, half*antiderivatives@projection


def piecewise_integral(nodes, edges):
    """Split at every compact-support edge before integrating steep jets."""
    parts = [collocation_integral(nodes, a, b) for a, b in zip(edges[:-1], edges[1:])]
    size = nodes*len(parts)
    primitive = np.zeros((size, size))
    all_weights = np.concatenate([part[1] for part in parts])
    for index, (_, _, local) in enumerate(parts):
        a, b = index*nodes, (index+1)*nodes
        primitive[a:b, :a] = all_weights[:a]
        primitive[a:b, a:b] = local
    return np.concatenate([part[0] for part in parts]), all_weights, primitive


_BUMP_MASS = quad(lambda s: math.exp(-1/(1-s*s)) if abs(s) < 1 else 0., -1, 1,
                  epsabs=1e-14, epsrel=1e-13)[0]


def pressure_bump(coords):
    r, _, z, time = coords
    scale = 1+.08*z+.03*time
    return bump((r/scale-1.4)/.12)/(scale*.12*_BUMP_MASS)


def integrated_mean_review(nodes=96, step=.002):
    """Pressure and integrated axial identities on full, compact source data.

    Only the pressure-source Z derivative uses a five-point diagnostic
    stencil: gr already uses third derivatives of a potential. Asking the
    inherited third-order Taylor object for one more derivative would be
    wrong. A tighter stencil/quadrature is checked separately in tests.
    """
    eps, z, time = .13, .2, .47
    scale = 1+.08*z+.03*time
    edges = [.8, .85, .92, scale*(1.4-.12), scale*(1.4+.12), 1.88, 1.95, 2.1]
    rvalues, weights, integral = piecewise_integral(nodes, edges)
    rows = []
    for radius in rvalues:
        coords = variables((radius, 0., z, time))
        r, _, zz, _ = coords
        base = [-eps*.17*r, waves.background(coords)[0][1], .34*zz]
        old = mean_correction(coords, eps)
        q, _ = waves.covariance(coords, 11)
        gr, et, ez = mean_balances(coords, base, old, q, epsilon=eps, normalized=True)
        shifted = []
        for dz in (-2*step, -step, step, 2*step):
            cc = variables((radius, 0., z+dz, time))
            bb = [-eps*.17*cc[0], waves.background(cc)[0][1], .34*cc[2]]
            mm = mean_correction(cc, eps)
            qq, _ = waves.covariance(cc, 11)
            shifted.append(mean_balances(cc, bb, mm, qq, epsilon=eps, normalized=True)[0].value)
        grz = (shifted[0]-8*shifted[1]+8*shifted[2]-shifted[3])/(12*step)
        rho = pressure_bump(coords)
        v, gamma = old[1:]
        flux_t = base[2]*v+base[1]*gamma+gamma*v+q[1][2]
        flux_z = 2*base[2]*gamma+gamma*gamma+q[2][2]
        arbitrary = bump((r-1.4)/.48)*(coords[3]+zz*zz)
        bad_mean = [old[0], old[1]+arbitrary, old[2]]
        bad_et = mean_balances(coords, base, bad_mean, q, epsilon=eps, normalized=True)[1]
        extra_flux = (base[2]+gamma)*arbitrary
        rows.append([gr.value, grz, et.value, ez.value, rho.value, rho.diff(2).value,
                     v.value, gamma.value, flux_t.value, flux_t.diff(2).value,
                     flux_z.value, flux_z.diff(2).value, (bad_et-et).value, extra_flux.diff(2).value])
    data = np.array(rows)
    gr, grz, et, ez0, rho, rhoz, v, gamma, ft, ftz, fz, fzz, bad_delta, bad_flux_z = data.T
    p, pz = weights@gr, weights@grz
    pressure = integral@(gr-rho*p)
    pressure_z = integral@(grz-rhoz*p-rho*pz)
    c = .5*weights@(rvalues**2*rho)
    cz = .5*weights@(rvalues**2*rhoz)
    jtheta = weights@(rvalues**2*ft)
    jtheta_z = weights@(rvalues**2*ftz)
    jz = weights@(rvalues*fz)-.5*weights@(rvalues**2*gr)
    jz_z = weights@(rvalues*fzz)-.5*weights@(rvalues**2*grz)
    lhs = np.array([weights@(rvalues**2*et), weights@(rvalues*(ez0+eps*pressure_z))])
    rhs = eps*np.array([jtheta_z, jz_z+cz*p+c*pz])
    pressure_moment = weights@(rvalues*pressure)
    pressure_moment_expected = -.5*weights@(rvalues**2*gr)+c*p
    # A zero-constraint control: an arbitrary extra azimuthal bump changes
    # its time/axial moment terms even though it stays divergence-free.
    arbitrary_moment = weights@(rvalues**2*np.array([
        bump((Jet(radius)-1.4)/.48).value for radius in rvalues]))
    zero_moment_omission_gap = eps*(1+2*eps**2)*arbitrary_moment
    observed_nonzero_moment_change = weights@(rvalues**2*(bad_delta-eps*bad_flux_z))
    off_constraint_error = abs(observed_nonzero_moment_change+zero_moment_omission_gap)
    return dict(nodes_per_segment=nodes, total_nodes=len(rvalues), radial_segment_edges=edges,
                pressure_z_stencil_step=step,
                radial_defect=float(p), angular_flux_defect=float(jtheta), axial_flux_defect=float(jz),
                c_rho=float(c), c_rho_z=float(cz), pressure_bump_mass=float(weights@rho),
                pressure_source_total=float(weights@(gr-rho*p)),
                pressure_moment=float(pressure_moment), expected_pressure_moment=float(pressure_moment_expected),
                integrated_mean_residual=lhs.tolist(), expected_integrated_flux_derivatives=rhs.tolist(),
                maximum_integrated_identity_error=float(np.max(np.abs(lhs-rhs))),
                maximum_zero_velocity_moment=float(max(abs(weights@(rvalues**2*v)), abs(weights@(rvalues*gamma)))),
                omitted_moving_c_rho_derivative_gap=float(abs(eps*cz*p)),
                omitted_axial_pressure_moment_gap=float(abs(eps*(c*pz+cz*p))),
                nonzero_velocity_moment_control_gap=float(abs(zero_moment_omission_gap)),
                nonzero_velocity_moment_observed_change=float(observed_nonzero_moment_change),
                nonzero_velocity_moment_identity_error=float(off_constraint_error),
                checks=dict(compact_pressure_source_has_zero_integral=abs(weights@(gr-rho*p)) < 1e-8,
                            pressure_moment_keeps_radial_source=abs(pressure_moment-pressure_moment_expected) < 1e-8,
                            complete_integrated_mean_identities=max(abs(lhs-rhs)) < 1e-8,
                            two_velocity_moments_are_preserved=max(abs(weights@(rvalues**2*v)), abs(weights@(rvalues*gamma))) < 1e-8,
                            moving_pressure_bump_derivative_required=abs(eps*cz*p) > 1e-5,
                            axial_pressure_moment_required=abs(eps*(c*pz+cz*p)) > 1e-5,
                            arbitrary_nonzero_velocity_moment_cannot_be_dropped=zero_moment_omission_gap > 1e-5,
                            off_constraint_time_and_axial_viscosity_identity=off_constraint_error < 1e-8))


def defect_increment_review():
    """Recompute every defect before/after an arbitrary compact mean curl.

    The linear rows and nonlinear remainders are integrated independently.
    Substituting the inverse's target rows then yields (22) in README.md.
    """
    rv, weights, _ = piecewise_integral(64, [.8, .85, .92, .99, 1.77, 1.88, 1.95, 2.1])
    eps, rows = .13, []
    for radius in rv:
        coords = variables((radius, 0., .2, .47))
        r, _, z, _ = coords
        base = [-eps*.17*r, waves.background(coords)[0][1], .34*z]
        old, delta = mean_correction(coords, eps), mean_increment(coords, eps)
        new = [a+b for a, b in zip(old, delta)]
        w, _ = waves.covariance(coords, 11)
        before = mean_balances(coords, base, old, w, epsilon=eps, normalized=True)[0]
        after = mean_balances(coords, base, new, w, epsilon=eps, normalized=True)[0]
        rg = radial_increment_remainder(coords, base, old, delta, epsilon=eps, normalized=True)
        b, vbase, gbase = base
        beta, v, gamma = old
        db, dv, dg = delta
        def jt(mean):
            return r*r*(gbase*mean[1]+vbase*mean[2]+mean[1]*mean[2]+w[1][2])
        def jz(mean, source):
            return r*(2*gbase*mean[2]+mean[2]*mean[2]+w[2][2])-r*r*source/2
        linear = [2*vbase*dv/r, r*r*(gbase*dv+vbase*dg), 2*r*gbase*dg-r*vbase*dv]
        remainder = [rg, r*r*(gamma*dv+v*dg+dg*dv), r*(2*gamma*dg+dg*dg)-r*r*rg/2]
        change = [after-before, jt(new)-jt(old), jz(new, after)-jz(old, before)]
        quadratic = [dv*dv/r, r*r*dg*dv, r*dg*dg-r*dv*dv/2]
        rows.append([x.value for x in change+linear+remainder+quadratic])
    data = np.array(rows)
    change, linear, remainder, quadratic = [weights@data[:, a:a+3] for a in (0, 3, 6, 9)]
    error = float(np.max(abs(change-linear-remainder)))
    return dict(integrated_actual_defect_changes=change.tolist(),
                integrated_base_linear_rows=linear.tolist(),
                integrated_nonlinear_and_radial_remainders=remainder.tolist(),
                quadratic_increment_contributions=quadratic.tolist(),
                maximum_recomputed_defect_increment_error=error,
                checks=dict(same_full_covariance_before_and_after_mean_update=error < 2e-10,
                            all_three_exact_defect_increment_identities=error < 2e-10,
                            nonlinear_defects_do_not_vanish_after_a_linear_solve=np.linalg.norm(quadratic) > 1e-4))


def torus_and_radial_review():
    vr, vt = torus_vectors()
    matrix = mp.matrix([[3, 1], [1, 5]])
    eigen_error = max(mp.norm(matrix*mp.matrix(vr)-(4-mp.sqrt(2))*mp.matrix(vr)),
                      mp.norm(matrix*mp.matrix(vt)-(4+mp.sqrt(2))*mp.matrix(vt)))
    lower_margin = mp.inf
    for k1 in range(-12, 13):
        for k2 in range(-12, 13):
            if k1 or k2:
                for temporal in (False, True):
                    lower_margin = min(lower_margin, abs(divisor((k1, k2), temporal))
                                       *(1+mp.sqrt(2))*(abs(k1)+abs(k2)))
    modes = {(1, -2): mp.mpc('.3', '.1'), (-1, 2): mp.mpc('.3', '-.1'),
             (2, 1): mp.mpc('-.2', '.4'), (-2, -1): mp.mpc('-.2', '-.4')}
    solved = temporal_inverse(modes, '.007')
    temporal_error = max(abs(2j*mp.pi*mp.mpf('.007')*divisor(k, True)*v-modes[k]) for k, v in solved.items())
    reality_error = max(abs(solved[(-k[0], -k[1])]-mp.conj(v)) for k, v in solved.items())
    zero_rejected = False
    try:
        temporal_inverse({(0, 0): 1}, '.007')
    except ValueError:
        zero_rejected = True
    chi = lambda r: (3*((r*r-1)/3)**2-2*((r*r-1)/3)**3)
    coefficients = polynomial_bump()
    radial_error = support_error = omitted_path_gap = mp.mpf(0)
    derivative_rows, flatness_rows = [], []
    mode = (1, 2)
    for scale in (1, 8, 64, 512):
        omega = 2*mp.pi*scale*divisor(mode)
        full = polynomial_oscillatory_integral(coefficients, 1, 4, omega)
        for e in (0, 1, 2):
            for rr in map(mp.mpf, ('1.17', '1.43', '1.81')):
                value = radial_mode(rr, mode, scale, e, chi, coefficients)
                ordinary_derivative = mp.diff(lambda r: radial_mode(r, mode, scale, e, chi, coefficients), rr)
                graph_derivative = ordinary_derivative+1j*omega*2*rr*value+e*value/rr
                source = rr**(-e)*2*rr*polynomial_value(coefficients, rr*rr)
                cutoff_error = rr**(-e)*mp.diff(chi, rr)*mp.exp(-1j*omega*rr*rr)*full
                radial_error = max(radial_error, abs(graph_derivative-source+cutoff_error))
                omitted_path_gap = max(omitted_path_gap, abs(1j*omega*2*rr*value))
        # Here chi has plateau values but only finite smoothness at endpoints;
        # this diagnostic tests values/first derivatives, not C-infinity.
        plateau_chi = lambda r: 0 if r <= 1 else (1 if r >= 2 else chi(r))
        for rr in ('.9', '2.1'):
            support_error = max(support_error, abs(radial_mode(rr, mode, scale, 1, plateau_chi, coefficients)))
        fixed_derivative = abs(mp.diff(lambda r: radial_mode(r, mode, scale, 0, chi, coefficients), mp.mpf('1.43')))
        derivative_rows.append(dict(scale=scale, coefficient_first_derivative=fixed_derivative,
                                    full_radial_integral_modulus=abs(full)))
        for p in (1, 2, 4):
            dp = polynomial_derivative(coefficients, p)
            integrated_dp = polynomial_oscillatory_integral(dp, 1, 4, omega)
            by_parts = (-1/(1j*omega))**p*integrated_dp
            flatness_rows.append(dict(scale=scale, integrations=p, identity_error=abs(full-by_parts)))
    # The zero Fourier mode is an exact obstruction unless its weighted
    # radial integral vanishes: this source's integral was normalized to 1.
    zero_mode_error = mp.diff(chi, mp.mpf('1.43'))
    return dict(maximum_torus_eigenvector_error=eigen_error,
                sampled_diophantine_lower_margin=lower_margin,
                maximum_temporal_inverse_error=temporal_error, maximum_reality_error=reality_error,
                maximum_weighted_graph_primitive_error=radial_error, outside_shell_value_error=support_error,
                omitted_radial_phase_path_gap=omitted_path_gap,
                nonzero_zero_mode_cutoff_remainder=zero_mode_error,
                derivative_rows=derivative_rows, flatness_rows=flatness_rows,
                checks=dict(actual_auxiliary_eigendirections=eigen_error < mp.mpf('1e-65'),
                            algebraic_divisor_bound_crosscheck=lower_margin >= 1,
                            zero_mean_temporal_inverse_and_reality=max(temporal_error, reality_error) < mp.mpf('1e-65'),
                            nonzero_auxiliary_mean_is_rejected=zero_rejected,
                            all_three_weighted_radial_graph_identities=radial_error < mp.mpf('1e-55'),
                            radial_path_cannot_be_omitted=omitted_path_gap > 1,
                            primitive_values_keep_shell_support=support_error == 0,
                            integration_by_parts_identity=max(row['identity_error'] for row in flatness_rows) < mp.mpf('1e-55'),
                            zero_radial_mean_condition_is_necessary=abs(zero_mode_error) > mp.mpf('.1')))


def moment_map_review():
    rows, worst, determinant_error = [], mp.mpf(0), mp.mpf(0)
    for lam in map(mp.mpf, ('.2', '.05', '.001')):
        for amplitude in map(mp.mpf, ('1.3', '-.8')):
            at, az = moment_blocks(lam, amplitude)
            p = [mp.mpf(2), -2-2*lam, -2*lam]
            z = [mp.mpf(1), 1-2*lam]
            rt = [mp.exp(mp.mpf('.12')*v) for v in p]
            rz = [mp.exp(mp.mpf('.12')*v) for v in z]
            det_t = (-2*amplitude**2*mp.fprod(bump_moment(mp.nstr(v, mp.mp.dps)) for v in p)
                     *mp.fprod(rt[j]-rt[i] for i in range(3) for j in range(i+1, 3)))
            det_z = amplitude*mp.fprod(bump_moment(mp.nstr(v, mp.mp.dps)) for v in z)*(rz[1]-rz[0])
            determinant_error = max(determinant_error, abs(mp.det(at)-det_t), abs(mp.det(az)-det_z))
            inverse_norm = max(mp.norm(at**-1, p=mp.inf), mp.norm(az**-1, p=mp.inf))
            for targets in ((1, -2, 3), (-3, 1, -1), (0, 4, 0)):
                pp, jt, jz = targets
                u, s = solve_moment_targets(pp, jt, jz, lam, amplitude)
                worst = max(worst, mp.norm(at*u-mp.matrix([0, -pp, -jz])), mp.norm(az*s-mp.matrix([0, -jt])))
            rows.append(dict(fixed_lambda=lam, patch_amplitude=amplitude,
                             angular_determinant=mp.det(at), axial_determinant=mp.det(az),
                             maximum_inverse_entry=inverse_norm))
    degeneracy_rejected = False
    try:
        moment_blocks(0, 1)
    except ValueError:
        degeneracy_rejected = True
    budget = mean_update_budget()
    return dict(rows=rows, maximum_five_equation_solve_error=worst,
                maximum_vandermonde_determinant_identity_error=determinant_error,
                epsilon_budget=budget,
                checks=dict(all_five_signed_target_rows=worst < mp.mpf('1e-60'),
                            exact_scaled_bump_determinants=determinant_error < mp.mpf('1e-65'),
                            fixed_lambda_degeneracy_is_rejected=degeneracy_rejected,
                            fixed_small_lambda_has_costly_inverse=rows[-2]['maximum_inverse_entry'] > 10*rows[0]['maximum_inverse_entry'],
                            mean_update_has_conditional_positive_gain=budget['moment_defect_gain'] > 0 and budget['minimum_margin'] >= 0,
                            full_cycle_not_inferred=budget['full_residual_cycle_gain_verified'] is False))
