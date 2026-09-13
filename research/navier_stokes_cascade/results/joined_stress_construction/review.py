"""Independent differentiation controls for the proposed stress construction.

The actual matching target is not replaced by these manufactured functions.
Quadrature tests an identity of the moment map; it does not certify its exact
integrals or a continuum derivative supremum.
"""
import sys
import mpmath as mp
from construction import RESULTS, module_at

sys.path.append(str(RESULTS/'axis_matching_audit'))
from moments import MomentMap

loop = module_at('joined_stress_loop', RESULTS/'stress_realization_audit'/'loop.py')


def curvature_identity(digits=80):
    with mp.workdps(digits):
        m = MomentMap(16, 4)
        eta, pressure = mp.mpf('.37'), mp.mpf(16)

        def c_at(e):
            return mp.matrix([mp.mpf('.01')*(j+1)*(1+e+(j+1)*e*e/10)
                              for j in range(5)])

        def target(e):
            c = c_at(e)
            return m.B*c+m.quadratic(c, c, e, pressure)

        def q_partial(c, d, order):
            # Differentiate the polynomial zeta itself, not a frozen Jet value.
            factor = mp.exp(mp.mpf('1.2'))/pressure**2
            factor *= 4*eta*(1+eta**2) if order == 1 else 4+12*eta**2
            q = mp.matrix(5, 1)
            q[3] = -factor*mp.fsum(m.QU[j]*c[j]*d[j] for j in range(2))
            return q

        c = c_at(eta)
        c1 = mp.matrix([mp.diff(lambda e: c_at(e)[j], eta) for j in range(5)])
        c2 = mp.matrix([mp.diff(lambda e: c_at(e)[j], eta, 2) for j in range(5)])
        z2 = mp.matrix([mp.diff(lambda e: target(e)[j], eta, 2) for j in range(5)])
        jac = m.B.copy()
        for j in range(5):
            unit = mp.matrix(5, 1)
            unit[j] = 1
            column = 2*m.quadratic(c, unit, eta, pressure)
            for i in range(5):
                jac[i, j] += column[i]
        common = z2-2*m.quadratic(c1, c1, eta, pressure)
        cross, curvature = 4*q_partial(c, c1, 1), q_partial(c, c, 2)
        correct = mp.lu_solve(jac, common-cross-curvature)
        missing_cross = mp.lu_solve(jac, common-curvature)
        frozen_zeta = mp.lu_solve(jac, common)
        norm = lambda v: max(abs(x) for x in v)
        return dict(manufactured_target=True, quadrature_order=16, quadrature_panels=4,
                    correct_absolute_error=norm(correct-c2),
                    omitted_cross_absolute_error=norm(missing_cross-c2),
                    frozen_zeta_absolute_error=norm(frozen_zeta-c2),
                    checks=dict(second_derivative_identity=norm(correct-c2) < mp.mpf('1e-65'),
                                missing_cross_detected=norm(missing_cross-c2) > mp.mpf('1e-9'),
                                frozen_zeta_detected=norm(frozen_zeta-c2) > mp.mpf('1e-9')),
                    actual_target_C2_verified=False)


def loop_derivative_controls(digits=80):
    with mp.workdps(digits):
        rows = []
        d0 = mp.mpf('.125')
        for mu, p in [('1', '0'), ('.25', '1e-20'), ('2', '.75'), ('2', '-.75')]:
            mu, p = mp.mpf(mu), mp.mpf(p)
            slope = mp.diff(lambda u: mp.sqrt(loop.variance(u, p, d0)), mu)
            bound = d0*mp.exp(-4*mu*abs(p))/32
            rows.append(dict(mu=mu, p=p, slope=slope, lower_bound=bound,
                             passed=bool(slope >= bound)))
        return dict(manufactured_input=True, rows=rows,
                    zero_amplitude_slope=d0/mp.sqrt(2),
                    check=all(row['passed'] for row in rows),
                    actual_loop_derivative_envelope_verified=False)


def rounded_endpoint_control():
    with mp.workdps(80):
        # A manageable example of the same loss; this is NOT the actual C.
        tc = mp.exp(-2000)/100
        base = mp.mpf(1)
        return dict(manufactured_scale=True, digits=80,
                    rounded_offsets_collapse=bool(base+tc/4 == base+tc/2 == base),
                    exact_offset_order_preserved='0 < tc/4 < tc/2 < tc, because tc>0')
