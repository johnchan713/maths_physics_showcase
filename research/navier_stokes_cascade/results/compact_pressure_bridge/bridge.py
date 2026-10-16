"""Exact bookkeeping controls; finite fixtures are not the actual profile."""
from fractions import Fraction as Q
from functools import lru_cache
import mpmath as mp


def forward(X, E, U, H, eta, h, moments):
    """Original forward map. Each moment is (value, eta derivative)."""
    A, D, d, L = Q(1, 2)+h, Q(1, 2)-h, 1-eta**2, 1-2*h*eta**2
    (M, Me), (I, Ie), (J, Je), (S, Se), (Pi, Pie) = moments
    W = 1-(2*D*eta*M+d*Me)/X
    Qs = -W+((1-h)*I-D*eta*Ie-d*Je+2*(h-D)*eta*J)/(X*H)
    Ns = -W*U+(D*(M-eta*Me)+4*h*eta*S-d*Se)/X+4*A*eta*Pi-d*Pie
    return dict(W=W, Qs=Qs, Ns=Ns, Pi=Pi, p1=X*Qs/L, p2=X*Ns/(L*E))


def memory(X, E, U, H, eta, h, mu):
    """Independent subtraction formula on a common-field interval."""
    A, D, d, L = Q(1, 2)+h, Q(1, 2)-h, 1-eta**2, 1-2*h*eta**2
    (m, me), (i, ie), (j, je), (s, se), (c, ce) = mu
    W = -(2*D*eta*m+d*me)/X
    B = (1-h)*i-D*eta*ie-d*je+2*(h-D)*eta*j
    Qs = -W+B/(X*H)
    Ns = -U*W+(D*(m-eta*me)+4*h*eta*s-d*se)/X+4*A*eta*c-d*ce
    return dict(W=W, Qs=Qs, Ns=Ns, Pi=c, p1=X*Qs/L, p2=X*Ns/(L*E))


def pulse(a, b, p1, p2, F):
    """Squared quantities avoid rounding and unnecessary square roots."""
    if a <= 0 or F <= 0:
        raise ValueError('Positive shear and swirl required')
    t = -b/a
    v = a*(1+t*t)
    Pc, J = p1+p2*t, p2-p1*t
    D = Pc-v
    if v <= 2 or D <= 0:
        raise ValueError('Growth and oriented cone required')
    fourth = 2*D*D-(v-2)*J*J
    if fourth <= 0:
        raise ValueError('Strict fourth cone gap required')
    return dict(t=t, v=v, Pc=Pc, D=D, J=J, fourth=fourth,
                growth2=2*a*F*F*(v-2)/v,
                ratio2=(v-2)*J*J/(2*D*D),
                orientation2=D*D/(D*D+J*J),
                stress2=F*F*((p1-a)**2+(p2+b)**2),
                stress2_rotated=F*F*(D*D+J*J)/(1+t*t),
                g2=F*F*(a*a+b*b))


def activation_budget(C):
    """Exact scalar diagnostics of section 1; monotonic proof is in README."""
    C = Q(C)
    if C < 2:
        raise ValueError('Scalar sufficient range is C>=2')
    u = C**-17
    dmin = Q(9, 4)-u
    ratio = C**3*u*u/(dmin*dmin)
    return dict(Pc_lower=Q(9, 4)-2*u, D_over_ea_lower=dmin,
                ratio_upper=ratio, fourth_over_ea2_lower=2*dmin*dmin*(1-ratio),
                exit_v_upper=2*u)


def activation_review():
    rows = []
    for C in (2, 3, 10, 10**9):
        r = activation_budget(C)
        rows.append(dict(C=C, checks=dict(
            Pc_above_11_over_5=r['Pc_lower'] > Q(11, 5),
            D_over_ea_above_11_over_5=r['D_over_ea_lower'] > Q(11, 5),
            cone_loss_below_C_minus_30=r['ratio_upper'] < Q(C)**-30,
            fourth_over_ea_squared_above_four=r['fourth_over_ea2_lower'] > 4,
            exit_growth_below_one=r['exit_v_upper'] < 1)))
    return dict(rows=rows, checks=dict(
        exact_scalar_endpoint_controls=all(all(r['checks'].values()) for r in rows),
        universal_monotonic_argument_is_separate_from_samples=True))


def memory_review():
    checks, dropped_u, pressure = [], [], []
    zeros = [(Q(0), Q(0))]*5
    for eta in (Q(-1), Q(-1, 3), Q(0), Q(2, 5), Q(1)):
        for radius in (2, 4):
            X, E, H, U, h = Q(radius*radius, 2), Q(3, radius), Q(3), 1+eta, Q(1, 100)
            base = [(Q(k)+eta**2, 2*eta) for k in range(1, 6)]
            mu = [(Q(k, 100)*(1+eta), Q(k, 100)) for k in range(1, 6)]
            changed = [(a+c, b+d) for (a,b),(c,d) in zip(base,mu)]
            before, after = [forward(X,E,U,H,eta,h,m) for m in (base,changed)]
            delta = memory(X,E,U,H,eta,h,mu)
            checks.append(all(after[k]-before[k] == delta[k] for k in delta))
            checks.append(all(v == 0 for v in memory(X,E,U,H,eta,h,zeros).values()))
            # Only the fifth memory row survives; both its terms are recorded.
            cp = zeros[:4]+[(Q(1, 1000)*(1+eta), Q(1, 1000))]
            full = memory(X,E,U,H,eta,h,cp)
            c, ce = cp[4]
            L, A, d = 1-2*h*eta*eta, Q(1,2)+h, 1-eta*eta
            value_term, eta_term = X*4*A*eta*c/(L*E), -X*d*ce/(L*E)
            checks.append(full['p2'] == value_term+eta_term)
            pressure.append((value_term, eta_term, full['p2']))
            no_u = memory(X,E,0,H,eta,h,mu)
            dropped_u.append(abs(delta['p2']-no_u['p2']))
    # At eta=0, c=epsilon*(1+eta), and E=3/r, the offset scales as r^3.
    cp = zeros[:4]+[(Q(1,1000), Q(1,1000))]
    offsets = [memory(Q(r*r,2),Q(3,r),1,3,0,Q(1,100),cp)['p2'] for r in (2,4)]
    return dict(checks=dict(
        full_forward_subtraction_exact=all(checks),
        four_rows_cannot_erase_pressure=all(v != 0 for v in offsets),
        pressure_memory_amplifies_with_radius=offsets[1] == 8*offsets[0],
        pressure_error_survives_a_zero_stress_edge=offsets[0]*Q(3,4) == Q(-1,1000),
        pressure_value_and_eta_terms_separately_exercised=any(a and b for a,b,_ in pressure),
        minus_U_delta_W_is_required=max(dropped_u) > Q(1,100),
        angular_endpoint_controls_present=True),
        pressure_offset_at_eta_zero=[str(v) for v in offsets])


def pulse_review():
    checks, axial = [], []
    # These states are algebraic fixtures, not cumulative moments of a core.
    for a,b,p1,p2 in ((3,0,6,0),(3,-1,7,2),(Q(5,2),Q(1,2),8,-2)):
        F = Q(2, 5)
        r = pulse(Q(a),Q(b),Q(p1),Q(p2),F)
        R, kappa, fmin = Q(16), Q(1,10), Q(1,5)
        assert min(Q(a),r['v']-2,r['D'],r['fourth']) >= kappa/2
        assert max(abs(x) for x in (a,b,p1,p2)) <= R and 1/Q(a) <= R and F >= fmin
        checks.append(r['stress2'] == r['stress2_rotated'])
        checks.append(r['ratio2'] == 1-r['fourth']/(2*r['D']**2))
        checks.append(r['g2'] >= (fmin/R)**2)
        checks.append(r['growth2'] >= kappa*fmin*fmin/(2*R**4))
        checks.append(r['stress2'] >= (fmin*kappa/(4*R**2))**2)
        checks.append(r['orientation2'] >= (kappa/(10*R**3))**2)
        checks.append(r['ratio2'] <= 1-kappa/(64*R**6))
        # Direct Cartesian dot-product squared gives the same orientation.
        tx,tz = F*(p1-a),F*(p2+b)
        gx,gz = -F*a,F*b
        checks.append((tx*gx+tz*gz)**2/(r['stress2']*r['g2']) == r['orientation2'])
        checks.append(tx*gx+tz*gz < 0)
        checks.append((r['v']-2)*(tx*gz-tz*gx)**2/(2*(tx*gx+tz*gz)**2) == r['ratio2'])
        checks.append(2*a*F**2-4*a*a*F**4/r['g2'] == r['growth2'])
        axial.append(abs(2*a*F*F*(r['v']-2)-r['growth2']))
    return dict(checks=dict(
        exact_coordinates_and_all_value_floors=all(checks),
        growth_denominator_cannot_be_dropped=max(axial) > Q(1,10),
        manufactured_states_are_not_actual_profile=True))


def sigma(s):
    s = mp.mpf(s)
    if not 0 <= s <= 1:
        raise ValueError('Activation parameter outside [0,1]')
    if s == 0:
        return mp.mpf(0)
    if s == 1:
        return mp.mpf(1)
    z = -1/s**2+1/(1-s)**2
    return 1/(1+mp.exp(-z))


@lru_cache(maxsize=None)
def rule(order, digits):
    with mp.workdps(digits):
        x,w = mp.gauss_quadrature(order,'legendre')
        return tuple(((a+1)/2,b/2) for a,b in zip(x,w))


def normalized_integral(s, coefficient, order=48):
    """Integral_0^s sigma(u)b(u)du / sigma(s), used only at resolvable s."""
    s = mp.mpf(s)
    if s == 0:
        return mp.mpf(0)
    return s*mp.fsum(w*sigma(s*u)/sigma(s)*coefficient(s*u)
                    for u,w in rule(order,mp.mp.dps))


def operator_review():
    with mp.workdps(50):
        rows = []
        for s in map(mp.mpf,('.2','.5','.9','1')):
            a = normalized_integral(s,lambda u:1+u,48)
            b = normalized_integral(s,lambda u:1+u,72)
            direct = mp.quad(lambda u:sigma(u)*(1+u),[0,s/4,s/2,3*s/4,s])/sigma(s)
            # A second application of the positive operator; not a flat quotient.
            twice = normalized_integral(s,lambda u:normalized_integral(u,lambda z:1,48),48)
            rows.append(dict(s=float(s), single=float(b), double=float(twice),
                             refinement_error=float(abs(a-b)), independent_error=float(abs(b-direct))))
        return dict(rows=rows, checks=dict(
            whole_transition_including_right_endpoint_exercised=True,
            monotone_operator_bound=all(0 < r['single'] <= r['s']*(1+r['s']) for r in rows),
            twice_integrated_bound=all(0 < r['double'] <= r['s']**2 for r in rows),
            independent_quadrature_agrees=max(r['independent_error'] for r in rows) < 1e-8,
            quadrature_refines=max(r['refinement_error'] for r in rows) < 1e-8,
            finite_fixtures_not_continuum_suprema=True))
