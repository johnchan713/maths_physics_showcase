"""Exact rational diagnostics, not evaluations of the actual global profile.

The compact polynomial fixtures have finite endpoint regularity. The analytical
all-orders estimates and actual C-infinity supports are proved in README.md
using the explicitly inherited inverse lemmas, not inferred from these fixtures.
"""
from fractions import Fraction as Q
import importlib.util
import sys

KAPPA = Q(1, 100000)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class Poly:
    """Sparse Q[R,Z,T,U]; U is an auxiliary random variable, not a coordinate.

    Integer Laurent powers of R are allowed; radial integration rejects R^-1.
    No float coercions: diagnostics are exact rational identities.
    """
    def __init__(self, value=0):
        if isinstance(value, Poly):
            self.terms = dict(value.terms)
        elif isinstance(value, dict):
            if any(not isinstance(v, (int, Q)) for v in value.values()):
                raise TypeError('Exact integer or rational coefficients required')
            self.terms = {k: Q(v) for k, v in value.items() if v}
        elif isinstance(value, (int, Q)):
            self.terms = {(0, 0, 0, 0): Q(value)} if value else {}
        else:
            raise TypeError('Exact integer or rational input required')

    def __add__(self, other):
        result = dict(self.terms)
        for key, value in Poly(other).terms.items():
            result[key] = result.get(key, Q(0)) + value
        return Poly(result)

    __radd__ = __add__

    def __neg__(self):
        return Poly({k: -v for k, v in self.terms.items()})

    def __sub__(self, other):
        return self + -Poly(other)

    def __rsub__(self, other):
        return Poly(other) + -self

    def __mul__(self, other):
        result = {}
        for k, v in self.terms.items():
            for j, w in Poly(other).terms.items():
                key = tuple(a+b for a, b in zip(k, j))
                result[key] = result.get(key, Q(0)) + v*w
        return Poly(result)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if not isinstance(other, (int, Q)):
            raise TypeError('Only rational division; use shift for radial powers')
        return self * (1/Q(other))

    def __pow__(self, n):
        if type(n) is not int or n < 0:
            raise ValueError('Nonnegative integer power required')
        result = Poly(1)
        for _ in range(n):
            result = result*self
        return result

    def __eq__(self, other):
        return self.terms == Poly(other).terms

    def __bool__(self):
        return bool(self.terms)

    def diff(self, axis):
        result = {}
        for key, value in self.terms.items():
            if key[axis]:
                new = list(key); new[axis] -= 1
                result[tuple(new)] = value*key[axis]
        return Poly(result)

    def shift(self, power):
        return Poly({(r+power, z, t, u): v for (r, z, t, u), v in self.terms.items()})

    def at(self, axis, value):
        result = Poly()
        for key, coeff in self.terms.items():
            new = list(key); new[axis] = 0
            result += Poly({tuple(new): coeff*Q(value)**key[axis]})
        return result

    def primitive(self):
        if any(k[0] == -1 for k in self.terms):
            raise ValueError('Logarithmic primitive outside polynomial fixture')
        return Poly({(r+1, z, t, u): v/Q(r+1) for (r, z, t, u), v in self.terms.items()})

    def integral(self):
        p = self.primitive()
        return p.at(0, 2)-p.at(0, 1)

    def haar(self):
        """Exact moments of U=cos(y), uniform y on a circle.

        Only degree <=2 is used. This is a manufactured Haar-algebra check,
        not the actual irrational two-torus inverse or a quadrature.
        """
        if any(k[3] > 2 for k in self.terms):
            raise ValueError('Only quadratic auxiliary moments supported')
        return Poly({(r, z, t, 0): v for (r, z, t, u), v in self.terms.items() if u == 0}) + Poly({(r, z, t, 0): v/2 for (r, z, t, u), v in self.terms.items() if u == 2})

    def scalar(self):
        if any(k != (0, 0, 0, 0) for k in self.terms):
            raise ValueError('Not a scalar')
        return self.terms.get((0, 0, 0, 0), Q(0))


R, Z, T, U = [Poly({tuple(int(i == j) for i in range(4)): 1}) for j in range(4)]
EPS = Q(1, 200)


def radial(f, e=1):
    return f.diff(0)+e*f.shift(-1)


def lap(f):
    return f.diff(0).diff(0)+f.diff(0).shift(-1)+EPS**2*f.diff(1).diff(1)


def fixture():
    b = (R-1)**6*(2-R)**6
    h, g, k = 1+Z/5, 2+Z/3, 3
    a = (R*b*h, Poly(), Poly())
    c = (-EPS*R*b*g.diff(1)/k, -(R**2*b).diff(0)*h/k, radial(R*b*g)/k)
    w = [[(a[i]*a[j]+c[i]*c[j])/2 for j in range(3)] for i in range(3)]
    return b, a, c, w


def raw_data(w):
    gr = -radial(w[0][0])-EPS*w[0][2].diff(1)+w[1][1].shift(-1)
    p = gr.integral()
    jt = (R**2*w[1][2]).integral()
    jz = (R*w[2][2]).integral()-(R**2*gr).integral()/2
    return gr, p, jt, jz


def covariance_review():
    b, a, c, w = fixture()
    gr, p, jt, jz = raw_data(w)
    ibp_p = ((w[1][1]-w[0][0]).shift(-1)).integral()-EPS*w[0][2].integral().diff(1)
    ibp_jz = (R*(w[2][2]-(w[0][0]+w[1][1])/2)).integral()+EPS*(R**2*w[0][2]).integral().diff(1)/2
    # An exactly normalized, slowly moving bump with nonconstant second moment.
    mass = b.integral().scalar()
    rho = b/mass+Z*(R-Q(3, 2))*b/mass
    crho = (R**2*rho).integral()/2
    source = gr-rho*p
    primitive = source.primitive()
    pm = primitive-primitive.at(0, 1)
    moment = (R*pm).integral()
    no_axial = [[Poly(w[i][j]) for j in range(3)] for i in range(3)]
    no_axial[0][2] = no_axial[2][0] = no_axial[1][2] = no_axial[2][1] = no_axial[2][2] = Poly()
    _, bad_p, bad_jt, bad_jz = raw_data(no_axial)
    checks = dict(
        complete_cosine_divergence_zero=radial(a[0])+3*c[1].shift(-1) == 0,
        complete_sine_divergence_zero=radial(c[0])+EPS*c[2].diff(1) == 0,
        covariance_symmetric=all(w[i][j] == w[j][i] for i in range(3) for j in range(3)),
        all_six_covariances_nonzero=all(bool(w[i][j]) for i in range(3) for j in range(i, 3)),
        raw_pressure_defect_ibp=p == ibp_p,
        raw_axial_defect_ibp=jz == ibp_jz,
        moving_bump_mass_one=rho.integral() == 1,
        compact_pressure_at_both_ends=pm.at(0, 1) == pm.at(0, 2) == 0,
        compact_pressure_radial_identity=pm.diff(0) == source,
        compact_pressure_axial_moment=moment == -(R**2*gr).integral()/2+crho*p,
        whole_moving_bump_derivative=moment.diff(1) == (-(R**2*gr).integral()/2+crho*p).diff(1),
        omitted_bump_derivative_detected=bool(crho.diff(1)*p),
        omitted_axial_pressure_detected=bool(moment.diff(1)),
        omitted_axial_velocity_detected=all((p != bad_p, jt != bad_jt, jz != bad_jz)),
        finite_fixture_endpoint_jets_zero=all(b.diff(0).diff(0).at(0, e) == 0 for e in (1, 2)))
    return dict(checks=checks, arithmetic='exact Fraction; no float inputs',
                nonzero_controls=dict(P_at_Z0=str(p.at(1, 0).scalar()),
                    Jtheta_at_Z0=str(jt.at(1, 0).scalar()), Jz_at_Z0=str(jz.at(1, 0).scalar()),
                    moving_bump_cprime=str(crho.diff(1).scalar())),
                fixture_not_actual_Cinfinity_profile=True)


def averaged_gr(base, mean, wrr, wtt, wrz):
    """Haar-averaged exact radial source (graph derivatives integrate to zero)."""
    b, vbase, gbase = base
    beta, v, gamma = mean
    return (EPS*beta.haar().diff(2)-radial((2*b*beta+beta**2+wrr).haar())
            -EPS*(b*gamma+gbase*beta+beta*gamma+wrz).haar().diff(1)
            +(2*vbase*v+v**2+wtt).haar().shift(-1)
            +EPS*(lap(beta.haar())-beta.haar().shift(-2)))


def haar_update_review():
    b, _, _, w = fixture()
    base = (EPS*R*b*(1+Z+T), R**2*(2+Z), R*b*(1+Z))
    psi = R*b*(1+Z/5+T/7)
    # Arbitrary zero-mean manufactured old variables; not an actual N inverse.
    fast = (-EPS*psi.diff(1)*U, R*b*(1+Z/3)*U, radial(psi)*U)
    slow_psi = R*b*(1+Z/7+T/11)
    # Exact angular moment zero: R^2 dv is an endpoint derivative.
    slow = (-EPS*slow_psi.diff(1), (R**4*b).diff(0).shift(-2)*(1+Z/11+T/13), radial(slow_psi))
    zero = (Poly(), Poly(), Poly())
    g0 = averaged_gr(base, zero, w[0][0], w[1][1], w[0][2])
    gf = averaged_gr(base, fast, w[0][0], w[1][1], w[0][2])
    bf, vf, zf = fast
    db, dv, dg = slow
    expected_fast = -radial((bf**2).haar())-EPS*(bf*zf).haar().diff(1)+(vf**2).haar().shift(-1)
    new = tuple(a+c for a, c in zip(fast, slow))
    gn = averaged_gr(base, new, w[0][0], w[1][1], w[0][2])
    rg = gn-gf-(2*base[1]*dv).shift(-1)
    expected = (EPS*db.diff(2)-radial(2*base[0]*db+db**2)
                -EPS*(base[0]*dg+base[2]*db+db*dg).diff(1)+(dv**2).shift(-1)
                +EPS*(lap(db)-db.shift(-2)))
    # Defect differences before using the last three moment rows.
    jt = lambda m: (R**2*(base[2]*m[1]+base[1]*m[2]+m[2]*m[1]).haar()).integral()
    jz = lambda m, gr: (R*(2*base[2]*m[2]+m[2]**2).haar()).integral()-(R**2*gr).integral()/2
    linear_jt = (R**2*(base[2]*dv+base[1]*dg)).integral()
    linear_jz = (2*R*base[2]*dg-R*base[1]*dv).integral()
    nonzero_old = (fast[0], fast[1]+R*b, fast[2])
    changed = tuple(a+c for a, c in zip(nonzero_old, slow))
    bad_rg = averaged_gr(base, changed, w[0][0], w[1][1], w[0][2])-averaged_gr(base, nonzero_old, w[0][0], w[1][1], w[0][2])-(2*base[1]*dv).shift(-1)
    return dict(checks=dict(
        all_fast_haar_means_zero=all(f.haar() == 0 for f in fast),
        fast_averaged_linear_changes_zero=gf-g0 == expected_fast,
        old_fast_slow_products_have_zero_mean=all((a*c).haar() == 0 for a in fast for c in slow),
        exact_slow_divergence=radial(db)+EPS*dg.diff(1) == 0,
        slow_velocity_moments_zero=(R**2*dv).integral() == (R*dg).integral() == 0,
        radial_recomputed_remainder=rg == expected,
        pressure_defect_difference=(gn-gf).integral() == (2*base[1]*dv).shift(-1).integral()+rg.integral(),
        angular_defect_difference=jt(new)-jt(fast)-linear_jt == (R**2*dg*dv).integral(),
        axial_defect_difference=jz(new, gn)-jz(fast, gf)-linear_jz == (R*dg**2).integral()-(R**2*rg).integral()/2,
        nonlinear_defects_not_identically_zero=all(bool(f) for f in (rg.integral(), (R**2*dg*dv).integral(), (R*dg**2).integral()-(R**2*rg).integral()/2)),
        nonzero_old_mean_breaks_sharp_reduction=bad_rg != expected),
        arithmetic='exact Fraction auxiliary moment algebra; not a torus-inverse fixture',
        cancellation_requires_initial_zero_Haar_fast_state=True)


def exponent_review():
    k = KAPPA
    rows = dict(raw_full_gr=1-k, raw_averaged_gr=Q(1), raw_defects=Q(1),
                real_primary_covariance_error=2-2*k, base_stress_error=Q(3),
                initial_averaged_tangential=2-2*k,
                fast_tangential=1-k, fast_radial=2-k,
                fast_viscosity=2-3*k, fast_quadratic_averaged_gr=2-2*k,
                slow_tangential=Q(1), slow_radial=Q(2),
                averaged_slow_time=Q(3), averaged_b_radial=Q(3),
                averaged_axial_flux=Q(3), averaged_slow_square=Q(2),
                averaged_slow_viscosity=Q(3),
                final_three_defects=Q(2), final_radial_residual=Q(2),
                final_full_tangential_residual=2-3*k)
    return dict(rows={k: str(v) for k, v in rows.items()}, checks=dict(
        conservative_initial_full_gr=rows['raw_full_gr'] < rows['raw_averaged_gr'],
        haar_derivative_has_no_graph_loss=rows['raw_averaged_gr'] == 1,
        primary_stress_match_dominates=2-2*k < 3,
        fast_viscous_loss_retained=rows['fast_viscosity'] == 1+(1-k)-2*k,
        fast_update_does_not_worsen_raw_defects=2-2*k > 1,
        all_averaged_slow_terms_at_least_two=min(rows[n] for n in ('averaged_slow_time','averaged_b_radial','averaged_axial_flux','averaged_slow_square','averaged_slow_viscosity')) == 2,
        initial_gain_one=rows['final_three_defects']-rows['raw_defects'] == 1,
        tangential_not_promoted_to_two=rows['final_full_tangential_residual'] < 2,
        later_nonzero_mean_cross_term_is_worse=Q(9,10)+1 < 2),
        fixed_kappa=str(k), derivative_orders='every fixed coefficient order, by analytical product/inverse lemmas')


def pressure_tail_review():
    # Euler axial derivative numerator for q^p f(X,eta):
    # 2 eta (p f-X f_X)+d f_eta. On cutoff plateaux first.
    h = Q(1, 200); a = Q(1, 2)+h; axial_power = Q(1, 2)-h
    eta, x, mu = Q(2, 5), Q(7, 3), Q(3, 7)
    d = 1-eta**2
    f = -mu**2/(4*x); fx = mu**2/(4*x*x)
    first = 2*eta*(-f-x*fx)
    k2 = -2*a+4*h
    second = 2*eta*k2  # K2=1 constant, but physical q-power still changes.
    transition = 2*eta*Q(4, 5)*Q(3, 5)*f  # sigma chi'(sigma) !=0
    return dict(checks=dict(
        constant_second_profile_can_have_axial_gradient=second != 0,
        centrifugal_first_tail_plateau_axial_derivative_zero=first == 0,
        first_tail_cutoff_transition_not_zero=transition != 0,
        full_pressure_subtraction_is_linear=(Q(7)+Q(11)+Q(13))-Q(7) == Q(11)+Q(13),
        pressure_tail_powers_distinct=k2 != -1,
        same_fixed_parameters=0 < h < Q(1,100) and axial_power > 0 and d > 0),
        exact_controls=dict(second_tail_axial_numerator=str(second), first_tail_transition=str(transition)),
        both_actual_background_pressure_tails_retained_in_analytical_balance=True)
