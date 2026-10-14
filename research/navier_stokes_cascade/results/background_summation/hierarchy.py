"""General-order background system with regular lower-order sources.

F_n=E_n/sqrt(2X) is the regular swirl profile; M_n is the integrated
axial profile. No enormous-frequency leading profile is sampled here.
"""
from fractions import Fraction
from pathlib import Path
import math
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
sys.path.insert(0, str(HERE.parent/'physical_residual_budget'))
from residual import Jet, parameters, profile_jets, axial_operator, axial_second
sys.path.insert(0, str(HERE))

DIAGONAL = np.array([0., 0., 2., 0., 3., 1.])


def regular_v(profile, X, eta, h, lam):
    """V_n/X and its two-derivative jet, regular also at X=0."""
    _, D = parameters(h)
    x, e = Jet.variable(X, 0), Jet.variable(eta, 1)
    average = profile.average_U(x, e)
    average_eta = profile.average_U.derivative(1)(x, e)
    return (2*e*profile.U(x, e)-2*(D+lam)*e*average
            -(1-e*e)*average_eta)/(1-2*h*e*e)


def lower_sources(profiles, X, eta, h, nu=1., omit_radial_axial_diffusion=False):
    """Known sources at order n=len(profiles), with all lower orders fixed.

    The regular identity Z_k(X*v)=X*Z_(k-1)(v) removes the apparent
    axis singularity in the shifted radial axial-viscosity source.
    """
    A, D = parameters(h)
    n = len(profiles)
    if n < 1 or X < 0 or abs(eta) > 1 or not math.isfinite(nu) or nu <= 0:
        raise ValueError('Known leading profile, nonnegative radius and positive viscosity required')
    b, c = -A-.5, -A
    jets = [profile_jets(p, X, eta, h, 2*h*j) for j, p in enumerate(profiles)]
    v = [regular_v(p, X, eta, h, 2*h*j) for j, p in enumerate(profiles)]
    lam_last = 2*h*(n-1)
    L, d = 1-2*h*eta*eta, 1-eta*eta
    last = v[-1]
    omega_over_X = (((1-lam_last)*last.v+X*last.g[0]+D*eta*last.g[1])/L
                    -2*nu*(2*last.g[0]+X*last.H[0, 0]))
    for i in range(n):
        j = n-1-i
        vi, vj, Ui = v[i], v[j], jets[i]['U'].v
        omega_over_X += vi.v*(vj.v/2+X*vj.g[0])
        omega_over_X += Ui*(2*(2*h*j-1)*eta*vj.v-2*eta*X*vj.g[0]+d*vj.g[1])/L
    radial_axial_term = 0.
    if n >= 2:
        radial_axial_term = nu*axial_second(2*h*(n-2)-1, v[n-2], X, eta, h)
        if not omit_radial_axial_diffusion:
            omega_over_X -= radial_axial_term

    pressure = -.5*omega_over_X
    theta, axial = 0., 0.
    for i in range(1, n):
        j = n-i
        Fi, Fj, Ui, Uj = jets[i]['F'], jets[j]['F'], jets[i]['U'], jets[j]['U']
        pressure += Fi.v*Fj.v
        theta += v[i].v*(X*Fj.g[0]+Fj.v)+Ui.v*axial_operator(b+2*h*j, Fj, X, eta, h)
        axial += X*v[i].v*Uj.g[0]+Ui.v*axial_operator(c+2*h*j, Uj, X, eta, h)
    theta -= nu*axial_second(b+lam_last, jets[-1]['F'], X, eta, h)
    axial -= nu*axial_second(c+lam_last, jets[-1]['U'], X, eta, h)
    return dict(pressure=pressure, theta=theta, axial=axial,
                omega_over_X=omega_over_X, radial_axial_term_over_X=radial_axial_term)


def matrices(profiles, xi, eta, h, nu=1., omit_radial_axial_diffusion=False):
    """A0_n, A1 and f_n for W_n=(F,U,M/X-U,Pi,F_xi,U_xi)."""
    A, D = parameters(h)
    if not math.isfinite(xi) or xi < 0:
        raise ValueError('Nonnegative finite xi required')
    n, X = len(profiles), xi*xi
    source = lower_sources(profiles, X, eta, h, nu, omit_radial_axial_diffusion)
    lam, b, c = 2*h*n, -A-.5, -A
    kp, L, d = -2*A+lam, 1-2*h*eta*eta, 1-eta*eta
    base = profile_jets(profiles[0], X, eta, h, 0.)
    F, U = base['F'], base['U']
    v0 = regular_v(profiles[0], X, eta, h, 0.).v
    Hc, JF, JU = D*eta+d*U.v, F.v+X*F.g[0], X*U.g[0]
    ZF, ZU = axial_operator(b, F, X, eta, h), axial_operator(c, U, X, eta, h)
    a0, a1 = np.zeros((6, 6)), np.zeros((6, 6))
    a0[0, 4] = a0[1, 5] = 1.
    a0[2, 5], a0[3, 0] = -1., 4*xi*F.v
    transport = xi*(.5-eta*U.v)/L+xi*v0/2
    a0[4, 0] = 2*((b+lam)*(2*eta*U.v-1)/L+v0)/nu
    a0[4, 1] = 2*(2*(A-lam)*eta*JF/L+ZF)/nu
    a0[4, 2] = -4*(D+lam)*eta*JF/(nu*L)
    a0[5, 0] = -8*eta*X*F.v/(nu*L)
    a0[5, 1] = 2*((c+lam)*(2*eta*U.v-1)/L+ZU+2*(A-lam)*eta*JU/L)/nu
    a0[5, 2] = -4*(D+lam)*eta*JU/(nu*L)
    a0[5, 3] = 4*kp*eta/(nu*L)
    a0[4, 4] = a0[5, 5] = 2*transport/nu
    a1[4, 0] = 2*Hc/(nu*L)
    a1[4, 1] = a1[4, 2] = -2*d*JF/(nu*L)
    a1[5, 1] = 2*(Hc-d*JU)/(nu*L)
    a1[5, 2] = -2*d*JU/(nu*L)
    a1[5, 3] = 2*d/(nu*L)
    f = np.array([0., 0., 0., 2*xi*source['pressure'],
                  2*source['theta']/nu,
                  2*(source['axial']-2*eta*X*source['pressure']/L)/nu])
    return a0, a1, f


def finite_residual_exponent(h, order, derivative_order=0):
    """Worst physical exponent after retaining orders 0,...,N and stress."""
    h = Fraction(h)
    if not 0 < h < Fraction(1, 100) or any(type(v) is not int or v < 0 for v in (order, derivative_order)):
        raise ValueError('Construction h and nonnegative integer orders required')
    return 2*h*(order+1)-(Fraction(3, 2)+2*h+derivative_order)


def induction_ledger():
    return dict(collar='0<=X<=Xa*exp(tc/4), the same at every order',
                source_uses_finalized_lower_orders_only=True,
                angular_strip_may_shrink_with_order=True,
                no_uniform_coefficient_bound_in_order_asserted=True,
                analytic_solution='sparse Picard series with a finite C_n on a smaller eta neighborhood',
                zero_axis_traces=['F_n', 'U_n', 'Pi_n'],
                current_moment_map='B*c=-z_n, with all lower-order products held fixed',
                radial_axial_diffusion_starts_in_pressure_order=2,
                all_velocity_and_pressure_support='[0,Xv]',
                signed_stress_order_one_support='[Xkeep,Xb]',
                signed_stress_orders_two_and_above_support='[Xkeep,Xv]',
                higher_stress_over_zeta_has_all_finite_derivative_bounds=True,
                velocities_vanish_on_Imean=True, pressure_vanishes_on_Imean_not_asserted=True,
                finite_residual_loss='K_m=3/2+2h+m, independent of truncation order',
                actual_coefficient_sequence_numerically_materialized=False)
