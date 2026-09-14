"""Manufactured controls, NOT samples or suprema of the actual joined profile."""
from fractions import Fraction as F
import mpmath as mp
from jets import FirstJet, SecondJet, pressure_from_moments


def ex(v):
    return v.exp() if isinstance(v, SecondJet) else mp.exp(v)


def manufactured(y, eta):
    """Regular axis F=f0+f1 X, U=u0+u1 X and their exact five primitives."""
    X = ex(y)
    f0, f1 = 2+eta**2, (1+eta)/10
    u0, u1 = eta+eta**3/3, (1-eta**2)/5
    swirl_factor = f0+f1*X
    E = mp.sqrt(2)*ex(y/2)*swirl_factor
    U = u0+u1*X
    M = u0*X+u1*X**2/2
    I = f0*X**2+2*f1*X**3/3
    J = u0*f0*X**2+2*(u0*f1+u1*f0)*X**3/3+u1*f1*X**4/2
    S = u0**2*X+u0*u1*X**2+u1**2*X**3/3 \
        -f0**2*X**2/2-2*f0*f1*X**3/3-f1**2*X**4/4
    Cp = f0**2*X+f0*f1*X**2+f1**2*X**3/3
    return E, U, [M, I, J, S, Cp], -2-eta**2


def original_sources(y, eta, h):
    """Independent (4.9) path from the two field polynomials, not from (4.16)."""
    X = mp.exp(y)
    f0, f1 = 2+eta**2, (1+eta)/10
    u0, u1 = eta+eta**3/3, (1-eta**2)/5
    f0e, f1e, u0e, u1e = 2*eta, mp.mpf('.1'), 1+eta**2, -2*eta/5
    factor, U = f0+f1*X, u0+u1*X
    Ue, Uy = u0e+u1e*X, u1*X
    E2 = 2*X*factor**2
    ell = 1+X*f1/factor
    A, D, d = mp.mpf('.5')+h, mp.mpf('.5')-h, 1-eta**2
    W = 1-2*D*eta*(u0+u1*X/2)-d*(u0e+u1e*X/2)
    Hc = D*eta+d*U
    Pi = -2-eta**2+f0**2*X+f0*f1*X**2+f1**2*X**3/3
    Pie = -2*eta+2*f0*f0e*X+(f0e*f1+f0*f1e)*X**2+2*f1*f1e*X**3/3
    Sq = -W*ell-h*(1-2*eta*U)-Hc*(f0e+f1e*X)/factor
    Sn = -W*Uy-A*(1-2*eta*U)*U-Hc*Ue-d*Pie+4*A*eta*Pi+eta*E2
    return Sq, Sn, ell


def scalar_integrated_stresses(y, eta, h):
    """Scalar (4.16) path; differentiate its scalar moments independently."""
    E,U,m,axis = manufactured(y,eta)
    M,I,J,S,Cp = m
    Me,Ie,Je,Se,Cpe = [mp.diff(lambda e: manufactured(y,e)[2][i],eta) for i in range(5)]
    Pie = mp.diff(lambda e: manufactured(y,e)[3],eta)+Cpe
    X,H = mp.exp(y),mp.sqrt(2)*mp.exp(y/2)*E
    A,D,d,L = mp.mpf('.5')+h,mp.mpf('.5')-h,1-eta**2,1-2*h*eta**2
    W = 1-(2*D*eta*M+d*Me)/X
    Q = -W+((1-h)*I-D*eta*Ie-d*Je+2*(h-D)*eta*J)/(X*H)
    N = -W*U+(D*(M-eta*Me)+4*h*eta*S-d*Se)/X+4*A*eta*(axis+Cp)-d*Pie
    return X*Q/L,X*N/(L*E)


def source_and_jet_review(digits=80):
    with mp.workdps(digits):
        h = mp.mpf('.003')
        jet_error, source_error, radial_error, pressure_jet_error = [mp.mpf(0) for _ in range(4)]
        for y in map(mp.mpf, ['-.4', '.7']):
            for eta in map(mp.mpf, ['-1', '-.37', '0', '.6', '1']):
                data = manufactured(SecondJet(y, y=1), SecondJet(eta, e=1))
                E, U, moments, axis = data
                for i, jet in enumerate([E, U, *moments, axis]):
                    def value(z, e):
                        a, b, m, p = manufactured(z, e)
                        return [a, b, *m, p][i]
                    expected = [value(y, eta), *[mp.diff(value, (y, eta), orders)
                                    for orders in ((1,0), (0,1), (2,0), (1,1), (0,2))]]
                    jet_error = max(jet_error, *(abs(a-b) for a,b in zip(jet.values(), expected)))
                result = pressure_from_moments(y, eta, h, E, U, moments, axis)
                for j,key in enumerate(('p1','p2')):
                    fun = lambda q,e: scalar_integrated_stresses(q,e,h)[j]
                    expected = [fun(y,eta),mp.diff(fun,(y,eta),(1,0)),mp.diff(fun,(y,eta),(0,1))]
                    pressure_jet_error = max(pressure_jet_error,*(abs(a-b) for a,b in
                        zip((result[key].v,result[key].y,result[key].e),expected)))
                Sq, Sn, ell = original_sources(y, eta, h)
                X, L = mp.exp(y), 1-2*h*eta**2
                source_error = max(source_error,
                    abs(result['p1'].y-X*Sq/L+ell*result['p1'].v),
                    abs(result['Ns'].y+result['Ns'].v-Sn))
                # Independent radial formulas (including the pressure-axis density).
                dy = [X*U.v, mp.sqrt(2)*X**mp.mpf('1.5')*E.v,
                      mp.sqrt(2)*X**mp.mpf('1.5')*U.v*E.v,
                      X*(U.v**2-E.v**2/2), E.v**2/2]
                dyy = [X*(U.v+U.y),
                       mp.sqrt(2)*X**mp.mpf('1.5')*(3*E.v/2+E.y),
                       mp.sqrt(2)*X**mp.mpf('1.5')*(3*U.v*E.v/2+U.y*E.v+U.v*E.y),
                       X*(U.v**2-E.v**2/2+2*U.v*U.y-E.v*E.y), E.v*E.y]
                for m, d1, d2 in zip(moments, dy, dyy):
                    radial_error = max(radial_error, abs(m.y-d1), abs(m.yy-d2))
        y, eta = mp.mpf('.7'), mp.mpf('.6')
        E, U, moments, axis = manufactured(SecondJet(y, y=1), SecondJet(eta, e=1))
        full = pressure_from_moments(y, eta, h, E, U, moments, axis)
        frozen = [SecondJet(m.v,m.y,m.e,m.yy,0,0) for m in moments]
        bad = pressure_from_moments(y, eta, h, E, U, frozen, axis)
        preserved_values = max(abs(full[k].v-bad[k].v) for k in ('p1','p2'))
        curvature_error = max(abs(full[k].e-bad[k].e) for k in ('p1','p2'))
        zero_axis = pressure_from_moments(y, eta, h, E, U, moments, SecondJet(mp.mpf(0)))
        axis_error = abs(full['p2'].v-zero_axis['p2'].v)
        # Drop exactly the pressure terms from Ns, retaining all other terms.
        ej = FirstJet(eta, e=1)
        Pie = (axis+moments[4]).partial_first('eta')
        force = 4*(mp.mpf('.5')+h)*ej*full['Pi']-(1-ej*ej)*Pie
        omitted = full['Ns']-force
        missing_pressure_residual = abs(omitted.y+omitted.v-original_sources(y,eta,h)[1])
        checks = dict(mixed_jets_agree_with_direct_differentiation=jet_error < mp.mpf('1e-65'),
                      all_five_original_radial_identities=radial_error < mp.mpf('1e-65'),
                      original_source_ODEs_agree=source_error < mp.mpf('1e-65'),
                      pressure_jets_agree_with_scalar_differentiation=pressure_jet_error < mp.mpf('1e-65'),
                      missing_curvature_preserves_values=preserved_values == 0,
                      missing_angular_curvature_detected=curvature_error > mp.mpf('.01'),
                      changed_axis_datum_detected=axis_error > mp.mpf('.01'),
                      dropped_pressure_force_detected=missing_pressure_residual > mp.mpf('.01'))
        return dict(manufactured_regular_axis=True, maximum_jet_error=jet_error,
                    maximum_radial_identity_error=radial_error, maximum_source_ODE_error=source_error,
                    maximum_pressure_jet_error=pressure_jet_error,
                    frozen_curvature_value_error=preserved_values,
                    frozen_curvature_derivative_error=curvature_error,
                    changed_axis_datum_error=axis_error,
                    omitted_pressure_source_residual=missing_pressure_residual, checks=checks)


def cutoff_review(digits=80):
    with mp.workdps(digits):
        def step(z):
            return 1/(1+mp.exp(1/z**2-1/(1-z)**2))
        z, eta = mp.mpf('.37'), mp.mpf('.3')
        records = []
        for width in map(mp.mpf, ['.01', '.0001']):
            y = width*z
            zj = SecondJet(z, y=1/width)
            sigma = zj.compose(step(z), mp.diff(step,z), mp.diff(step,z,2))
            ej = SecondJet(eta, e=1)
            logfield = ej*ej+(1+ej)*sigma
            field = logfield.exp()
            scalar = lambda q,e: mp.exp(e*e+(1+e)*step(q/width))
            ye = mp.diff(scalar,(y,eta),(1,1))
            yy = mp.diff(scalar,(y,eta),(2,0))
            missing_product = field.v*logfield.ye
            records.append(dict(width=width, F=logfield, mixed_error=abs(field.ye-ye),
                                radial_error=abs(field.yy-yy),
                                omitted_exponential_product_error=abs(missing_product-ye)))
        a,b = records
        ratio_error = abs(b['F'].yy/a['F'].yy-10000)
        return dict(manufactured_cutoff=True, records=records,
                    checks=dict(equal_values_and_angular_derivatives=a['F'].v == b['F'].v
                                    and a['F'].e == b['F'].e and a['F'].ee == b['F'].ee,
                                second_radial_derivative_has_width_squared_loss=ratio_error < mp.mpf('1e-65'),
                                mixed_and_radial_exponential_rules=all(max(v['mixed_error'],v['radial_error'])
                                    < mp.mpf('1e-60') for v in records),
                                omitted_exponential_product_detected=all(v['omitted_exponential_product_error'] > 1
                                    for v in records)))


def shear_sign_review():
    a, bs = FirstJet(F(3),F(1),F(2)), FirstJet(F(2),F(3),F(-1))
    p1,p2 = a,-bs  # Exactly zero leading stress: ps=(a,-bs).
    ts = -bs/a
    vs = a*(1+ts*ts)
    Pc,Jc = p1+p2*ts,p2-p1*ts
    wrong = bs/a
    return dict(manufactured_zero_stress_direction=True,
                checks=dict(source_sign_gives_Pc_equal_vs=Pc == vs,
                            source_sign_gives_zero_Jc=Jc == FirstJet(F(0)),
                            reversed_shear_sign_detected=(p2-p1*wrong).v != 0))
