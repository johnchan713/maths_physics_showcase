"""Manufactured controls of the actual bump and original moment equations."""
from functools import lru_cache
import mpmath as mp
from bounds import repair


def coeff(eta,amplitude):
    return mp.matrix([amplitude*(1+eta/10),amplitude*(-mp.mpf('.5')+3*eta/10),
                      amplitude*(mp.mpf('.3')-eta/5),
                      amplitude*(-mp.mpf('.2')+2*eta**2/5),
                      amplitude*(mp.mpf('.1')+eta/7)])


def bump(y,center,width):
    return repair.step_prime((y-center)/width+mp.mpf('.5'))/width


@lru_cache(maxsize=None)
def bump_review(digits=70):
    with mp.workdps(digits):
        eta,center,z = mp.mpf('.3'),mp.mpf('1.5'),mp.mpf('.31')
        lam,amplitude = mp.mpf('.003'),mp.mpf('1e-5')
        e = lambda t: mp.mpf('1.2')+t+t*t/20
        c = lambda t: amplitude*(mp.mpf('.3')-t/5)
        records=[]
        for width in map(mp.mpf,('.1','.05')):
            y=center+width*(z-mp.mpf('.5'))
            E0=lambda q,t:e(t)*mp.exp((-mp.mpf('.5')-lam)*q)
            de=lambda q,t:e(t)*c(t)*bump(q,center,width)
            E=E0(y,eta)+de(y,eta)
            Ey=mp.diff(lambda q:E0(q,eta)+de(q,eta),y)
            expected=1-2*Ey/E
            frozen_denominator=1-2*Ey/E0(y,eta)
            angular=mp.diff(lambda t:de(y,t),eta)
            frozen_e=e(eta)*mp.diff(c,eta)*bump(y,center,width)
            records.append(dict(width=width,radial_bump_derivative=mp.diff(lambda q:bump(q,center,width),y),
                                denominator_omission_error=abs(expected-frozen_denominator),
                                normalization_derivative_omission_error=abs(angular-frozen_e)))
        ratio=abs(records[1]['radial_bump_derivative']/records[0]['radial_bump_derivative']-4)
        peak=bump(center,center,mp.mpf('.1'))
        # A correction beyond the small neighborhood can make swirl negative.
        bad_c=-2*mp.exp((-mp.mpf('.5')-lam)*center)/peak
        bad_E=e(eta)*(mp.exp((-mp.mpf('.5')-lam)*center)+bad_c*peak)
        return dict(manufactured=True,records=records,
                    checks=dict(width_squared_cost=ratio < mp.mpf('1e-60'),
                                perturbed_denominator_required=all(r['denominator_omission_error'] > mp.mpf('1e-7') for r in records),
                                e_angular_derivative_required=all(r['normalization_derivative_omission_error'] > mp.mpf('1e-7') for r in records),
                                unrestricted_coefficients_break_positivity=bad_E < 0))


def evalue(eta):
    return mp.mpf('1.2')+eta/10+eta**2/20


def moments(y,eta,amplitude,model,drop_cp=False):
    """Exact power primitives plus a diagnostic five-bump moment change.

    The incoming moment values are manufactured and angular dependent. The
    normalization factors and all their derivatives are evaluated as functions.
    The evaluation point y>=3 is after every bump support.
    """
    rho,e,alpha = mp.mpf(3),evalue(eta),-mp.mpf('.5')-model.lam
    baseline=mp.matrix([1+eta,2+eta**2,eta/5,1+eta**2,mp.mpf('.4')+eta/10])
    baseline[1]+=mp.sqrt(2)*rho**mp.mpf('1.5')*e*mp.expm1((mp.mpf('1.5')+alpha)*y)/(mp.mpf('1.5')+alpha)
    baseline[3]-=rho*e**2*mp.expm1((1+2*alpha)*y)/(2*(1+2*alpha))
    baseline[4]+=e**2*mp.expm1(2*alpha*y)/(4*alpha)
    factors=[rho*e,mp.sqrt(2)*rho**mp.mpf('1.5')*e,
             mp.sqrt(2)*rho**mp.mpf('1.5')*e**2,rho*e**2,e**2]
    changes=model.ordinary_changes(coeff(eta,amplitude))
    if drop_cp:
        changes[4]=0
    return baseline+mp.matrix([factors[j]*changes[j] for j in range(5)])


def pressure_state(y,eta,amplitude,model,drop_cp=False,freeze_correction_eta=False):
    rho=mp.mpf(3)
    X=rho*mp.exp(y)
    E=evalue(eta)*mp.exp((-mp.mpf('.5')-model.lam)*y)
    h=model.lam**2
    A,D,d,L=mp.mpf('.5')+h,mp.mpf('.5')-h,1-eta**2,1-2*h*eta**2
    M,I,J,S,Cp=moments(y,eta,amplitude,model,drop_cp)
    amplitude_for_eta=0 if freeze_correction_eta else amplitude
    Me,Ie,Je,Se,Cpe=[mp.diff(lambda t:moments(y,t,amplitude_for_eta,model,drop_cp)[j],eta) for j in range(5)]
    Pi=-5-eta**2+Cp
    Pie=-2*eta+Cpe
    W=1-(2*D*eta*M+d*Me)/X
    Hx=mp.sqrt(2*X)*E
    Q=-W+((1-h)*I-D*eta*Ie-d*Je+2*(h-D)*eta*J)/(X*Hx)
    Ns=(D*(M-eta*Me)+4*h*eta*S-d*Se)/X+4*A*eta*Pi-d*Pie
    return dict(p1=X*Q/L,p2=X*Ns/(L*E),Ns=Ns,W=W,Pi=Pi,Pie=Pie,E=E,X=X)


@lru_cache(maxsize=None)
def moment_pressure_review(digits=70):
    with mp.workdps(digits):
        model=repair.Repair('.003',order=24,panels=8)
        eta,y,amplitude=mp.mpf('.3'),mp.mpf('3.1'),mp.mpf('1e-5')
        c=coeff(eta,amplitude)
        direct=model.direct_changes(c)
        mapped=model.ordinary_changes(c)
        moment_error=max(abs(x-z) for x,z in zip(direct,mapped))
        full=pressure_state(y,eta,amplitude,model)
        dropped=pressure_state(y,eta,amplitude,model,drop_cp=True)
        frozen=pressure_state(y,eta,amplitude,model,freeze_correction_eta=True)
        h=model.lam**2
        A,D,d,L=mp.mpf('.5')+h,mp.mpf('.5')-h,1-eta**2,1-2*h*eta**2
        ell=-model.lam
        Sq=-full['W']*ell-h-D*eta*mp.diff(evalue,eta)/evalue(eta)
        Sn=-d*full['Pie']+4*A*eta*full['Pi']+eta*full['E']**2
        radial_p=mp.diff(lambda q:pressure_state(q,eta,amplitude,model)['p1'],y)
        radial_n=mp.diff(lambda q:pressure_state(q,eta,amplitude,model)['Ns'],y)
        source_error=max(abs(radial_p+ell*full['p1']-full['X']*Sq/L),
                         abs(radial_n+full['Ns']-Sn))
        u_square=mp.fsum(model.squares[j]*c[j]**2 for j in range(2))
        return dict(manufactured=True,ordinary_moment_map_error=moment_error,
                    original_source_ODE_error=source_error,
                    dropped_pressure_moment_error=abs(full['p2']-dropped['p2']),
                    frozen_target_angular_error=max(abs(full[k]-frozen[k]) for k in ('p1','p2')),
                    U_square_contribution=u_square,
                    checks=dict(all_five_moment_changes=moment_error < mp.mpf('1e-55'),
                                original_source_ODEs=source_error < mp.mpf('1e-55'),
                                fifth_pressure_moment_required=abs(full['p2']-dropped['p2']) > mp.mpf('1e-7'),
                                target_angular_derivatives_required=max(abs(full[k]-frozen[k]) for k in ('p1','p2')) > mp.mpf('1e-6'),
                                U_square_term_required=u_square > mp.mpf('1e-9')))
