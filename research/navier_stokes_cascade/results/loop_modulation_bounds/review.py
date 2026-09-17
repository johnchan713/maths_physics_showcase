"""Independent differentiation checks on MANUFACTURED inputs, not the actual core."""
from functools import lru_cache
import importlib.util
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    'modulation_previous_loop', HERE.parent/'stress_realization_audit'/'loop.py')
loop = importlib.util.module_from_spec(spec)
spec.loader.exec_module(loop)


def signed_F(mu, p, d0):
    """The signed square root is smooth at zero; sqrt(V) alone has a cusp."""
    if mu == 0:
        return mp.mpf(0)
    return mp.sign(mu)*mp.sqrt(loop.variance(abs(mu), p, d0))


def solve_mu(s, p, d0):
    """A diagnostic bisection, deliberately independent of implicit differentiation."""
    if s == 0:
        return mp.mpf(0)
    if p == 0:
        return mp.sqrt(2)*s/d0
    target = s*s
    lo, hi = mp.mpf(0), max(mp.mpf(1), 2*mp.sqrt(2)*abs(s)/d0)
    for _ in range(30):
        if loop.variance(hi, p, d0) > target:
            break
        hi *= 2
    else:
        raise ArithmeticError('Manufactured root did not bracket')
    for _ in range(240):
        mid = (lo+hi)/2
        if loop.variance(mid, p, d0) < target:
            lo = mid
        else:
            hi = mid
    return mp.sign(s)*(lo+hi)/2


@lru_cache(maxsize=None)
def root_review(digits=70):
    with mp.workdps(digits):
        d0, x = mp.mpf('.5'), mp.mpf('.2')
        steps = list(map(mp.mpf, ('1e-5','1e-7','1e-9')))
        records = []
        for p0 in map(mp.mpf, ('0','.3','-.3')):
            s = lambda q: mp.mpf('.8')+mp.mpf('.05')*q
            p = lambda q: p0+mp.mpf('.2')*(q-x)
            mu = solve_mu(s(x),p(x),d0)
            fmu = mp.diff(lambda u: signed_F(u,p(x),d0),mu)
            fp = mp.diff(lambda z: signed_F(mu,z,d0),p(x))
            predicted = (mp.mpf('.05')-fp*mp.mpf('.2'))/fmu
            errors = []
            for h in steps:
                numerical = (solve_mu(s(x+h),p(x+h),d0)
                             -solve_mu(s(x-h),p(x-h),d0))/(2*h)
                errors.append(abs(numerical-predicted))
            records.append(dict(p=p0, mu=mu, derivative=predicted,
                                errors=errors,
                                omitted_pressure_derivative_error=abs(fp*mp.mpf('.2')/fmu),
                                fmu=fmu,
                                fmu_lower=d0/32*mp.exp(-4*abs(mu*p(x))),
                                fp_bound=3*d0*mu**2*mp.exp(2*abs(mu*p(x)))))
        fmu_zero = mp.diff(lambda u: signed_F(u,mp.mpf('.3'),d0),mp.mpf(0))
        fp_zero = mp.diff(lambda p: signed_F(mp.mpf(0),p,d0),mp.mpf('.3'))
        return dict(manufactured=True, records=records,
                    checks=dict(
                        implicit_derivative_refines=all(r['errors'][2] < r['errors'][0]/10**6
                                                       for r in records),
                        implicit_derivative_agrees=all(r['errors'][2] < mp.mpf('1e-16') for r in records),
                        root_denominator_bound=all(r['fmu'] >= r['fmu_lower'] for r in records),
                        pressure_derivative_bound=all(abs(mp.diff(
                            lambda z: signed_F(r['mu'],z,d0),r['p'])) <= r['fp_bound'] for r in records),
                        omitted_pressure_term_detected=all(r['omitted_pressure_derivative_error'] > mp.mpf('.001')
                                                          for r in records if r['p'] != 0),
                        zero_root_limit=abs(fmu_zero-d0/mp.sqrt(2)) < mp.mpf('1e-60') and fp_zero == 0))


def fixture(y, eta):
    """Positive swirl with its exact shear a; U=bs=ts=p2=0 in this fixture.

    Integrating (1-a)/2 in y fixes E, so the modulation check starts from
    compatible fields rather than choosing their shears independently.
    """
    a = mp.mpf('.8')+mp.mpf('.05')*eta+mp.mpf('.03')*mp.sin(y)
    E = (mp.mpf('1.2')+mp.mpf('.1')*eta)*mp.exp(
        (mp.mpf('.1')-mp.mpf('.025')*eta)*y+mp.mpf('.015')*mp.cos(y))
    v = mp.mpf('2.1')+mp.mpf('.02')*eta+mp.mpf('.01')*mp.sin(y)
    return a,v,E


def theta_for_phase(y, eta, phi):
    """Invert the explicit strictly increasing lift on the circle."""
    a,v,E = fixture(y,eta)
    turns = mp.floor(phi)
    phi0 = phi-turns
    target = 2*mp.pi*phi0
    coefficient = (v-a)/(2*v)
    theta = mp.findroot(lambda t: t-coefficient*mp.sin(2*t)-target,
                        target, df=lambda t: 1-2*coefficient*mp.cos(2*t),
                        solver='newton', maxsteps=100)
    return theta+2*mp.pi*turns


def primitive(y, eta, phi, key):
    a,v,E = fixture(y,eta)
    theta = theta_for_phase(y,eta,phi)
    return loop.zero_pressure_primitives(theta,a,v,E)[key]


@lru_cache(maxsize=None)
def phase_review(digits=70):
    with mp.workdps(digits):
        y,eta = mp.mpf('.37'),mp.mpf('.2')
        steps = list(map(mp.mpf, ('1e-5','1e-7','1e-9')))
        records = []
        for phi in map(mp.mpf, ('.19','.61','.87')):
            theta = theta_for_phase(y,eta,phi)
            a,v,E = fixture(y,eta)
            data = loop.zero_pressure_primitives(theta,a,v,E)
            for axis in (0,1):
                def at(q,angle,key):
                    args = (q,eta) if axis == 0 else (y,q)
                    return loop.zero_pressure_primitives(angle,*fixture(*args))[key]
                x = (y,eta)[axis]
                phi_x = mp.diff(lambda q: at(q,theta,'phi'),x)
                theta_x = -phi_x/data['phase_derivative']
                for key in ('A','B'):
                    frozen = mp.diff(lambda q: at(q,theta,key),x)
                    angular = mp.diff(lambda t: at(x,t,key),theta)
                    expected = frozen+angular*theta_x
                    def actual(q):
                        args = (q,eta) if axis == 0 else (y,q)
                        return primitive(*args,phi,key)
                    errors = [abs((actual(x+h)-actual(x-h))/(2*h)-expected) for h in steps]
                    records.append(dict(phase=phi,axis=axis,primitive=key,errors=errors,
                                        frozen_phase_error=abs(angular*theta_x)))
        return dict(manufactured=True,records=records,
                    checks=dict(
                        inverse_phase_chain_rule_refines=all(r['errors'][2] < r['errors'][0]/10**6 for r in records),
                        inverse_phase_chain_rule_agrees=all(r['errors'][2] < mp.mpf('1e-18') for r in records),
                        frozen_inverse_phase_rejected=max(r['frozen_phase_error'] for r in records) > mp.mpf('1e-4'),
                        both_slow_coordinates_checked={r['axis'] for r in records} == {0,1}))


@lru_cache(maxsize=None)
def modulation_review(digits=70):
    with mp.workdps(digits):
        y,eta = mp.mpf('.37'),mp.mpf('.2')
        records = []
        for n in (32,64,128):
            phi = n*y
            a,v,E = fixture(y,eta)
            A,B = [primitive(y,eta,phi,k) for k in ('A','B')]
            Ay,By = [mp.diff(lambda q: primitive(q,eta,phi,k),y) for k in ('A','B')]
            al,bl = [primitive(y,eta,phi,k) for k in ('a','b')]
            def fields(q,e):
                base = fixture(q,e)[2]
                return (base*mp.exp(primitive(q,e,n*q,'A')/n),
                        primitive(q,e,n*q,'B')/n)
            en,un = fields(y,eta)
            ey = mp.diff(lambda q: fields(q,eta)[0],y)
            uy = mp.diff(lambda q: fields(q,eta)[1],y)
            direct_a,direct_b = 1-2*ey/en,2*uy/en
            predicted_a = al-2*Ay/n
            inside = bl+2*By/(n*E)
            predicted_b = mp.exp(-A/n)*inside
            angular_error = abs(mp.diff(lambda e: fields(y,e)[0]-fixture(y,e)[2],eta))
            records.append(dict(N=n,
                exact_shear_identity_error=max(abs(direct_a-predicted_a),abs(direct_b-predicted_b)),
                omitted_exponential_error=abs(direct_b-inside),
                N_times_angular_field_error=n*angular_error,
                raw_radial_field_error=abs(ey-mp.diff(lambda q: fixture(q,eta)[2],y))))
        # An eta-dependent phase destroys the small angular error, even in
        # this simple independent control with value norm at most 1/N.
        wrong_phase_derivative = mp.diff(lambda e: mp.sin(2*mp.pi*128*(y+e))/128,eta)
        return dict(manufactured=True,records=records,
                    checks=dict(
                        exact_shear_identities=all(r['exact_shear_identity_error'] < mp.mpf('1e-55') for r in records),
                        omitted_exponential_rejected=max(r['omitted_exponential_error'] for r in records) > mp.mpf('1e-6'),
                        angular_error_scaled_bounded_in_fixture=max(r['N_times_angular_field_error'] for r in records) < 1,
                        radial_smallness_not_assumed=max(r['raw_radial_field_error'] for r in records) > mp.mpf('.01'),
                        eta_dependent_phase_rejected=abs(wrong_phase_derivative) > 1))
