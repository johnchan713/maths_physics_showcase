"""Actual phase/frame FORMULAS, tested with manufactured smooth background data.

No function here evaluates the actual huge-frequency global leading profile
or the infinite background sum. The continuum proof uses their inherited
normalized bounds; the numerical data only test this interface's algebra.
"""
from fractions import Fraction
import importlib.util
from pathlib import Path
import mpmath as mp

PARENT = Path(__file__).resolve().parent.parent/'physical_wave_interface'
_spec = importlib.util.spec_from_file_location('frozen_wave_taylor3', PARENT/'jets.py')
_jets = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_jets)
Jet, variables = _jets.Jet, _jets.variables

H = Fraction(1, 200)
KAPPA = Fraction(1, 100000)
REP = ('1.4', '.12', '.62')
U = mp.mpf('1.5')


def rational(value):
    return mp.mpf(value.numerator)/value.denominator if isinstance(value, Fraction) else mp.mpf(value)


def nearest_nonzero_integer(x):
    answer = int(mp.floor(x+mp.mpf('.5')))
    return answer if answer else (1 if x >= 0 else -1)


def band(ell, h=H, r0='.04'):
    if type(ell) is not int or ell < 16:
        raise ValueError('A fixed integer band >=16 is required')
    h_input = h
    h = rational(h)
    if not 0 < h < mp.mpf('.01') or not 0 < mp.mpf(r0) < mp.mpf('.1'):
        raise ValueError('Construction h and a small fixed rectangle radius required')
    tg, lg = 4+mp.sqrt(2), 4-mp.sqrt(2)
    rho = mp.log(lg)/mp.log(tg)
    scale, logq = mp.mpf(ell)**2, -ell*mp.log(2)
    # Avoid a one-unit ceil error at exact powers of two.
    exponent = rational(h_input*ell/2) if isinstance(h_input, Fraction) else h*ell/2
    eps = mp.power(2, -2*exponent)
    if mp.floor(exponent) == exponent:
        k = 1 << int(exponent)
    else:
        k = int(mp.ceil(mp.power(2, exponent)))
    level = int(mp.floor((-(1+h)*logq-mp.log(scale))/mp.log(tg)))
    ci = mp.exp(level*mp.log(tg)+(1+h)*logq)
    dr = 2*((1+h)*rho-h*rational(KAPPA))
    mi = mp.exp(level*mp.log(lg)+dr*logq/2)
    return dict(ell=ell, h=h, scale=scale, epsilon=eps, k=k, level=level,
                ci=ci, length=2*mp.mpf(r0)/ci, mi=mi, rho=rho, tg=tg, lg=lg, dr=dr)


def qstar(z, time, h):
    """The bounded implicit coordinate q/Q, on the positive chart branch."""
    q = time+z*z
    for _ in range(12):
        q -= (q-z*z*q**(2*h)-time)/(1-2*h*z*z*q**(2*h-1))
    return q


def background(r, z, time, epsilon, h):
    """Smooth profile fixture plus O(eps^2) tangential perturbations."""
    q = qstar(z, time, h)
    d, a = mp.mpf('.5')-h, mp.mpf('.5')+h
    eta = z/q**d
    factor = q**d*(1+mp.mpf('.1')*eta+mp.mpf('.03')*eta*eta)
    leading_f = factor/r**3
    leading_g = -mp.mpf('1.5')*factor/r**2
    f = leading_f+epsilon**2*(mp.mpf('.2')*r*r+mp.mpf('.1')*z*time)
    g = leading_g+epsilon**2*(mp.mpf('.13')*r*z+mp.mpf('.07')*time*r*r)
    b = epsilon*(mp.mpf('.17')*r+mp.mpf('.08')*z*time)
    return b, f, g, leading_f, leading_g


def label(data, sign):
    if sign not in (1, -1):
        raise ValueError('A fixed pulse sign is required')
    r0, z0, t0 = map(mp.mpf, REP)
    h, eps = data['h'], data['epsilon']
    leading_f = background(r0, z0, t0, 0, h)[1]
    gr = mp.matrix([r0*mp.diff(lambda r: background(r, z0, t0, 0, h)[1], r0),
                    mp.diff(lambda r: background(r, z0, t0, 0, h)[2], r0)])
    magnitude = mp.norm(gr)
    nt = gr/magnitude
    kt = mp.matrix([-nt[1], nt[0]])
    lam0 = mp.sqrt(-2*leading_f*nt[0]*(2*leading_f*nt[0]+magnitude))
    c0 = lam0/(2*leading_f*nt[0])
    bs = mp.sqrt(lam0/(eps*data['k']**2*(1+U**2)**mp.mpf('1.5')))
    tan = bs*(kt-sign*U*gr/(data['length']*magnitude**2))
    ptilde = r0*tan[0]
    angular = nearest_nonzero_integer(data['k']*ptilde)
    return dict(sign=sign, f0=leading_f, shear=gr, nt=nt, kt=kt, lambda0=lam0,
                c0=c0, bs=bs, p=mp.mpf(angular)/data['k'], pz=tan[1],
                angular_integer=angular, ptilde=ptilde, x0=sign*bs*U/2)


def local_coefficients(data, lab, point):
    r, z, time = map(mp.mpf, point)
    h, eps = data['h'], data['epsilon']
    b, f, g = background(r, z, time, eps, h)[:3]
    fr = mp.diff(lambda x: background(x, z, time, eps, h)[1], r)
    gr = mp.diff(lambda x: background(x, z, time, eps, h)[2], r)
    fz = mp.diff(lambda x: background(r, x, time, eps, h)[1], z)
    gz = mp.diff(lambda x: background(r, x, time, eps, h)[2], z)
    ft = mp.diff(lambda x: background(r, z, x, eps, h)[1], time)
    gt = mp.diff(lambda x: background(r, z, x, eps, h)[2], time)
    p, pz = lab['p'], lab['pz']
    normal0 = mp.matrix([lab['x0'], p/r, pz])
    normal_prime = mp.matrix([-(p*fr+pz*gr), 0, -eps*(p*fz+pz*gz)])
    kmat = mp.matrix([[0, -2*f, 0], [2*f+r*fr, 0, 0], [gr, 0, 0]])
    return dict(normal0=normal0, normal_prime=normal_prime, kmat=kmat, b=b,
                transport_derivative=(p*ft+pz*gt)-g*(p*fz+pz*gz))


def outer(a, b):
    return mp.matrix([[a[i]*b[j] for j in range(len(b))] for i in range(len(a))])


def frame(data, lab, coefficients, v):
    n = coefficients['normal0']+v*coefficients['normal_prime']
    np_ = coefficients['normal_prime']
    tan, tanprime = mp.matrix([n[1], n[2]]), mp.matrix([np_[1], np_[2]])
    mag = mp.norm(tan)
    magprime = (tan.T*tanprime)[0]/mag
    kt = tan/mag
    ktprime = tanprime/mag-kt*magprime/mag
    nt, ntprime = mp.matrix([kt[1], -kt[0]]), mp.matrix([ktprime[1], -ktprime[0]])
    sa = n[0]/mag
    saprime = np_[0]/mag-sa*magprime/mag
    root = mp.sqrt(1+sa*sa)
    rootprime = sa*saprime/root
    b = mp.matrix(3, 2)
    bp = mp.matrix(3, 2)
    for col, sign in ((0, 1), (1, -1)):
        b[0, col] = 1
        for j in range(2):
            b[j+1, col] = -sa*kt[j]+sign*lab['c0']*root*nt[j]
            bp[j+1, col] = -saprime*kt[j]-sa*ktprime[j]+sign*lab['c0']*(rootprime*nt[j]+root*ntprime[j])
    left = mp.matrix([[mp.mpf('.5'), nt[0]/(2*lab['c0']*root), nt[1]/(2*lab['c0']*root)],
                      [mp.mpf('.5'), -nt[0]/(2*lab['c0']*root), -nt[1]/(2*lab['c0']*root)]])
    k = coefficients['kmat']
    projected = -k+outer(n, n.T*k-np_.T)/(n.T*n)[0]
    s = lab['sign']*U*(mp.mpf('.5')+v/data['length'])
    growth = lab['lambda0']/mp.sqrt(1+s*s)
    damping = data['epsilon']*data['k']**2*(n.T*n)[0]
    dref = lab['lambda0']*(1+s*s)/(1+U*U)**mp.mpf('1.5')
    error = left*(projected*b-bp)-mp.diag([growth, -growth])
    return dict(normal=n, normal_prime=np_, basis=b, basis_prime=bp, left=left,
                projected=projected, error=error, growth=growth, damping=damping,
                damping_error=damping-dref, reference_normal=lab['bs']*mp.matrix([s, *lab['kt']]))


def qstar_jet(z, time, h):
    val = float(qstar(mp.mpf(z.value), mp.mpf(time.value), mp.mpf(h)))
    q = Jet(val)
    for _ in range(4):
        q = q-(q-z*z*q**(2*h)-time)/(1-2*h*z*z*q**(2*h-1))
    return q


def background_jet(r, z, time, epsilon, h):
    q = qstar_jet(z, time, h)
    d = .5-h
    eta = z/q**d
    factor = q**d*(1+.1*eta+.03*eta*eta)
    f = factor/r**3+epsilon**2*(.2*r*r+.1*z*time)
    g = -1.5*factor/r**2+epsilon**2*(.13*r*z+.07*time*r*r)
    b = epsilon*(.17*r+.08*z*time)
    return b, f, g
