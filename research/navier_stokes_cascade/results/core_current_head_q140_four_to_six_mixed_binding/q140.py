"""Bind four actual PREMOD mixed fields to Q139's six scalar jets."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util

HERE = Path(__file__).resolve().parent

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, HERE.parent/path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

previous = load('q140_previous', 'core_current_head_q139_mixed_primitive_binding/q139.py')
cutoff = load('q140_cutoff', 'core_current_head_q128_real_cutoff_signed_root/q128.py')
exact = previous.exact
order = previous.order
lg = previous.lg

IDENTITY = 'same-actual-four-PREMOD-fields-full-p2-selected-cutoff-and-signed-root'
KINDS = {'a', 'ts', 'p', 'mu', 'E', 'v'}

def four_to_six(H, Ay, Be, Acap, Bcap, afloor, delta, G, Eexp, d0, identity):
    """Close the selected cutoff/root algebra under a four-field mixed input.

    The input is the summed raw mixed bound for (a,b,p2,E), divided by
    i!k!.  The output is a common six-scalar raw bound in Q139's convention.
    """
    vals = tuple(map(exact, (H, Ay, Be, Acap, Bcap, afloor, delta, G, Eexp, d0)))
    H, Ay, Be, Acap, Bcap, afloor, delta, G, Eexp, d0 = vals
    if identity != IDENTITY:
        raise ValueError('Same actual PREMOD fields, full p2, selected cutoff and signed root required')
    if min(H, Ay, Be, Acap, afloor, delta) <= 0 or Bcap < 0 or G < 2:
        raise ValueError('Positive mixed caps/floor/gap and original real value cap required')

    # Analytic maps (a,b)->(-b/a, a+b^2/a) on the same real-center disk.
    r = min(Q(1, 4), afloor/4)
    Tcap = 2*(Bcap+1)/afloor
    Vcap = Acap+1+2*(Bcap+1)**2/afloor
    scalar_factor = 8*(1+H/r)
    A0, B0 = scalar_factor*Ay, scalar_factor*Be

    # The original real cutoff.  Only its one-variable Gevrey-3 derivative
    # bound is used; no analytic continuation through a flat edge is claimed.
    gap = delta/4
    Hstar = 2+delta/2
    q = min(Q(1, 4), afloor/4, gap/4)
    Hav = H+Vcap
    Aav, Bav = max(Ay, A0), max(Be, B0)
    rad_factor = 8*(1+Hav/q)
    Arad, Brad = rad_factor*Aav, rad_factor*Bav
    Crad = 2*max(Q(1), Hstar/afloor)
    cutoff_factor = 8*(1+64*Vcap/delta)
    Ak, Bk = cutoff_factor*A0, cutoff_factor*B0
    Ms = 729*Crad
    As, Bs = 2*max(Ak, Arad), 2*max(Bk, Brad)

    # Same Q65 signed removable inverse, including zero p and zero variance.
    inverse = cutoff.signed.radii(G, Eexp, d0)
    Hps = H+Ms
    Aps, Bps = max(Ay, As), max(Be, Bs)
    mu_factor = 8*(1+(2/inverse['sigma'])*Hps)
    Amu, Bmu = mu_factor*Aps, mu_factor*Bps
    Mmu = max(G, inverse['a'])

    # v is the augmented variance used by Q130/Q139, not the unaugmented v0.
    Mv = Vcap+H*Ms**2
    Av = max(A0, 2*max(Ay, 2*As))
    Bv = max(B0, 2*max(Be, 2*Bs))

    Hsix = 3*H+Tcap+Mmu+Mv
    Asix = max(Ay, A0, Amu, Av)
    Bsix = max(Be, B0, Bmu, Bv)
    return dict(
        stage='actual-four-PREMOD-to-six-mixed-scalars',
        H=Hsix, Ay=Asix, Be=Bsix,
        components=dict(a=H, ts=Tcap, p=H, mu=Mmu, E=H, v=Mv),
        unaugmented_variance_cap=Vcap, amplitude_prefactor=Ms,
        amplitude_y_base=As, amplitude_eta_base=Bs,
        signed_root_radius=inverse['sigma'], signed_root_prefactor=Mmu,
        scalar_disk=r, radicand_disk=q,
        full_pressure_identity_required=True,
        flat_cutoff_slow_analytic_extension_claimed=False,
        actual_regional_four_field_mixed_bound_certified=False)

def bind(recipe, disk, V0, N, identity):
    if identity != IDENTITY or recipe.get('stage') != 'actual-four-PREMOD-to-six-mixed-scalars':
        raise ValueError('Same four-to-six recipe required')
    primitive = previous.bind(recipe['H'], recipe['Ay'], recipe['Be'], disk, previous.IDENTITY)
    pullback = previous.pullback(primitive, V0, N, previous.IDENTITY)
    return dict(stage='conditional-four-field-to-Ny-pullback', primitive=primitive,
                pullback=pullback,
                actual_regional_four_field_mixed_bound_certified=False,
                actual_complete_profile_cone_or_scale_certified=False)

def bound(recipe, kind, i, k, binary=False):
    order(i); order(k)
    if recipe.get('stage') != 'actual-four-PREMOD-to-six-mixed-scalars' or kind not in KINDS:
        raise ValueError('Six-scalar mixed recipe and component required')
    P, A, B = recipe['H'], recipe['Ay'], recipe['Be']
    if binary:
        value = lg(max(Q(1), P))+i*lg(A)+k*lg(B)
        value += 3*i*lg(Q(i)) if i else 0
        value += 3*k*lg(Q(k)) if k else 0
    else:
        if max(i, k) > 12:
            raise ValueError('Use binary mixed budget at high order')
        value = P*A**i*B**k*(factorial(i)*factorial(k))**3
    return dict(component=kind, y_order=i, eta_order=k, upper=value,
                binary=binary, actual_regional_input_certified=False)

def q72_log2_lower(t, c, m):
    """Strict lower bound for log2(exp(c*T_m)) from Q72's last T_m term."""
    t, c = exact(t), exact(c); order(m)
    if t <= 0 or c <= 0 or m == 0:
        raise ValueError('Positive inherited width/coefficient and positive order required')
    last = 162*Q(3260)**m*factorial(m)**2/t**m
    return (c*last).numerator//(c*last).denominator

def gevrey3_log2_budget(H, A, m):
    H, A = exact(H), exact(A); order(m)
    if min(H, A) < 1:
        raise ValueError('Positive comparison constants at least one required')
    return lg(H)+m*lg(A)+(3*m*lg(Q(m)) if m else 0)

def record():
    return dict(
        status='four-PREMOD-mixed-fields-to-six-scalar-binding-and-short-taper-obstruction',
        input='sum over a,b,p2,E of raw mixed jets/(i!k!) <= H Ay^i Be^k(i!k!)^2; p2 retains Pi and Pi_eta',
        composition='finite two-variable Gevrey-3 composition: outer Q C^n(n!)^2 gives bases 8(Ay,Be)(1+CH)',
        cutoff='same Q66 flat cutoff, same positive radicand gap and same Q65 signed removable inverse; no slow analytic extension',
        output='six actual scalars a,-b/a,p2,mu,E,a+b^2/a+a*s^2 satisfy Q139 with explicit Hsix,Asix,Bsix',
        correction='Q72 finite-order cap exp(c*T_m), T_m>=162*3260^m*(m!)^2/t^m, cannot itself be compressed into fixed Gevrey-3 constants',
        obstruction_scope='obstruction is in the inherited cap, not a non-Gevrey claim about the actual short-taper profile',
        next='derive a value-based differentiated short-taper recurrence, then assemble actual four-field regional constants and apply Q139/Q138',
        actual_complete_profile_pulse_covariance_residual_or_blowup_certified=False,
        historical_frozen_records_modified=False)
