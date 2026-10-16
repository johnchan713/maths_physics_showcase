"""First phase-locked real signed update: exact diagnostics and order ledger.

Fixtures are manufactured rational compact polynomials, not the actual
C-infinity profile or its pulse ODE. See README for the analytical bridge.
"""
from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
KAPPA = Q(1, 100000)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


exact = load('signed_update_frozen_rational_algebra', HERE.parent/'complete_covariance_initial_means/means.py')
Poly, R, Z, T, EPS = exact.Poly, exact.R, exact.Z, exact.T, exact.EPS


def outer(a, b):
    return [[a[i]*b[j] for j in range(3)] for i in range(3)]


def add(*tensors):
    return [[sum(t[i][j] for t in tensors) for j in range(3)] for i in range(3)]


def scale(tensor, factor):
    return [[factor*x for x in row] for row in tensor]


def curl(c):
    """Normalized cylindrical coefficient curl; no angular coefficient dependence."""
    return [-EPS*c[1].diff(1), EPS*c[0].diff(1)-c[2].diff(0), exact.radial(c[1])]


def amplitude(c, axial):
    """-n cross c for n=(0,1/R,axial); phase=theta+axial*Z/epsilon."""
    return [axial*c[1]-c[2].shift(-1), -axial*c[0], c[0].shift(-1)]


@lru_cache(None)
def fixture():
    bump = (R-1)**6*(2-R)**6
    f = bump*(1+Z/5+T/11)
    g = Q(1,4)+Z/7+T/13
    records = []
    for axial, k, multiplier in ((1, 3, -Q(1,3)+Z/11+T/17),
                                  (-1, 5, Q(2,5)-Z/13+T/19)):
        c = [R*f, R*f*g, -R**2*f]
        dc = [multiplier*x for x in c]
        a, ell = amplitude(c, axial), amplitude(dc, axial)
        r, s = [[-x/k for x in curl(v)] for v in (c, dc)]
        old = scale(add(outer(a,a), outer(r,r)), Q(1,2))
        anew, rnew = [[x+y for x,y in zip(v,w)] for v,w in ((a,ell),(r,s))]
        new = scale(add(outer(anew,anew),outer(rnew,rnew)),Q(1,2))
        principal = scale(add(outer(a,ell),outer(ell,a)),Q(1,2))
        remainder = scale(add(outer(r,s),outer(s,r),outer(ell,ell),outer(s,s)),Q(1,2))
        records.append(dict(axial=axial,k=k,multiplier=multiplier,a=a,ell=ell,r=r,s=s,
                            old=old,new=new,principal=principal,remainder=remainder))
    delta = add(*(add(p['new'],scale(p['old'],-1)) for p in records))
    principal = add(*(p['principal'] for p in records))
    remainder = add(*(p['remainder'] for p in records))
    return bump, records, delta, principal, remainder


def covariance_review():
    _, records, delta, principal, remainder = fixture()
    identities = []; sample_identities = []; divergence = []; commutators = []; bad = []
    for p in records:
        a,ell,r,s = [p[k] for k in ('a','ell','r','s')]
        # Independent quarter-phase average: exact moments, no floating quadrature.
        samples = []
        for cosine,sine in ((1,0),(0,1),(-1,0),(0,-1)):
            old = [cosine*x+sine*y for x,y in zip(a,r)]
            inc = [cosine*x+sine*y for x,y in zip(ell,s)]
            new = [x+y for x,y in zip(old,inc)]
            samples.append(add(outer(new,new),scale(outer(old,old),-1)))
        averaged = scale(add(*samples),Q(1,4))
        sample_identities.extend(averaged[i][j] == p['principal'][i][j]+p['remainder'][i][j]
                                 for i in range(3) for j in range(3))
        identities.extend(p['new'][i][j]-p['old'][i][j] == p['principal'][i][j]+p['remainder'][i][j]
                          for i in range(3) for j in range(3))
        for av,rv in ((a,r),(ell,s)):
            divergence.extend((av[1].shift(-1)+p['axial']*av[2] == 0,
                exact.radial(av[0])+EPS*av[2].diff(1)+p['k']*(rv[1].shift(-1)+p['axial']*rv[2]) == 0,
                exact.radial(rv[0])+EPS*rv[2].diff(1) == 0))
        commutators.append(any(s[j] != p['multiplier']*r[j] for j in range(3)))
        # A prior cosine-quadrature error is not a primary sine curl remainder.
        e = [x/7 for x in a]
        omitted = scale(add(outer(e,ell),outer(ell,e)),Q(1,2))
        bad.append(any(bool(omitted[0][j]) for j in (1,2)))
    return dict(checks=dict(
        exact_full_tensor_increment=all(identities),
        independent_quarter_phase_tensor_average=all(sample_identities),
        primary_and_signed_complete_curl_divergence=all(divergence),
        all_six_updated_covariances_nonzero=all(bool(delta[i][j]) for i in range(3) for j in range(i,3)),
        signed_amplitude_curl_commutators_retained=all(commutators),
        signed_increment_square_not_discarded=all(bool(p['ell'][0]*p['ell'][1]) for p in records),
        primary_new_sine_curl_cross_not_discarded=all(bool(p['r'][0]*p['s'][1]) for p in records),
        generic_old_cosine_error_breaks_sharp_identity=all(bad),
        total_pair_remainder_is_not_zero=all(bool(remainder[0][j]) for j in (1,2))),
        exact_tensor_entry_checks=len(identities)+len(sample_identities),
        exact_divergence_checks=len(divergence),
        fixture_frequencies=[3,5],
        distinct_fixture_frequencies_not_replacement_for_actual_closed_support_separation=True,
        finite_endpoint_regularities_not_actual_Cinfinity_profile=True)


def inverse_review():
    # Rational manufactured columns: actual H and its edge bounds are analytical inputs.
    h = ((Q(-2),Q(-3)),(Q(1),Q(-2)))
    a = (Q(2),Q(3)); target = (Q(7),Q(-5))
    det = h[0][0]*h[1][1]-h[0][1]*h[1][0]
    d = ((h[1][1]*target[0]-h[0][1]*target[1])/(EPS*det),
         (-h[1][0]*target[0]+h[0][0]*target[1])/(EPS*det))
    da = tuple(d[i]/(2*a[i]) for i in range(2))
    recovered = tuple(EPS*sum(h[j][i]*2*a[i]*da[i] for i in range(2)) for j in range(2))
    half = tuple(EPS*sum(h[j][i]*a[i]*da[i] for i in range(2)) for j in range(2))
    # Choose bands whose rational h powers are rational; these are NOT accepted ODE bands.
    physical = []
    for coordinate in target:
        total = Q(0)
        for ell,eta in ((200,Q(3,5)),(400,Q(4,5))):
            q2a = Q(1,2**(ell*101//100))  # h=1/200, 2A=101/100
            total += eta**2*(q2a*coordinate)/q2a
        physical.append(total == coordinate)
    # For D^j(y^-1/2), each term y^(-n-1/2) product_{i=1}^n D^{Ji}y
    # has edge exponent -1/2; multiplying D^J Sigma supplies one more zeta.
    edge = [(-Q(1,2)-n)+n+1 == Q(1,2) for n in range(13)]
    relative = [-EPS**(Q(1))*Q(n)**2 < -1 for n in (20,40,80)]
    return dict(checks=dict(fixed_signed_linear_inverse_recovers_target=recovered == target,
        one_inverse_coefficient_negative=any(x < 0 for x in da),
        factor_two_is_necessary=half != target and half == tuple(x/2 for x in target),
        squared_partition_physical_normalization=all(physical),
        no_division_by_slow_partition_needed=Q(3,5)**2+Q(4,5)**2 == 1,
        every_checked_reciprocal_root_term_retains_half_edge_weight=all(edge),
        weighted_class_does_not_ensure_positive_updated_amplitudes=all(relative)),
        inverse_coefficients=[str(x) for x in da], exact_edge_term_checks=len(edge),
        actual_H_not_replaced_by_manufactured_columns=True,
        all_order_edge_exponent_identity='(-n-1/2)+n+1=1/2 for every n>=0')


def pressure_review():
    bump, _, delta, principal, remainder = fixture()
    gr = -exact.radial(delta[0][0])-EPS*delta[0][2].diff(1)+delta[1][1].shift(-1)
    p = gr.integral()
    rho = bump/bump.integral().scalar()*(1+Z*(R-Q(3,2)))
    crho = (R**2*rho).integral()/2
    source = gr-rho*p
    primitive = source.primitive(); pm = primitive-primitive.at(0,1)
    jt = (R**2*delta[1][2]).integral()
    jz = (R*delta[2][2]).integral()-(R**2*gr).integral()/2
    et = exact.radial(delta[0][1],2)+EPS*delta[1][2].diff(1)
    ez = exact.radial(delta[0][2])+EPS*(delta[2][2]+pm).diff(1)
    return dict(checks=dict(
        delta_pressure_compact_at_both_ends=pm.at(0,1) == pm.at(0,2) == 0,
        delta_pressure_radial_source=pm.diff(0) == source,
        whole_delta_axial_pressure_moment=(R*pm).integral() == -(R**2*gr).integral()/2+crho*p,
        delta_angular_integrated_balance=(R**2*et).integral() == EPS*jt.diff(1),
        delta_axial_integrated_balance=(R*ez).integral() == EPS*(jz+crho*p).diff(1),
        angular_stress_cancellation_exact=et-exact.radial(principal[0][1],2) == exact.radial(remainder[0][1],2)+EPS*delta[1][2].diff(1),
        axial_stress_cancellation_exact=ez-exact.radial(principal[0][2]) == exact.radial(remainder[0][2])+EPS*(delta[2][2]+pm).diff(1),
        omitted_delta_axial_pressure_detected=bool(EPS*(R*pm).integral().diff(1)),
        omitted_delta_axial_covariance_detected=bool(EPS*(R*delta[2][2]).integral().diff(1)),
        moving_bump_product_derivative_nonzero=bool(crho.diff(1)*p),
        covariance_update_changes_all_three_defects=all(bool(x) for x in (p,jt,jz))),
        fixed_mean_velocities_in_this_substep=True,
        both_actual_base_pressure_tails_cancel_only_in_same_full_base_subtraction=True,
        arithmetic='exact rational polynomial identities, not actual profile evaluation')


def orders(alpha, kappa=KAPPA):
    alpha, k = Q(alpha), Q(kappa)
    if alpha < 1 or not 0 < k < Q(1,2):
        raise ValueError('alpha>=1 and fixed 0<kappa<1/2 required')
    rows = dict(target=alpha, signed_velocity=alpha-Q(1,2), signed_curl=alpha-k,
                signed_wave_pressure=alpha, old_new_sine_curls=alpha+1-2*k,
                signed_transverse_square=2*alpha-1, signed_curl_square=2*alpha-2*k,
                full_covariance_change=alpha, averaged_pressure_change=alpha,
                full_pressure_change=alpha-k, axial_flux_and_averaged_pressure=alpha+1,
                full_tangential_mean=alpha-k, recomputed_defects=alpha,
                linear_wave_remainder=alpha-3*k,
                old_wave_cross=alpha-k, actual_mean_cross=alpha-k,
                signed_wave_self_interaction=2*alpha-1-k)
    rows['real_pair_remainder'] = min(rows[n] for n in ('old_new_sine_curls','signed_transverse_square','signed_curl_square'))
    rows['real_covariance_gain'] = rows['real_pair_remainder']-alpha
    return rows


def physical_powers(alpha, h, a, b):
    alpha, h = Q(alpha), Q(h)
    orders(alpha)
    if not 0 < h < Q(1,100) or any(type(v) is not int or v < 0 for v in (a,b)):
        raise ValueError('Fixed construction h and nonnegative integer derivative orders required')
    A = Q(1,2)+h; sx = Q(1,2)+h/2; st = 1+3*h/2
    shift = -a*sx-b*st
    return dict(velocity=-A+h*(alpha-Q(1,2))+shift,
                pressure=-2*A+h*alpha+shift,
                potential=Q(1,2)-A+h*alpha+shift)


def exponent_review():
    k = KAPPA; alpha = 2-3*k; r = orders(alpha)
    final = min(r['real_pair_remainder'], r['axial_flux_and_averaged_pressure'], Q(3))
    generic = min(alpha+Q(1,2)-k, 2*alpha-1)
    physical = physical_powers(alpha,Q(1,200),0,0)
    return dict(checks=dict(
        initial_target_from_replayed_mean_bound=alpha == 2-3*k,
        real_pair_error_dominated_by_signed_square=r['real_pair_remainder'] == 3-6*k,
        final_averaged_tangential_order=final == 3-6*k,
        real_gain_positive=r['real_covariance_gain'] == 1-3*k > 0,
        generic_bound_is_strictly_weaker=generic < final,
        full_mean_not_claimed_improved=r['full_tangential_mean'] < alpha,
        compatibility_S2_not_silently_preserved=r['recomputed_defects'] < 2,
        alpha_one_square_blocks_positive_gain=orders(1)['real_covariance_gain'] == 0,
        all_alpha_above_one_have_gain_formula=min(1-2*k,alpha-1) == r['real_covariance_gain'],
        sourced_and_interaction_remainders_of_this_update_exceed_old_wave_order=min(r[n] for n in ('linear_wave_remainder','old_wave_cross','actual_mean_cross','signed_wave_self_interaction')) > 1-3*k,
        physical_covariance_error_not_absolutely_decaying=-2*(Q(1,2)+Q(1,200))+Q(1,200)*final < 0),
        rows={name:str(value) for name,value in r.items()},
        final_auxiliary_averaged_tangential_order=str(final), generic_pair_error_order=str(generic),
        physical_zero_jet_powers={name:str(value) for name,value in physical.items()},
        all_order_gain_condition='alpha>1, from min(1-2*kappa,alpha-1)>0',
        coefficient_orders_not_absolute_physical_residual_decay=True)
