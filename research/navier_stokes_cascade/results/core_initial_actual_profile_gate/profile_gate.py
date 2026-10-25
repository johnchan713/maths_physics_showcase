"""Exact Fourier gate for the actual first-real covariance profile.

This is an algebraic control.  The continuum input that the selected
covariance has the inherited strict transverse support is recorded by the
checkpoint, rather than inferred from this finite calculation.
"""
from fractions import Fraction as Q


def exponents(kappa):
    kappa=Q(kappa)
    if not 0<kappa<Q(1,4):
        raise ValueError('need 0<kappa<1/4')
    alpha=2-3*kappa
    D=alpha-kappa
    nu=D+1-2*kappa
    return dict(kappa=kappa,alpha=alpha,D=D,nu=nu)


def mode_divisors(k1,k2):
    """Return exact quadratic data for vr=(1,1-sqrt(2)), vt=(sqrt(2)-1,1)."""
    if not isinstance(k1,int) or not isinstance(k2,int) or (k1,k2)==(0,0):
        raise ValueError('a nonzero integer Fourier mode is required')
    a=k1+k2
    pell_norm=a*a-2*k2*k2
    # alpha*tau=(sqrt(2)-1)*pell_norm.  An integer solution of
    # a^2=2*k2^2 is necessarily zero, so the norm cannot vanish here.
    return dict(mode=(k1,k2),a=a,pell_norm=pell_norm,
                radial_divisor_nonzero=pell_norm!=0,
                temporal_divisor_nonzero=pell_norm!=0,
                source_to_viscosity_multiplier_nonzero=pell_norm!=0)


def fourier_injectivity(modes):
    """Check injectivity of L^2 Nt^-1 L on a finite exact mean-zero fixture."""
    clean={tuple(k):Q(v) for k,v in modes.items() if Q(v)}
    if (0,0) in clean:
        raise ValueError('the zero-Haar profile may not contain mode zero')
    rows=[mode_divisors(*k) for k in clean]
    return dict(input_modes=tuple(sorted(clean)),output_modes=tuple(sorted(clean)),
                all_radial_divisors_nonzero=all(r['radial_divisor_nonzero'] for r in rows),
                all_temporal_divisors_nonzero=all(r['temporal_divisor_nonzero'] for r in rows),
                transfer_is_injective=bool(clean) and all(r['source_to_viscosity_multiplier_nonzero'] for r in rows),
                pell_norms={str(r['mode']):r['pell_norm'] for r in rows})


def actual_profile_implication(selected_target_nonzero=True,
                               strict_transverse_support=True,
                               multiplier_auxiliary_independent=True):
    """Logical implication supplied by the inherited actual-profile facts."""
    profile_nonzero=bool(selected_target_nonzero)
    profile_nonconstant=profile_nonzero and bool(strict_transverse_support)
    graph_source_nonzero=profile_nonconstant
    viscosity_leading_coefficient_nonzero=profile_nonconstant
    return dict(
        selected_target_nonzero=bool(selected_target_nonzero),
        multiplier_auxiliary_independent=bool(multiplier_auxiliary_independent),
        physical_profile_auxiliary_independent=False,
        profile_nonzero=profile_nonzero,
        profile_nonconstant=profile_nonconstant,
        graph_source_nonzero=graph_source_nonzero,
        viscosity_leading_coefficient_nonzero=viscosity_leading_coefficient_nonzero,
        exact_profile_cancellation_excluded=viscosity_leading_coefficient_nonzero,
        quantitative_epsilon_gain_proved=False)


def rate_gate(kappa=Q(1,100000)):
    o=exponents(kappa)
    return dict(
        orders=o,
        leading_source='M L Lambda_re',
        leading_fast='-c_i^-1 M Nt^-1 L Lambda_re',
        leading_radial_viscosity='epsilon c_i^-1 M^3 L^2 Nt^-1 L Lambda_re',
        current_covariance_order=o['alpha'],
        current_source_order=o['D'],
        current_viscosity_order=o['nu'],
        strict_sufficient_inequality=(
            'there exist delta>0 and finite positive P_n such that '
            'N_n(L^2 Nt^-1 L Lambda_re)<=epsilon^(alpha+delta) P_n(S) '
            'for every selected component and accepted band'),
        target_alternative=(
            'Xi_re is identically zero, or has the same positive epsilon-order gain, '
            'on every active selected component'),
        strict_gain_from_current_evidence=False)
