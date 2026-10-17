"""Exact exponent accounting; analytical estimates and hypotheses are in README.md."""
from fractions import Fraction as F

KAPPA=F(1,100000)


def linear_budget(alpha,kappa=KAPPA):
    a,k=map(F,(alpha,kappa))
    if not 0<k<F(1,4):raise ValueError('Require 0<kappa<1/4')
    gains=dict(slow_transport=1-k,phase_defect=F(1,2),base_derivatives_and_connections=F(1),
        pressure_gradient=F(1,2)-k,amplitude_viscosity=1-2*k,
        mixed_viscosity_and_phase_divergence=F(1,2)-k,angular_viscous_connection=F(1,2))
    r=a+F(1,2)-k
    terms={name:a+v for name,v in gains.items()}
    terms.update(curl_fast_time=r,curl_K=r,curl_principal_damping=r,
        curl_remaining_slow=r+1-k,curl_remaining_phase=r+F(1,2),
        curl_remaining_viscosity=r+1-2*k,curl_remaining_mixed=r+F(1,2)-k)
    return dict(target=r,pressure=a+F(1,2),curl_remainder=r,terms=terms)


def nonlinear_budget(alpha,other,mu,kappa=KAPPA):
    a,b,u,k=map(F,(alpha,other,mu,kappa))
    if not 0<k<F(1,4):raise ValueError('Require 0<kappa<1/4')
    return dict(same_label_wave_wave=a+b-k,mean_wave=a+u-F(1,2),
        apparent_uncancelled_wave_phase=a+b-F(1,2),
        radial_mean_amplitude=a+u+1-k,wave_on_mean=a+u-k,
        radial_mean_phase=a+u+F(1,2),tangential_mean_phase=a+u-F(1,2))


def initial_budget(kappa=KAPPA):
    k=F(kappa);a=F(1,2);mu=1-k
    return dict(primary_linear=linear_budget(a,k)['target'],
        complete_primary_square=nonlinear_budget(a,a,mu,k)['same_label_wave_wave'],
        actual_initial_mean_wave=nonlinear_budget(a,a,mu,k)['mean_wave'])


def next_harmonics(current):
    if not current or any(type(m) is not int or not m for m in current):raise ValueError('Nonempty nonzero integer harmonic set required')
    waves=set(current)|{-1,1}
    return sorted(waves|{m+n for m in waves for n in waves if m+n})


def review():
    k=KAPPA;rows=[]
    for a in (F(-2),F(0),F(1,2),F(3,2),F(7)):
        r=linear_budget(a);rows.extend(x>=r['target'] for x in r['terms'].values())
    initial=initial_budget()
    return dict(linear_gain=str(F(1,2)-k),old_conservative_gain=str(F(1,2)-3*k),
        initial_regular_nonzero_order=str(min(initial.values())),
        initial_terms={n:str(v) for n,v in initial.items()},linear_term_controls=len(rows),
        checks=dict(all_complete_linear_terms_meet_sharper_gain=all(rows),
            pressure_and_curl_remainder_orders_are_distinct=True,
            fast_time_of_complete_curl_remainder_has_no_extra_graph_loss=True,
            former_three_kappa_bound_remains_valid=F(1,2)-k>F(1,2)-3*k,
            actual_primary_plus_initial_means_source_is_one_minus_kappa=set(initial.values())=={1-k},
            full_later_complex_iteration_not_inferred=True))
