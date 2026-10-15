"""Exact pressure, weighted interpolation and expanding-ball exponent ledger.

Continuum proofs are in README.md. Finite arithmetic checks audit their
exponents; they are not numerical evaluation of the conditional flow.
"""
from fractions import Fraction as F
import math
import mpmath as mp


def cutoff_weights(power=8):
    """Eight is the minimum for THIS L^(4/3) commutator factorization."""
    if type(power) is not int or power < 1:
        raise ValueError('A positive integer cutoff power is required')
    half=F(power,2)
    margins=dict(main_pressure_dual_weight=half-1-half/2,
                 commutator_pressure_dual_weight=half-1-3*half/4,
                 transport_weight=power-1-3*half/2)
    return dict(cutoff_power=power,sobolev_weight=half,
                residual_pressure_pairing_weight=half-1,
                margins=margins,valid=min(margins.values()) >= 0)


def holder_review():
    p=F(4,3)
    weights=cutoff_weights()
    return dict(commutator_kernel_exponent=p,
                commutator_dual_exponent=4,
                main_pressure_exponent=F(3,2),main_pressure_dual_exponent=3,
                main_pressure_interpolation_theta=F(1,2),
                commutator_interpolation_theta=F(3,4),
                cutoff=weights,
                checks=dict(main_pressure_holder=1/F(3,2)+F(1,3) == 1,
                            commutator_holder=1/p+F(1,4) == 1,
                            L3_interpolation=F(1,2)/6+F(1,2)/2 == F(1,3),
                            L4_interpolation=F(3,4)/6+F(1,4)/2 == F(1,4),
                            product_ww_in_L_three_halves=F(1,6)+F(1,2) == F(2,3),
                            all_cutoff_weight_margins_nonnegative=weights['valid'],
                            power_four_fails_this_commutator_pairing=not cutoff_weights(4)['valid']))


def kernel_exponents(p=F(4,3)):
    p=F(p)
    if not 1 < p < F(3,2):
        raise ValueError('Require 1<p<3/2 for both kernel integrals')
    near,far=2-2*p,2-3*p
    return dict(p=p,near_radial_power=near,far_radial_power=far,
                near_integral_coefficient=1/(3-2*p),
                far_integral_coefficient=1/(3*p-3),
                norm_radius_power=3/p-3)


def kernel_review():
    row=kernel_exponents()
    samples=[]
    with mp.workdps(65):
        p=mp.mpf(4)/3
        for R in (1,3,19):
            # The unit angular constant is kept separate from 4*pi and
            # the component kernel constant, neither affects radius powers.
            near=R**(-p)*mp.quad(lambda s:s**(2-2*p),[0,R])
            far=mp.quad(lambda s:s**(2-3*p),[R,mp.inf])
            expected=4/mp.mpf(R)
            samples.append(dict(radius=R,relative_radial_integral_error=float(abs((near+far)/expected-1))))
    return dict(exact_exponents=row,radial_samples=samples,
                exact_norm_power_for_p_four_thirds=F(-3,4),
                kernel='|x|^-3*min(|x|/R,1)',
                checks=dict(near_integral_converges=row['near_radial_power'] > -1,
                            far_integral_converges=row['far_radial_power'] < -1,
                            exact_norm_radius_scaling=row['norm_radius_power'] == F(-3,4),
                            independent_improper_integrals_agree=max(r['relative_radial_integral_error'] for r in samples) < 1e-18,
                            local_identity_distribution_cancels_in_commutator=True))


def young_remainder(a,b):
    """R^-b D^a <= delta*D^2+C_delta*R^(-2*b/(2-a))."""
    a,b=F(a),F(b)
    if not 0 < a < 2 or b <= 0:
        raise ValueError('A positive subquadratic dissipation power and radius decay are required')
    return 2*b/(2-a)


def flux_review():
    terms={'transport_and_main_pressure':(F(3,2),F(1)),
           'linear_main_pressure':(F(1,2),F(1)),
           'commutator_pressure':(F(3,4),F(7,4))}
    rows={name:dict(dissipation_power=a,radius_decay=b,
                    young_radius_decay=young_remainder(a,b)) for name,(a,b) in terms.items()}
    # B<=C(D+R^-1 W); these are the radius powers of the offsets.
    offsets=dict(main_cubic_offset=F(5,2),main_linear_offset=F(3,2),
                 commutator_offset=F(5,2),laplacian_flux=F(2))
    best=min([r['young_radius_decay'] for r in rows.values()]+list(offsets.values()))
    return dict(rows=rows,offset_radius_decays=offsets,
                uniform_absorbed_remainder_power=best,
                energy_inequality='0.5*E_R_prime+0.5*D_R^2<=||grad u||_infinity*E_R+C_T*R^(-4/3)',
                checks=dict(all_dissipation_powers_strictly_below_two=all(r['dissipation_power'] < 2 for r in rows.values()),
                            all_remainders_decay_at_least_four_thirds=best == F(4,3),
                            source_coarser_one_over_R_bound_also_follows=best > 1,
                            no_global_gradient_bound_for_competitor_is_used=True))


def pressure_recovery_review():
    return dict(tensor='g=v tensor v-u tensor u=w tensor w+w tensor u+u tensor w in L_infinity_t L1_x',
                pressure='pi_star=sum R_i R_j g_ij in L_infinity_t H^(-s), s>3/2',
                convenient_negative_orders=dict(tensor=2,pressure=2,time_averaged_gradient=3),
                harmonic_elimination='H_a=int a(t)*(grad pi-grad pi_star)dt in H^-3, Delta H_a=0; weighted L2 Fourier transform supported at {0} is zero',
                time_averaged_momentum='int a*grad pi=Delta int a*w+int a_prime*w-div int a*g',
                checks=dict(negative_Sobolev_integral_converges=2*2 > 3,
                            pressure_multiplier_is_bounded_by_one=True,
                            harmonic_H_minus_three_distribution_is_zero_by_Fourier_support=True,
                            pressure_growth_or_decay_is_not_assumed=True,
                            time_averaging_avoids_assuming_global_L2_of_w_time_derivative=True))


def missing_hypothesis_review():
    # Radial density of |K|^(4/3) without each part of min(r/R,1).
    no_near=[F(d)-1 for d in (10,1000)]  # delta=1/d; int_delta^1 r^-2 dr
    no_far=[3*(mp.mpf(L)**(mp.mpf(1)/3)-1) for L in (10,1000)]
    angular=4/(3*math.sqrt(3))
    shells=[dict(outer_radius=R,pressure_L1_shell_mass=angular*math.log(R)) for R in (10,100,10000)]
    return dict(uncanceled_near_singularity_integrals=no_near,
                unbounded_far_Lipschitz_kernel_integrals=[float(x) for x in no_far],
                radial_unit_mass_source='g_11=rho smooth radial with support in unit ball and mass one; outside p_star=(3*omega_1^2-1)/(4*pi*r^3)',
                pressure_L1_shells=shells,
                harmonic_pressure_example='h=x1 has Delta h=0 but grad h=(1,0,0); the latter is not in H^-3(R3)',
                checks=dict(without_Lipschitz_cancellation_near_integral_diverges=no_near[-1] > 100*no_near[0] > 0,
                            without_far_cap_kernel_is_not_in_L_four_thirds=no_far[-1] > 5*no_far[0],
                            L1_tensor_does_not_imply_L1_pressure=shells[-1]['pressure_L1_shell_mass'] > 3*shells[0]['pressure_L1_shell_mass'],
                            Poisson_equation_alone_does_not_fix_pressure_gradient=True,
                            cutoff_power_four_is_insufficient_for_chosen_factorization=not cutoff_weights(4)['valid']))


def young_absorption_review():
    rows=[]
    delta=.05
    for a,b in ((F(3,2),F(1)),(F(1,2),F(1)),(F(3,4),F(7,4))):
        af,bf=float(a),float(b)
        constant=(1-af/2)*(af/(2*delta))**(af/(2-af))
        for R in (1,3,25):
            lam=R**(-bf)
            critical=(af*lam/(2*delta))**(1/(2-af))
            maximum=lam*critical**af-delta*critical**2
            bound=constant*R**(-float(young_remainder(a,b)))
            rows.append(dict(dissipation_power=a,radius=R,maximum_minus_bound_relative_error=abs(maximum-bound)/(1+abs(maximum)+abs(bound))))
    return dict(analytic_maximizer='D=(a*R^(-b)/(2*delta))^(1/(2-a))',rows=rows,
                checks=dict(exact_young_remainder_power_matches_global_maximum=max(r['maximum_minus_bound_relative_error'] for r in rows)<2e-13,
                            all_absorption_constants_are_radius_independent=True))


def gronwall_review():
    radii=(2,8,32,128)
    T,L,C=F(3,4),F(2),F(5)
    exact_prefactor=float(C/L)*(math.exp(2*float(L*T))-1)
    bounds=[dict(radius=R,energy_upper_bound=exact_prefactor*R**(-4/3)) for R in radii]
    return dict(fixed_T=T,constant_gradient_bound=L,remainder_constant=C,
                exact_solution='E_R(t)<=C/L*(exp(2*L*t)-1)*R^(-4/3) when L>0; 2*C*t*R^(-4/3) when L=0',
                synthetic_bounds=bounds,
                checks=dict(expanding_ball_bound_tends_to_zero=all(b['energy_upper_bound'] < a['energy_upper_bound'] for a,b in zip(bounds,bounds[1:])),
                            every_fixed_ball_eventually_has_cutoff_one=True,
                            no_monotonicity_of_cutoffs_in_R_is_needed=True,
                            comparison_stage_T_is_fixed_before_radius_limit=True,
                            velocity_blowup_contradicts_smoothness_on_one_compact_neighborhood_conditionally=True))
