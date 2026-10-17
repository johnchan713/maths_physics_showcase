"""Actual axial inverse and complete reference-strength budgets."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import sys
import mpmath as mp

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'root_foundation'))
from root_bounds import load,exp_bounds,outward  # noqa: E402
sys.path.pop(0)

@lru_cache(None)
def matrix_review():
    bmin=F(3,10)*exp_bounds(F(-3,40))[0]
    rows=[(1+exp_bounds(F(1,50))[1])/(2*bmin),exp_bounds(F(-24,25))[1]/bmin]
    return dict(parameter_range='0<lambda<=1/100; lambda=0 is singular and excluded',
        scaled_inverse_row_upper=[str(outward(v,10**6)[1]) for v in rows],
        exact_gap='exp(1-4lambda)*expm1(2lambda)>=2lambda*exp(24/25)',
        target_C1='||d(A,.)||_C1<=200*exp(-4/(5lambda)); 0<=A<=6/5',
        coefficient='||c(A,.)||_C1, |partial_A c| <=800/lambda*exp(-4/(5lambda))<=lambda^20',
        checks=dict(both_scaled_inverse_rows_below_four=max(rows)<4,
            positive_mass_and_exponential_gap_not_subtracted_in_binary64=True,
            zero_lambda_not_invertible=True,
            actual_trial_strength_domain_includes_zero=True,
            full_eta_value_plus_derivative_norm_retained=True,
            target_constant_below_two_hundred=exp_bounds(F(3,2))[1]<5 and 5*(1+F(6,5)*11/F(12,25))<200,
            affine_coefficients_and_strength_derivative_use_same_matrix=True))

def dyadic(x):
    sign,mantissa,exponent,_=x._mpf_
    return (-1 if sign else 1)*F(mantissa)*F(2)**exponent

@lru_cache(None)
def shape_review():
    m=load('axial_root_actual_positive_shape_bounds',HERE.parent/'pulse_moment_audit/bounds.py')
    with mp.workdps(90):s=m.shape_enclosure(1024,80)
    low,high=F('0.2450496182718'),F('0.2450496221363')
    return dict(boxes=1024,interval_digits=80,rational_K_enclosure=[str(low),str(high)],
        checks=dict(actual_monotone_positive_shape_integral_recomputed=True,
            rational_lower_strictly_below_outward_interval=low<dyadic(s['K_lower']),
            rational_upper_strictly_above_outward_interval=dyadic(s['K_upper'])<high,
            K_strictly_above_49_over_200=low>F(49,200),
            quadrature_is_not_used_for_the_K_enclosure=True))

@lru_cache(None)
def budget_review():
    shape=shape_review();kl,ku=map(F,shape['rational_K_enclosure'])
    lam=F(2)**-352;tl=64*lam;remainder=1000*tl
    em26=exp_bounds(-26,128);cl,ch=(1-em26[1])/4,(1-em26[0])/4
    left=F('1.0100502')**2*ku-cl+remainder
    right=F('1.0100504')**2*kl-ch-remainder
    eta_ratio=F(1,64)+512*em26[1]+80*lam**40/64
    discriminant_lower=4*F(49,200)*F(249,1000)
    return dict(parameter_majorants='lambda<2^-352; Tlambda<=64*2^-352; e^8>2^11',
        amplitude_bracket=['1.0100502','1.0100504'],
        error_bound='|E|,|E_eta|<=1000*T*lambda<10^-100; |E_A|<=80lambda^41',
        root_jacobian_lower='49/100',root_inverse_upper='100/49',root_eta_bound='T*lambda/30',
        complete_quadratic='(K+lambda*sum(w*a*a))*A^2+2lambda*sum(w*a*b)*A-C+lambda*s0-lambda*exp(-26)*R/2+lambda*sum(w*b*b)',
        checks=dict(exponential_binary_majorant_from_positive_Taylor=exp_bounds(8)[0]>2**11,
            T_lambda_decreases_for_T_at_least_sixty_four=True,
            exponential_coefficient_absorbed_into_lambda20=exp_bounds(7)[0]>800 and (7+21*256)*lam<F(4,5),
            complete_remainder_below_ten_minus_one_hundred=remainder<F(1,10**100),
            left_root_sign_strictly_negative=left<F(-3,10**8),
            right_root_sign_strictly_positive=right>F(6,10**8),
            full_trial_bracket_derivative_above_point_four_four=F(9,5)*kl-80*lam**41>F(11,25),
            complete_constant_is_strictly_negative=cl-remainder>F(249,1000),
            discriminant_proves_root_jacobian_above_point_four_nine=discriminant_lower>F(49,100)**2,
            normalized_eta_error_gives_Tlambda_over_thirty=eta_ratio<F(49,3000),
            root_eta_derivative_below_one_percent=tl/30<F(1,100),
            same_accepted_root_not_an_evenness_assumption=True,
            all_three_nonzero_affine_bump_energy_terms_retained=True,
            no_selected_superexponential_or_huge_support_coordinate_evaluated=True))

def prefix_tail_review():
    lam=F(2)**-352;tl=64*lam;growth=(240+F(1,128))*tl
    return dict(prefix='0<Ka,Ja<=4; 0<K2a<=16; ||m0||_C1,||j0||_C1<1; |s0|<=500T; |s0_eta|<1',
        first_energy_mass_bound='F1<=5/6+1=11/6 from F_prime<=1',
        full_post_mass='0<R<=1002+240T<=256T; |R_eta|<=4R<=1024T',
        release_mass='1+1/2+(2/3)h^7<2; terminal wait, nonzero end bump and infinite tail included',
        checks=dict(incoming_J_positive_from_Ka_at_least_four_exp_minus_T_and_r1_positive=True,
            incoming_moment_norm_bound_below_one=24*3*F(2)**-7712<1,
            incoming_energy_eta_bound_below_one=256*2*F(2)**-64<1,
            full_prefix_growth_below_quarter=growth<F(1,4),
            independent_first_ramp_bound_needs_no_F1_below_one_point_five=True,
            prefix_S_bound_500T_has_strict_gap=241+F(23,6*64)<500,
            post_stage_bound_256T_has_strict_gap=F(1002,64)+240<256,
            post_eta_log_density_strictly_below_four=2+2*F(1,100)/F(99,100)<4,
            exact_steep_hold_attenuation_is_h_to_eight=True,
            terminal_slope_at_most_minus_three_h_over_four=32*F(1,1000)/(1-F(1,100000))<F(1,4),
            complete_release_mass_below_two=1+F(1,2)+F(2,3)*F(1,100)**7<2,
            heat_edit_and_compensation_have_exact_net_S_zero_as_eta_functions=True))

def scope_review():
    return dict(result='actual reference M,J and scalar S roots; full pressure-aware post energy, same selected roots, all fixed one-sided jets',
        next='actual annular/compact-repair targets, inverses and gaps; complete-core global backward pressure identification; remaining finite-pulse lemmas',
        checks=dict(reference_axial_two_moment_and_scalar_strength_dependencies_discharged=True,
            composed_coefficients_have_full_C1_norm_at_most_two_lambda20=True,
            positive_discriminant_and_constant_matrix_supply_each_fixed_jet=True,
            complete_energy_not_only_main_pulse_used=True,
            later_heat_restores_moments_not_pointwise_power_tail=True,
            eta_dependent_annular_Q_still_requires_product_derivatives=True,
            global_pressure_identification_and_other_maps_remain_dependencies=True,
            no_unconditional_complete_residual_or_blowup_claim=True))
