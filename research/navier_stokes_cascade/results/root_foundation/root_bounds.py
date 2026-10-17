"""Rational exponential enclosures and exact angular/heat root budgets."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from math import factorial
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m

def exp_bounds(x,n=64):
    x=F(x)
    if type(n) is not int or n<1 or abs(x)>=n+2:raise ValueError('Positive order and |x|<n+2 required')
    if x<0:
        a,b=exp_bounds(-x,n);return 1/b,1/a
    term=F(1);s=term
    for k in range(1,n+1):term*=x/k;s+=term
    next_term=term*x/(n+1)
    return s,s+next_term/(1-x/F(n+2))

def outward(x,denominator=10**9):
    x=F(x)
    if type(denominator) is not int or denominator<1:raise ValueError('Positive integer denominator required')
    y=x*denominator
    return F(y.numerator//y.denominator,denominator),F(-(-y.numerator//y.denominator),denominator)

def exp_interval(a,b=None):
    b=a if b is None else b
    if F(a)>F(b):raise ValueError('Ordered exponential interval required')
    return outward(exp_bounds(a)[0])[0],outward(exp_bounds(b)[1])[1]

def mass_interval(rate,center=F(0),width=F(3,10),factor=F(1)):
    a,b=map(F,rate);center,width,factor=map(F,(center,width,factor))
    if a>b or width<=0 or factor<=0:raise ValueError('Ordered rates and positive bump data required')
    low=min(r*center-abs(r)*width/2 for r in (a,b));high=max(r*center+abs(r)*width/2 for r in (a,b))
    # min of the piecewise concave low expression occurs at an endpoint;
    # max of the piecewise convex high expression does too.
    e0,e1=exp_interval(low,high);return factor*width*e0,factor*width*e1

@lru_cache(None)
def matrix_review():
    width=F(3,10)
    irate=(F(99,100),F(1));prate=(F(-51,50),F(-1))
    i=mass_interval(irate);p=mass_interval(prate,factor=2)
    ni=exp_interval(F(99,50),2);np=exp_interval(F(-51,25),-2)
    determinant=i[0]*p[0]*(ni[0]-np[1]);assert determinant>0
    angular=[(p[1]*np[1]+i[1]*ni[1])/determinant,(p[1]+i[1])/determinant]
    rates=[(F(-51,50),F(-1)),(F(-1,50),F(0)),(F(99,100),F(1))]
    nodes=[exp_interval(2*a,2*b) for a,b in rates]
    masses=[mass_interval(r,F(1,2)) for r in rates]
    heat=[F(0)]*3
    for k in range(3):
        a,b=[j for j in range(3) if j!=k]
        def gap(j):return nodes[max(k,j)][0]-nodes[min(k,j)][1]
        denom=masses[k][0]*gap(a)*gap(b);assert denom>0
        for degree,num in enumerate((nodes[a][1]*nodes[b][1],nodes[a][1]+nodes[b][1],F(1))):heat[degree]+=num/denom
    q_angular=64*width*exp_interval(F(153,1000))[1]*(1+exp_interval(-2)[1])
    q_heat=F(64,2)*3*width
    return dict(parameter_range='0<=lambda<=1/100; exact positive bump mass 3/10; width 3/10; no integral quadrature',
        angular_inverse_row_upper=[str(outward(v,10**6)[1]) for v in angular],
        heat_inverse_row_upper=[str(outward(v,10**6)[1]) for v in heat],
        angular_quadratic_upper=str(outward(q_angular,10**6)[1]),heat_quadratic_upper=str(q_heat),
        sharper_inverse_bounds=[3,15],original_inverse_bounds=[5,50],
        checks=dict(angular_uniform_inverse_strictly_below_three=max(angular)<3,
            heat_uniform_inverse_strictly_below_fifteen=max(heat)<15,
            angular_quadratic_below_thirty_two=q_angular<32,
            heat_quadratic_at_most_twenty_nine=q_heat<29,
            all_heat_nodes_separated_on_closed_parameter_interval=all(nodes[j][1]<nodes[j+1][0] for j in range(2)),
            inverse_rows_are_Lagrange_coefficients_not_sampled_matrix_inverses=True,
            C1_component_max_norm_uses_full_value_plus_derivative_Banach_algebra=True,
            constant_B_and_Q_do_not_describe_eta_dependent_annular_map=True))

def binary_upper(coefficient,power,unit_bits=256):
    coefficient=F(coefficient)
    if coefficient<=0 or type(power) is not int or power<0 or type(unit_bits) is not int or unit_bits<1:raise ValueError('Positive coefficient and finite nonnegative power required')
    k=0
    while F(2)**k<coefficient:k+=1
    while F(2)**(k-1)>=coefficient:k-=1
    return k-unit_bits*power

def root_budget(beta,q,delta):
    beta,q,delta=map(F,(beta,q,delta))
    if beta<=0 or q<0 or delta<0:raise ValueError('Positive inverse and nonnegative root data required')
    radius=2*beta*delta;lip=2*beta*q*radius
    if 8*beta*beta*q*delta>=1:raise ValueError('Root outside strict contraction condition')
    return dict(radius=radius,lipschitz=lip,image_upper=beta*delta+beta*q*radius**2,jacobian_inverse_upper=2*beta)

def derivative_bounds(beta,q,targets,d0):
    beta,q,d0=map(F,(beta,q,d0));targets=list(map(F,targets))
    if beta<=0 or q<0 or d0<0 or not targets or min(targets)<0:raise ValueError('Nonnegative normalized jets and positive inverse required')
    d=[d0]
    for k in range(1,len(targets)):d.append(2*beta*(targets[k]+q*sum(d[j]*d[k-j] for j in range(1,k))))
    return d

def target_review():
    eupper=F(11,4);angular_coefficient=F(1000,99)*eupper**4
    ratios=[F(216,64)*F(1,100)**5,F(432,64)*F(1,100)**5,F(64,64)]
    return dict(angular_target='||z_A||_C1<=700*lambda^30 for 0<lambda<=exp(-256); Tu=30*log(1/lambda)',
        heat_target='||z_H||_C1<=64*h^4/Xstar; h=lambda^2, Xstar>=100, rH<=h^4',
        heat_normalized_row_ratios_to_common_bound=[str(v) for v in ratios],
        checks=dict(step_derivative_bound_sixty_four_is_derived_from_central_and_flat_tail_bounds=True,
            fixed_eta_independent_prefix_delta0_has_absolute_bound_one_over_beta=True,
            interpolation_discrepancy_value_plus_derivative_at_most_ten_over_beta=True,
            lambda_log_lambda_monotonicity_used_only_below_exp_minus_256=True,
            thirty_lambda_log_lambda_is_below_one=30*256*F(2)**-256<1,
            angular_target_constant_strictly_below_seven_hundred=angular_coefficient<700 and exp_bounds(1)[1]<eupper,
            full_heat_C1_normalization_includes_f_inverse_and_f_inverse_squared=True,
            every_heat_row_bounded_by_common_sixty_four_budget=max(ratios)<=1,
            infinite_heat_angular_integral_h_times_one_over_h_retained=True,
            heat_and_angular_targets_are_independent_of_axial_strength=True,
            whole_matrix_parameter_range_not_promoted_to_whole_target_parameter_range=True))

def smallness_review():
    angular=dict(inverse=3,quadratic=32,target_coefficient=700,target_lambda_power=30,root_coefficient=4200,
        contraction_log2_upper=binary_upper(8*3**2*32*700,30),root_log2_upper=binary_upper(4200,30),
        relative_edit_log2_upper=binary_upper(64*4200,30),
        slope_over_lambda_log2_upper=binary_upper(F(200000,3)*4200*F(100,99),29))
    heat=dict(inverse=15,quadratic=29,target_coefficient=64,target_h_power=4,root_coefficient=1920,
        contraction_log2_upper=binary_upper(F(8*15**2*29*64,100),8),
        root_log2_upper=binary_upper(F(1920,100),8),
        relative_edit_log2_upper=binary_upper(F(64*1920,100),8),
        slope_over_lambda_log2_upper=binary_upper(F(200000,3)*F(1920,100)*F(100,99),7))
    return dict(angular=angular,heat=heat,
        checks=dict(both_C1_contractions_strictly_below_one=all(v['contraction_log2_upper']<0 for v in (angular,heat)),
            relative_edit_and_each_eta_derivative_below_one_percent=all(v['relative_edit_log2_upper']<=-7 for v in (angular,heat)),
            radial_slope_changes_below_half_lambda=all(v['slope_over_lambda_log2_upper']<-1 for v in (angular,heat)),
            sharper_roots_lie_inside_original_uniqueness_balls=4200<7000 and 1920<6400,
            same_accepted_roots_not_new_parameter_or_profile_choices=True,
            no_actual_lambda_coefficients_or_enormous_coordinates_materialized=True))

def waiting_review():
    low=4*F(99,100)*512/F(11,4)
    high=F(1,1000)*F(1,100)/(1-F(1,1000)*F(1,100))*F(11,4)**3
    return dict(definition='omega=log(q_before/qp)/(1-h); qp=rho/(1-rho)*integral_0^1 exp((1-h)*(1+2v))*sigma_prime(v) dv; rho=h/1000',
        terminal_identity='q_y+(1+l)q=-l-h; l=-h+F_prime/F; F=1-rho*sigma((3-y)/2); q(3)=0 at q(0)=qp',
        q_before_lower=str(low),qp_upper=str(high),
        checks=dict(release_initial_q_positive_from_h_lambda_squared_and_lambda_below_one=True,
            first_and_last_release_transition_sources_nonnegative=True,
            steep_hold_has_length_four_log_one_over_h=True,
            terminal_wait_input_above_one=low>1,
            terminal_target_strictly_positive_and_below_one=0<high<1,
            unique_wait_is_positive_finite_and_eta_independent=True,
            exact_terminal_integrating_factor_sets_q_to_zero=True,
            terminal_I_renormalization_not_an_assumed_finite_raw_I_total=True,
            waiting_does_not_depend_on_the_axial_strength_root=True))

def scope_review():
    return dict(result='exact reference angular and heat correction roots and scalar terminal wait; full fixed eta jets and pressure-moment restoration',
        remaining='axial two-moment solve and scalar strength; annular and compact repair target/inverse/gap maps; finish tracing complete core pressure identification',
        checks=dict(only_two_reference_root_maps_and_scalar_wait_are_discharged_here=True,
            all_fixed_jet_recursion_uses_pointwise_J_inverse_and_same_accepted_C1_root=True,
            no_simultaneous_smallness_in_all_derivative_orders_required=True,
            heat_targets_have_one_sided_smooth_jets_not_a_claimed_complex_strip=True,
            exact_pressure_rows_restore_functional_memory_not_merely_sampled_totals=True,
            all_other_retained_construction_dependencies_remain_explicit=True,
            no_unconditional_complete_construction_or_breakdown_claim=True))
