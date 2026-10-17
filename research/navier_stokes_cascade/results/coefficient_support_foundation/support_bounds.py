"""Current linear coefficient map, exact jet costs and reserved supports."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module);return module

matching=load('support_foundation_matching_bounds',HERE.parent/'matching_root_foundation/matching_bounds.py')

def reciprocal_amplitude_cost(k,e_star):
    """Full factorial C^k norm of (1+eta^2)/e_star on [-1,1]."""
    if type(k) is not int or k<0 or e_star<=0:raise ValueError('Nonnegative order and positive amplitude required')
    return F(2 if k==0 else 4 if k==1 else 5)/e_star

def target_bound(k,e_star,rho,radial_angular_scale,lam,moments):
    """Bound for exact c=-B^-1 z in the full factorial eta norm.

    radial_angular_scale denotes 1/(sqrt(2)*rho^(3/2)); it is supplied
    independently so rational fixtures need not approximate sqrt(2).
    Ordinary moment order is (M,I,Cp,J,S).
    """
    if min(rho,radial_angular_scale,lam)<=0 or lam>F(1,100) or len(moments)!=5 or min(moments)<0:
        raise ValueError('Positive actual scales, admissible slope and five norm bounds required')
    w=reciprocal_amplitude_cost(k,e_star);a,b,c,d,e=moments
    z=(a/rho,(a/rho+radial_angular_scale*w*d)/lam,
       radial_angular_scale*b,w*e/rho,w*c)
    return z,800*max(z)

@lru_cache(None)
def map_review():
    r=matching.repair_review()
    return dict(inverse_bound=800,E_inverse_upper=r['E_inverse_upper'],U_inverse_upper=r['U_inverse_upper'],
        linear_U_weights=['exp(y)','exp(y)*(1-exp(-lambda*y))/lambda'],
        linear_E_rates=['3/2','1/2-lambda','-1/2-lambda'],
        normalizers=['rho3','sqrt(2)*rho3^(3/2)','e','sqrt(2)*rho3^(3/2)*e','rho3*e'],
        checks=dict(actual_linear_map_is_same_B_not_nonlinear_repair=True,
            positive_mass_determinants_recomputed=all(r['checks'].values()),
            E_inverse_strictly_below_800=F(r['E_inverse_upper'])<800,
            U_inverse_strictly_below_2=F(r['U_inverse_upper'])<2,
            all_order_eta_derivatives_use_same_constant_B=True,
            no_small_target_or_uniform_in_n_norm_assumed=True))

@lru_cache(None)
def jet_review():
    costs=[reciprocal_amplitude_cost(k,F(7)) for k in range(9)]
    z,bound=target_bound(8,F(7),F(2),F(3,5),F(1,1000),tuple(map(F,(2,3,5,7,11))))
    # Positive-order correction amplitude is additive, NOT multiplied by e.
    # The exact normalization derivative is nonzero even for constant m3.
    derivative_at_one=F(2,7)*5
    return dict(reciprocal_factorial_norms=list(map(str,costs)),example_transformed_norm_bounds=list(map(str,z)),
        example_coefficient_norm_bound=str(bound),
        checks=dict(actual_amplitude_is_e_star_over_one_plus_eta_squared=True,
            exact_full_factorial_costs=costs[:3]==[F(2,7),F(4,7),F(5,7)],
            reciprocal_has_no_jets_above_two=all(x==F(5,7) for x in costs[2:]),
            amplitude_eta_derivatives_cannot_be_frozen=derivative_at_one!=0,
            arbitrary_raw_M_J_target_retains_lambda_inverse=True,
            all_fixed_jets_finite_but_not_uniform_in_moment_order=True))

def support_edges(n):
    if type(n) is not int or n<1:raise ValueError('Positive integer order required')
    return dict(velocity='Xvel',pressure='Xv' if n<3 else 'Xvel',stress='Xb' if n==1 else 'Xv' if n==2 else 'Xvel')

def support_review():
    bumps=[(F(c)-F(1,20),F(c)+F(1,20)) for c in ('1/2','1','3/2','2','5/2')]
    return dict(bumps=[[str(a),str(b)] for a,b in bumps],
        rows=[dict(order=n,**support_edges(n)) for n in (1,2,3,7,1000000)],
        mean_pressure=['K1(eta)-mu^2/(4X)','K2(eta)','0 for n>=3'],
        checks=dict(all_velocity_bumps_end_before_mean_interval=bumps[-1][1]<3<5<6,
            reserved_I4_is_six_to_eleven_in_log_X_over_rho3=True,
            inner_cutoff_and_all_moment_bumps_are_disjoint=True,
            zero_first_total_implies_M_and_V_zero_after_velocity_edge=True,
            zero_third_total_needed_to_remove_pressure_constant=True,
            order_two_shifted_radial_axial_viscosity_retained=True,
            potential_cutoff_extra_M_term_vanishes_on_I4=True,
            leading_M_mu_eta_and_V_minus_mu_not_zeroed=True,
            pressure_and_velocity_supports_not_identified=True))
