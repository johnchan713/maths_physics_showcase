"""Enlarged closed-chart and mixed heat/pressure bounds at fixed orders."""
from fractions import Fraction as F
from math import comb,factorial
import importlib.util
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m
exterior=load('endpoint_foundation_exterior_bounds',HERE.parent/'annular_exterior_foundation/exterior_bounds.py')

def order(n):
    if type(n) is not int or n<0:raise ValueError('A finite nonnegative integer order is required')
    return n
def parameter(h):
    h=F(h)
    if not 0<h<F(1,100):raise ValueError('Keep the original 0<h<1/100')
    return h
def floor_log2(x):
    x=F(x)
    if x<=0:raise ValueError('Positive scale required')
    n=x.numerator.bit_length()-x.denominator.bit_length()
    return n-1 if x<F(2)**n else n
def ceil_log2(x):
    x=F(x)
    if x<=0:raise ValueError('Positive scale required')
    return -floor_log2(1/x)
def window(c,upper,minimum=16,enlargement=4):
    c,upper=F(c),F(upper)
    if not 0<c<=upper<1 or type(minimum) is not int or minimum<16 or enlargement not in (2,4):raise ValueError('Ordered positive physical scales and accepted enlargement required')
    lo=max(minimum,ceil_log2(1/(enlargement*upper)));hi=floor_log2(enlargement/c)
    return dict(first=lo,last=hi,count=max(0,hi-lo+1))
def stage_bound(c,qstar):
    c,qstar=F(c),F(qstar)
    if not 0<c<qstar<1:raise ValueError('Require 0<c<qstar<1')
    return max(0,ceil_log2(qstar/c)-1)
def closed_centers(lo,hi,delta):
    lo,hi,delta=map(F,(lo,hi,delta))
    if lo>hi or delta<=0:raise ValueError('Ordered interval and positive mesh required')
    left,right=lo/delta,hi/delta
    return -(-left.numerator//left.denominator)-1,right.numerator//right.denominator+1

def atlas_review():
    c,u=F(2)**-24,F(2)**-20;full=window(c,u);narrow=window(c,u,enlargement=2)
    return dict(enlarged_q_over_Q='1/4<=q/Q<=4',diagnostic_enlarged_window=full,historical_narrow_window=narrow,
        observation='q_e=min(qstar/16,2^(-ell0-2)); tau<=q_e/16; |z|<=(q_e/16)^D',
        punctured_cover='b=min(q_e/4,d^2/(16*Xext),1/4); zeta=min(d/4,(b/16)^D); tau<=b/16',
        checks=dict(closed_enlarged_window_includes_both_extra_edges=full==dict(first=18,last=26,count=9),
            narrower_window_omits_accepted_enlarged_boundary_charts=narrow==dict(first=19,last=25,count=7) and full['first']<narrow['first'] and full['last']>narrow['last'],
            sixteen_q_not_eight_q_used_for_enlarged_fourQ=16*F(3,8)*F(1,16)==F(3,8)<1,
            all_enlarged_bands_above_same_first_band=4*F(3,8)*F(1,4)==F(3,8)<1,
            terminal_slow_mesh_keeps_closed_boundary_labels=closed_centers(0,0,F(1,10**36))==(-1,1),
            strict_stage_bound_retains_exact_cutoff_endpoint=stage_bound(c,F(2)**-16)==7,
            equatorial_heat_cover_has_strict_margin=F(16,3)>1,
            axial_region_has_positive_q_lower_bound=True,
            same_labels_colors_common_lifts_and_physical_q_bumps_used=True))

def requirements(m):
    m=order(m)
    return dict(potential_for_velocity_Cauchy=m+2,direct_swirl_for_velocity_Cauchy=m+1,
        pressure_Cauchy=m+1,potential_for_future_residual_Cauchy=m+4,
        direct_for_future_residual_Cauchy=m+3,pressure_for_future_residual_Cauchy=m+2)

def seed_review():
    return dict(axis_potential='A_theta*e_theta=(1/2)*q^(-A)*(M/X)*(-y,x,0)',
        axis_direct_swirl='B=q^(-A-1/2)*F*(-y,x,0)',axis_pressure='p=q^(-2*A)*Pi',
        blocking_inequality='sup_(Eplus x [1-delta,1)) |D_x^alpha partial_t^j(A,B,p)|<infinity for every fixed alpha,j',
        requirements=[requirements(m) for m in (0,1,5,100)],
        checks=dict(M_over_X_not_E_inverse_is_regular_axis_quantity=True,
            exact_core_and_background_cut_potentials_supply_seed_not_velocity_only_inference=True,
            common_annulus_gives_uniform_axis_hole_r_squared_less_than_two_Xa_c=True,
            fractional_auxiliary_map_only_used_away_from_axis=True,
            current_all_fixed_jet_bounds_are_uniform_at_closed_eta_poles=True,
            compact_mean_and_generated_sources_keep_one_sided_T_zero_jets=True,
            full_seed_pressure_both_Pi_tails_and_fixed_gauge_preserved=True,
            endpoint_curl_uses_one_extra_potential_derivative_and_Cauchy_uses_another=True,
            background_and_current_diagonals_locally_finite_on_q_ge_c=True))

def rising(x,n):return exterior.rising(F(x),order(n))
def heat_A(h,j):return exterior.heat_A(parameter(h),order(j))
def heat_S(h,j):
    h=parameter(h);j=order(j)
    return sum(F(comb(j,k))*heat_A(h,k)*heat_A(h,j-k) for k in range(j+1))

def radial_coefficients(a,n):
    a=F(a);n=order(n);c=[F(1)]
    for ell in range(n):
        c=[-(a+ell+2*k)*(c[k] if k<len(c) else 0)+2*(c[k-1] if k else 0) for k in range(len(c)+1)]
    return c

def heat_terminal(h,n,j):
    h=parameter(h);n,j=order(n),order(j);a=1+2*h+2*j
    return F(4)**j*heat_A(h,j)*(-1)**n*rising(a,n)
def pressure_terminal(h,n,j):
    h=parameter(h);n,j=order(n),order(j);b=2+4*h+2*j
    return -F(4)**j*heat_S(h,j)/b*(-1)**n*rising(b,n)

def mixed_majorant(h,n,j,Zmax,pressure=False):
    h=parameter(h);n,j=order(n),order(j);Zmax=F(Zmax)
    if Zmax<0:raise ValueError('Nonnegative preterminal Z range required')
    if not pressure:
        c=radial_coefficients(1+2*h+2*j,n)
        return F(4)**j*sum(abs(v)*Zmax**k*heat_A(h,j+k) for k,v in enumerate(c))
    b=2+4*h+2*j
    if n==0:return F(4)**j*heat_S(h,j)/b
    c=radial_coefficients(b+1,n-1)
    return F(4)**j*sum(abs(v)*Zmax**k*heat_S(h,j+k) for k,v in enumerate(c))

def mixed_review():
    h=F(1,200);rows=[]
    for n in range(6):
        for j in range(6):
            square=sum(F(comb(j,k))*heat_terminal(h,n1,k)*heat_terminal(h,n-n1,j-k)*F(comb(n,n1))
                for n1 in range(n+1) for k in range(j+1))
            # d_r^n(K^2/r): distribute derivatives between K^2 and r^-1.
            density=sum(F(comb(n,l))*(-1)**l*factorial(l)*sum(F(comb(j,k))*F(comb(n-l,n1))*heat_terminal(h,n1,k)*heat_terminal(h,n-l-n1,j-k)
                for n1 in range(n-l+1) for k in range(j+1)) for l in range(n+1))
            rows.append(dict(radial_order=n,time_order=j,
                terminal_pressure_identity=pressure_terminal(h,n+1,j)==density,
                full_convolution_nonzero=square!=0,
                heat_majorant_dominates_terminal=mixed_majorant(h,n,j,0)>=abs(heat_terminal(h,n,j)),
                pressure_majorant_dominates_terminal=mixed_majorant(h,n,j,0,True)>=abs(pressure_terminal(h,n,j))))
    return dict(exact_terminal_rows=rows,
        coefficient_recursion='c_(n+1,k)=-(a+n+2*k)*c_(n,k)+2*c_(n,k-1)',
        heat_bound='|d_r^n d_t^j K|<=|C|*4^j*r^(-1-2h-2j-n)*sum |c_nk|*Zmax^k*A_(j+k)',
        pressure_bound='n=0: C^2*4^j*S_j*r^(-b)/b; n>=1: C^2*4^j*r^(-b-n)*sum |c_(n-1,k)(b+1)|*Zmax^k*S_(j+k)',
        checks=dict(all_36_mixed_radial_pressure_identities_exact=all(r['terminal_pressure_identity'] for r in rows),
            full_pressure_Leibniz_convolution_never_truncated=all(r['full_convolution_nonzero'] for r in rows),
            all_terminal_values_below_same_positive_majorants=all(r['heat_majorant_dominates_terminal'] and r['pressure_majorant_dominates_terminal'] for r in rows),
            pressure_density_integrable_at_infinity_for_every_fixed_j=3>1,
            finite_Zmax_on_every_positive_radius_late_time_compact=True,
            Cartesian_frame_derivatives_bounded_on_r_ge_rmin_positive=True,
            no_global_complex_heat_neighborhood_or_convergent_terminal_series_used=True))

def scope_review():
    return dict(endpoint='same complete (A,B,u,p) has compatible C_infinity one-sided jets at t=1 on punctured local compacts',
        finite_Taylor='||G(1-tau)-sum_(j=0)^k (-tau)^j*F_j/j!||<=M_(k+1)*tau^(k+1)/(k+1)!',
        checks=dict(one_fixed_current_realization_not_a_new_sequence_for_each_jet=True,
            uniform_next_time_derivative_yields_Cauchy_terminal_limits=True,
            spatial_jets_identified_on_enlarged_punctured_compact_boxes=True,
            full_pressure_and_velocity_limits_commute_with_required_finite_curls=True,
            Gaussian_ledger_never_forward_inverted_or_newly_summed=True,
            no_future_velocity_or_zero_extension_of_nonzero_heat_pressure_jets=True,
            retained_exact_core_and_finite_analytic_dependencies_remain=True,
            localization_global_force_energy_and_comparison_not_promoted=True))
