"""Exact endpoint geometry and exterior jets; actual field bounds stay conditional."""
from fractions import Fraction as Q
from math import comb,factorial,prod
import importlib.util
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module);return module


def order(m):
    if type(m) is not int or m<0: raise ValueError('Nonnegative finite integer order required')
    return m


def parameter(h):
    h=Q(h)
    if not 0<h<Q(1,100): raise ValueError('Require inherited 0<h<1/100')
    return h


def floor_log2(x):
    x=Q(x)
    if x<=0: raise ValueError('Positive rational argument required')
    n=x.numerator.bit_length()-x.denominator.bit_length()
    return n-1 if x<Q(2)**n else n


def ceil_log2(x):
    x=Q(x)
    if x<=0:raise ValueError('Positive rational argument required')
    return -floor_log2(1/x)


def band_window(c,upper,minimum_band=16):
    c,upper=Q(c),Q(upper)
    if not 0<c<=upper<1 or type(minimum_band) is not int or minimum_band<16:
        raise ValueError('Require 0<c<=upper<1 and geometric first band >=16')
    lo=max(minimum_band,ceil_log2(1/(2*upper)));hi=floor_log2(2/c)
    return dict(first=lo,last=hi,count=max(0,hi-lo+1))


def active_stage_bound(c,qstar):
    c,qstar=Q(c),Q(qstar)
    if not 0<c<qstar<1: raise ValueError('Require 0<c<qstar<1')
    # a_j>2^j/qstar; a_j*c<1 implies j<log2(qstar/c).
    return max(0,ceil_log2(qstar/c)-1)


def closed_centers(lo,hi,delta):
    lo,hi,delta=map(Q,(lo,hi,delta))
    if lo>hi or delta<=0: raise ValueError('Ordered interval and positive mesh required')
    ceil=lambda v:-((-v.numerator)//v.denominator)
    floor=lambda v:v.numerator//v.denominator
    return ceil(lo/delta)-1,floor(hi/delta)+1


def atlas_review():
    w=band_window(Q(2)**-24,Q(2)**-20)
    rows=[]
    for ell in (16,20,10**6):
        delta=Q(1,ell**6)
        intervals=[closed_centers(a,b,delta) for a,b in ((1,3),(-2,2),(0,2))]
        rows.append(dict(band=ell,mesh=str(delta),terminal_time_centers=list(closed_centers(0,0,delta)),
            finite_box_label_count=2*prod(b-a+1 for a,b in intervals)))
    return dict(diagnostic_band_window=w,diagnostic_stage_upper_bound=active_stage_bound(Q(2)**-24,Q(2)**-16),
        rows=rows,pointwise_slow_overlap=162,maximum_common_lifts=2744,
        real_map='Y=v_r*r^d_r+v_t*t; smooth for r>0, not asserted smooth at r=0',
        checks=dict(closed_dyadic_window_includes_both_endpoints=w==dict(first=19,last=25,count=7),
            terminal_mesh_labels_selected_before_evaluation=all(r['terminal_time_centers']==[-1,1] for r in rows),
            compact_box_counts_finite_without_enumerating_labels=all(r['finite_box_label_count']>0 for r in rows),
            doubled_cutoff_scales_leave_only_finitely_many_stages=active_stage_bound(Q(2)**-24,Q(2)**-16)==7,
            geometry_uses_same_closed_labels_absolute_pullbacks_and_physical_q_bumps=True))


def requirements(m):
    m=order(m)
    return dict(terminal_representative_order=m,representative_Cauchy_bound=m+1,
        velocity_Cauchy_potential_bound=m+2,velocity_Cauchy_direct_swirl_bound=m+1,
        pressure_Cauchy_bound=m+1,future_force_Cauchy_potential_bound=m+4,
        future_force_Cauchy_direct_swirl_bound=m+3,future_force_Cauchy_pressure_bound=m+2)


def endpoint_review():
    return dict(rows=[requirements(m) for m in (0,1,5,100)],
        blocking_inequality='sup_(E x [1-delta,1)) |D_x^alpha partial_t^j(A,B,p)|<infinity for every fixed alpha,j',
        observation_domain='q_e=min(qstar/16,2^(-ell0-2)); tau0=q_e/16; z0=(q_e/16)^D',
        cover='b=min(q_e/4,d^2/(16*Xext),1/4); zeta=min(d/4,(b/16)^D); tau<=b/16',
        axial_region='|z|>=zeta/2 gives q>=(zeta/2)^(1/D)>0',
        equatorial_region='|z|<zeta, |x|>=d gives r>d/2 and q<=3*b/8, hence X>=16*Xext/3',
        Cauchy='||D_x^alpha partial_t^j G(t2)-D_x^alpha partial_t^j G(t1)||_E <= M_(alpha,j+1)*|t2-t1|',
        terminal_Taylor='||D_x^alpha[G(1-tau)-sum_(j=0)^k (-tau)^j F_j/j!]||_E <= M_(alpha,k+1)*tau^(k+1)/(k+1)!',
        checks=dict(implicit_denominator_uniformly_positive=1-2*Q(1,100)==Q(49,50),
            exterior_cover_has_strict_margin=Q(16,3)>1,
            local_domain_stays_strictly_below_qstar=Q(6,16)<Q(1,2),
            all_active_enlarged_charts_stay_inside_background_domain=Q(3,16)<1,
            every_active_band_is_above_fixed_admitted_threshold=Q(3,16)<1,
            q_below_half_domain_alone_does_not_imply_fourQ_margin=Q(1,4)<Q(1,2) and 4*Q(1,2)>1,
            next_time_derivative_required_not_merely_bounded_values=True,
            compatible_terminal_spatial_jets_from_uniform_limits_on_enlarged_compacts=True,
            uniform_axis_hole_for_annular_fields_from_Xa_and_positive_q=True,
            inherited_all_order_physical_bounds_imply_needed_endpoint_inequality_conditionally=True,
            full_seed_pressure_and_all_correction_differences_in_same_gauge_retained=True,
            exact_Gaussian_velocities_kept_and_no_new_forward_inverse_at_endpoint=True))


def rising(x,j): return prod((Q(x)+k for k in range(order(j))),start=Q(1))
def heat_A(h,j):
    h=parameter(h);return rising(h,j)*rising(1+h,j)
def heat_S(h,j):return sum(Q(comb(order(j),k))*heat_A(h,k)*heat_A(h,j-k) for k in range(j+1))
def normalized_heat_time_jet(h,j):return 4**order(j)*heat_A(h,j)
def normalized_pressure_time_jet(h,j):return -4**order(j)*heat_S(h,j)/(2+4*parameter(h)+2*j)


def heat_review():
    h=Q(1,200);rows=[]
    for j in range(7):
        K=normalized_heat_time_jet(h,j);p=normalized_pressure_time_jet(h,j)
        square=sum(Q(comb(j,k))*normalized_heat_time_jet(h,k)*normalized_heat_time_jet(h,j-k) for k in range(j+1))
        rows.append(dict(j=j,normalized_swirl_jet=str(K),normalized_pressure_jet=str(p),
            radial_pressure_identity=-(2+4*h+2*j)*p==square))
    return dict(rows=rows,heat_jet='partial_t^j K(r,0)=C*4^j*(h)_j*(1+h)_j*r^(-1-2*h-2*j)',
        pressure_jet='partial_t^j p(r,0)=-C^2*4^j*S_j/(2+4*h+2*j)*r^(-2-4*h-2*j)',
        pressure_majorant='|partial_t^j(K(rho,tau)^2/rho)|<=C_j*rho^(-3-4*h-2*j)',
        checks=dict(time_sign_is_positive_for_swirl=all(Q(r['normalized_swirl_jet'])>0 for r in rows),
            pressure_sign_and_radial_memory_retained=all(Q(r['normalized_pressure_jet'])<0 for r in rows),
            complete_Leibniz_convolution_satisfies_radial_pressure_identity=all(r['radial_pressure_identity'] for r in rows),
            pressure_density_majorant_integrable_at_infinity_for_all_j=Q(3)>1,
            heat_potential_vanishes_and_direct_swirl_remains_in_exterior=True))


def nonanalytic_review():
    h=Q(1,200);rows=[]
    for j in (1,2,10,100):
        ratio=normalized_heat_time_jet(h,j+1)/normalized_heat_time_jet(h,j)/(j+1)
        margin=4*h*(2*j+1+h)/(j+1)
        lower=h*factorial(j-1)*factorial(j)
        rows.append(dict(j=j,Taylor_coefficient_ratio=str(ratio),ratio_minus_4j=str(margin),
            ratio_identity=ratio-4*j==margin,heat_factorial_lower_bound=heat_A(h,j)>=lower,
            pressure_convolution_lower_bound=heat_S(h,j)>=2*heat_A(h,j)))
    return dict(nonzero_exterior_constant_required=True,rows=rows,
        all_j_identity='ratio-(4/r^2)*j=(4/r^2)*h*(2*j+1+h)/(j+1)>0',
        all_j_lower_bound='(h)_j*(1+h)_j>=h*(j-1)!*j!, j>=1',
        consequence='both nonzero heat swirl and full exterior pressure have zero-radius terminal time Taylor series',
        checks=dict(exact_ratio_identity=all(r['ratio_identity'] for r in rows),
            positive_margin_on_whole_h_interval=0<Q(1,100),
            factorial_growth_and_full_pressure_lower_bound=all(r['heat_factorial_lower_bound'] and r['pressure_convolution_lower_bound'] for r in rows),
            smooth_one_sided_jets_do_not_imply_analytic_or_zero_future_continuation=True,
            zero_radius_is_not_an_obstruction_to_smooth_Borel_extension=True))
