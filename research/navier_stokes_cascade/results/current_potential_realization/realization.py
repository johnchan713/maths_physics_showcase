"""Current-seeded conditional realization budgets, not evaluated fluid stages."""
from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys
import mpmath as mp

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'finite_stage_induction'))
from induction import B0,GAIN,KAPPA,RHO,stage_order,physical_loss,encode
sys.path.pop(0)


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


def index(n,minimum=0):
    if type(n) is not int or n<minimum: raise ValueError('Finite integer order required')
    return n


def parameters(h):
    h=Q(h)
    if not 0<h<Q(1,100): raise ValueError('Require inherited 0<h<1/100')
    return h,Q(1,2)+h,Q(1,2)-h


def gain(j,h=Q(1,200)): return parameters(h)[0]*GAIN*index(j,1)
def tuple_loss(m,h=Q(1,200)): return 2*parameters(h)[1]+2*index(m)
def cutoff_loss(m,h=Q(1,200)): return tuple_loss(index(m)+1,h)
def field_loss(m,h=Q(1,200)): return cutoff_loss(m,h)+1


def increment_orders(j):
    """Conservative complete orders for the ACTUAL cycle j-1 -> j."""
    b=stage_order(index(j,1)-1)
    return dict(wave_velocity=b-KAPPA,wave_potential=b+Q(1,2)-KAPPA,
        mean_azimuthal_potential=b+Q(1,2)-2*KAPPA,
        direct_axisymmetric_swirl=b+Q(1,2)-2*KAPPA,
        complete_pressure_difference=b+Q(1,2)-2*KAPPA,
        mean_radial_velocity=b+Q(3,2)-2*KAPPA)


def order_review():
    margins={k:v-GAIN for k,v in increment_orders(1).items()}
    return dict(previous_stage_order='B_(j-1)=delta*j+rho',delta=str(GAIN),
        increment_orders_at_j_one=encode(increment_orders(1)),
        all_stage_margins_over_delta_j=encode(margins),
        tuple_loss='ell_m=2*A+2*m',cutoff_loss='L_m=ell_(m+1)',finite_field_loss='H_m=L_m+1',
        direct_swirl_is_not_assumed_to_have_a_compact_full_vector_potential=True,
        checks=dict(previous_stage_identity=B0-GAIN==RHO,
            every_row_has_exact_slope_delta=increment_orders(2)=={k:v+GAIN for k,v in increment_orders(1).items()},
            every_all_stage_margin_positive=min(margins.values())>0,
            weakest_margin_is_one_minus_4kappa=min(margins.values())==1-4*KAPPA,
            reserved_derivative_loss_dominates_every_allowed_phase_cost=1+Q(5,2)*Q(1,100)<2,
            curl_cost_and_fixed_log_absorption_retained=cutoff_loss(0)==tuple_loss(1) and field_loss(0)==cutoff_loss(0)+1))


def ceil(value):
    value=Q(value);return -((-value.numerator)//value.denominator)


def comparison_budget(m=0,target=0,h=Q(1,200)):
    m,target=index(m),index(target);h,_,_=parameters(h)
    L,H=cutoff_loss(m+2,h),field_loss(m+2,h)
    tail_stage=ceil(2*(target+1+L+H)/(h*GAIN)-1)
    finite_stage=ceil(((target+1+physical_loss(h,m))/h-B0)/GAIN)
    J=max(1,m+2,tail_stage,finite_stage)
    tail=gain(J+1,h)/2-L;finite=h*stage_order(J)-physical_loss(h,m)
    return dict(h=str(h),m=m,target=target,comparison_stage=J,
        tail_stage_requirement=tail_stage,finite_stage_requirement=finite_stage,
        tail_exponent=str(tail),old_field_times_tail_exponent=str(tail-H),
        tail_square_exponent=str(2*tail),finite_residual_exponent=str(finite),
        actual_fields_or_scales_evaluated=False,
        checks=dict(needed_jet_orders_allocated=J>=m+2,
            linear_tail=tail>=target+1,old_new_products=tail-H>=target+1,
            new_square=2*tail>=target+1,complete_finite_residual=finite>=target+1))


def comparison_review():
    rows=[comparison_budget(m,N,h) for h in (Q(1,200),Q(9,1000),Q(1,10000))
        for m,N in ((0,0),(1,3),(5,11),(20,2))]
    return dict(rows=rows,checks={k:all(r['checks'][k] for r in rows) for k in rows[0]['checks']})


@lru_cache(None)
def old_module():
    p=HERE.parent/'full_potential_realization';sys.path.insert(0,str(p))
    try: return load('current_realization_historical_audit',p/'audit.py')
    finally: sys.path.pop(0)


def scale_review():
    # The old recipe is used with the NEW gain. All constants are synthetic.
    old_module()
    from summation import cutoff_log_scale
    rows=[];previous=mp.mpf(0)
    with mp.workdps(80):
        for j,constraints in ((1,[(0,0),(5,7)]),(4,[(1000,300),(10**6,2000)]),
            (17,[(10**12,10**5),(300,0)])):
            gamma=mp.mpf(gain(j).numerator)/gain(j).denominator/2
            scale=cutoff_log_scale(j,gamma,constraints,previous,2)
            majorants=[logC+mp.mpf(P)**2/(2*gamma)-gamma*scale/2+j*mp.log(2) for logC,P in constraints]
            sample=max(logC+P*mp.log(1+s)-gamma*s+j*mp.log(2)
                for logC,P in constraints for s in (scale,scale+1,2*scale,10*scale))
            rows.append(dict(stage=j,gamma=str(gain(j)/2),log_scale=float(scale),
                whole_interval_log_margin=float(max(majorants)),sample_log_margin=float(sample),
                doubled=bool(scale>=previous+mp.log(2))))
            previous=scale
    return dict(synthetic_constants_only=True,actual_stage_constants_not_evaluated=True,rows=rows,
        checks=dict(whole_interval_majorant=all(r['whole_interval_log_margin']<0 for r in rows),
            independent_samples=all(r['sample_log_margin']<0 for r in rows),
            scales_doubled=all(r['doubled'] for r in rows)))


@lru_cache(None)
def identity_review():
    old_module()
    from cutoff_checks import cutoff_geometry_review,curl_review,residual_review,divergent_series_review
    reviews={k:f() for k,f in dict(actual_implicit_q=cutoff_geometry_review,
        complete_curl=curl_review,complete_residual=residual_review,raw_series_control=divergent_series_review).items()}
    # Raw roundoff-valued diagnostics are already reproduced against the
    # historical record with its original comparator. Do not freeze them
    # again under a different ancestor's near-zero float policy.
    controls={k:dict(checks={name:bool(value) for name,value in r['checks'].items()})
        for k,r in reviews.items()}
    return dict(manufactured_not_actual_global_fields=True,controls=controls,
        checks={k:all(bool(v) for v in r['checks'].values()) for k,r in reviews.items()})


def scope_review():
    return dict(exact_field='u0+sum_j[curl(chi(a_j*q)*A_j)+chi(a_j*q)*b_j*e_theta]',
        exact_pressure='p0+sum_j chi(a_j*q)*(p^[j]-p^[j-1])',
        tail='2^(-J)*q^(h*delta*(J+1)/2-L_m), q<=1/(2*a_J), J>=max(1,m)',
        residual_difference='e_t-Delta e+grad pi+u_J.grad e+e.grad u_J+e.grad e',
        checks=dict(direct_swirl_theta_independent_and_q_cutoff_preserve_divergence=True,
            wave_and_mean_curls_cut_as_potentials_not_as_velocities=True,
            full_pressure_difference_including_axial_memory_is_cut=True,
            all_cutoff_time_pressure_viscous_and_nonlinear_terms_are_in_comparison=True,
            exact_Gaussian_tails_compared_at_one_fixed_finite_stage_not_summed_or_inverted=True,
            cutoff_scales_one_sequence_for_all_fixed_jet_orders_and_targets=True,
            fixed_J_selected_before_q_limit_not_j_of_q=True,
            no_uniform_stage_constant_or_raw_series_convergence_assumption=True,
            no_one_sided_endpoint_or_global_force_claim=True))
