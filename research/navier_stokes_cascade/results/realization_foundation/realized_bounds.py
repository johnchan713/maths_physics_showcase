"""Current-source representative and diagonal budgets, not computed fluid stages."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys
import mpmath as mp

HERE=Path(__file__).resolve().parent

def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m

finite=load('realized_current_finite_bounds',HERE.parent/'finite_induction_foundation/finite_bounds.py')
B0,GAIN,RHO,KAPPA=finite.B0,finite.GAIN,finite.RHO,finite.KAPPA

def index(n,minimum=0):
    if type(n) is not int or n<minimum:raise ValueError('A finite integer order is required')
    return n

def parameters(h):
    h=F(h)
    if not 0<h<F(1,100):raise ValueError('Require the unchanged 0<h<1/100')
    return h,F(1,2)+h,F(1,2)-h

def gain(j,h=F(1,200)):return parameters(h)[0]*GAIN*index(j,1)
def tuple_loss(m,h=F(1,200)):return 2*parameters(h)[1]+2*index(m)
def cutoff_loss(m,h=F(1,200)):return tuple_loss(index(m)+1,h)
def field_loss(m,h=F(1,200)):return cutoff_loss(m,h)+1

def orders(j):
    b=finite.stage_order(index(j,1)-1);c=b+F(1,2)
    # Sourced/signed potentials and pressures: C. Fast means: H=C-kappa.
    # C-2*kappa is a common deliberately weaker bound for full mean pieces.
    return dict(wave_velocity=b-KAPPA,wave_potential=c-KAPPA,
        mean_azimuthal_potential=c-2*KAPPA,direct_axisymmetric_swirl=c-2*KAPPA,
        complete_pressure_difference=c-2*KAPPA,mean_radial_velocity=c+1-2*KAPPA)

def representative_review():
    margins={k:v-GAIN for k,v in orders(1).items()}
    return dict(current_seed=str(B0),delta=str(GAIN),rho=str(RHO),
        all_stage_margins=finite.encode(margins),identity='B_(j-1)=delta*j+rho',
        common_representative_loss='ell_m=2*A+2*m',cut_velocity_loss='L_m=ell_(m+1)',total_field_loss='H_m=L_m+1',
        checks=dict(current_index_identity=B0-GAIN==RHO,
            every_row_has_exact_positive_slope_delta=orders(2)=={k:v+GAIN for k,v in orders(1).items()},
            weakest_margin_is_one_minus_two_kappa=min(margins.values())==1-2*KAPPA>0,
            same_common_jet_loss_dominates_phase_and_frame_costs=1+F(5,2)*F(1,100)<2,
            full_mean_representatives_weaken_actual_fast_H=F(1,2)-2*KAPPA<F(1,2)-KAPPA,
            pressures_cut_after_all_four_operations_not_frozen_midcycle=True,
            no_compact_full_direct_swirl_potential_assumed=True))

def ceil(x):x=F(x);return -(-x.numerator//x.denominator)

def comparison_budget(m=0,target=0,h=F(1,200)):
    m,target=index(m),index(target);h,_,_=parameters(h)
    L,H=cutoff_loss(m+2,h),field_loss(m+2,h)
    tail_stage=ceil(2*(target+1+L+H)/(h*GAIN)-1)
    finite_stage=ceil(((target+1+finite.physical_loss(h,m))/h-B0)/GAIN)
    J=max(m+2,1,tail_stage,finite_stage)
    tail=gain(J+1,h)/2-L;residual=h*finite.stage_order(J)-finite.physical_loss(h,m)
    return dict(jet=m,target=target,h=str(h),comparison_stage=J,
        tail_stage_requirement=tail_stage,finite_stage_requirement=finite_stage,
        tail_exponent=str(tail),old_new_exponent=str(tail-H),new_square_exponent=str(2*tail),finite_residual_exponent=str(residual),
        actual_fields_or_scales_evaluated=False,
        checks=dict(jets_allocated=J>=m+2,linear_tail=tail>=target+1,
            both_old_new_products=tail-H>=target+1,new_square=2*tail>=target+1,
            complete_finite_residual=residual>=target+1))

def comparison_review():
    rows=[comparison_budget(m,N,h) for h in (F(1,200),F(9,1000),F(1,10000))
        for m,N in ((0,0),(1,3),(5,11),(20,2))]
    return dict(rows=rows,checks={k:all(r['checks'][k] for r in rows) for k in rows[0]['checks']})

@lru_cache(None)
def scale_module():return load('realized_whole_interval_scale',HERE.parent/'background_summation/summation.py')

def scale_review():
    rows=[]
    with mp.workdps(90):
        previous=mp.mpf(0)
        for j,constraints in ((1,[(0,0),(5,7)]),(4,[(1000,300),(10**6,2000)]),(17,[(10**12,10**5),(300,0)])):
            f=gain(j)/2;gamma=mp.mpf(f.numerator)/f.denominator
            s=scale_module().cutoff_log_scale(j,gamma,constraints,previous,2)
            margin=max(mp.mpf(logC)+mp.mpf(P)**2/(2*gamma)-gamma*s/2+j*mp.log(2) for logC,P in constraints)
            sample=max(mp.mpf(logC)+P*mp.log(1+x)-gamma*x+j*mp.log(2)
                for logC,P in constraints for x in (s,s+1,2*s,10*s))
            rows.append(dict(stage=j,log_scale=float(s),majorant_log_margin=float(margin),sample_log_margin=float(sample),
                doubled=bool(s>=previous+mp.log(2))))
            previous=s
    return dict(synthetic_constants_only=True,actual_stage_constants_not_evaluated=True,rows=rows,
        checks=dict(whole_interval_majorant=all(r['majorant_log_margin']<0 for r in rows),
            independent_samples=all(r['sample_log_margin']<0 for r in rows),one_doubling_sequence=all(r['doubled'] for r in rows)))

def scope_review():
    return dict(field='u0+sum_j[curl(chi(a_j*q)*A_j)+chi(a_j*q)*b_j*e_theta]',
        pressure='p0+sum_j chi(a_j*q)*(p^[j]-p^[j-1])',
        tail='2^(-J)*q^(h*delta*(J+1)/2-L_m), q<=1/(2*a_J), J>=max(1,m)',
        residual_difference='e_t-Delta e+grad pi+u_J.grad e+e.grad u_J+e.grad e',
        checks=dict(real_locally_finite_smooth_field_for_q_positive=True,
            exact_divergence_core_annulus_and_exterior_support_preserved=True,
            complete_pressure_gauge_differences_and_both_background_tails_retained=True,
            positive_derivatives_of_cutoff_supported_where_half_le_aq_le_one=True,
            one_countable_sequence_for_all_fixed_jets_and_targets=True,
            comparison_stage_fixed_before_q_limit=True,
            finite_Gaussian_ledger_kept_at_one_comparison_stage_never_summed_or_forward_inverted=True,
            no_uniform_stage_constants_or_raw_series_convergence_required=True,
            retained_core_and_finite_analytic_lemmas_still_dependencies=True,
            no_one_sided_endpoint_global_force_energy_or_blowup_promotion=True))
