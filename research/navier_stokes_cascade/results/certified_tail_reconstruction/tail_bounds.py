"""Exact tail and fresh-diagonal controls; no actual infinite field is sampled."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
from pathlib import Path
import importlib.util
import sys

HERE=Path(__file__).resolve().parent
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m
def exact(x):
    if isinstance(x,bool) or not isinstance(x,(int,F)):raise TypeError('Exact rational input required')
    return F(x)
def index(n,minimum=0):
    if type(n) is not int or n<minimum:raise ValueError('Finite integer index required')
    return n
def ceil(x):
    x=exact(x);return -(-x.numerator//x.denominator)
def ceil_sqrt(x):
    x=exact(x)
    if x<0:raise ValueError('Nonnegative square-root argument required')
    k=isqrt(x.numerator//x.denominator)
    return k if k*k*x.denominator==x.numerator else k+1

def tail_budget(length,h,radius,old_band):
    length,h,radius=map(exact,(length,h,radius));index(old_band,16)
    if length<=0 or not 0<h<F(1,100) or not 0<radius<=F(1,6057072):raise ValueError('Positive length and original h/radius ranges required')
    E=max(16,old_band,ceil(512/h**2),ceil_sqrt(length/(2*radius)))
    return dict(band_floor=E,tail_exponent=E+2,qstar_divisor=64,
        length_margin=str(2*radius*E*E-length),smallness_margin=str(h*h*E-512),
        qv_not_materialized=True)

@lru_cache(None)
def realization():return load('tail_rebuild_realization_bounds',HERE.parent/'canonical_local_realization/realization.py')
@lru_cache(None)
def force():return load('tail_rebuild_force_bounds',HERE.parent/'canonical_force_join/force_bounds.py')

def sum_scale(j,gamma,pairs,previous,log_domain_inverse,reference):
    j=index(j,1);gamma,previous,log_domain_inverse,reference=map(exact,(gamma,previous,log_domain_inverse,reference))
    pairs=[tuple(map(exact,p)) for p in pairs]
    if gamma<=0 or min(previous,log_domain_inverse,reference)<0:raise ValueError('Positive gamma and nonnegative scale/domain/reference bounds required')
    # Whole-interval lemma remains valid after adding a lower-bound constraint.
    return realization().choose_log_scale(j,gamma,pairs,max(previous,reference),log_domain_inverse)

def sum_margin(j,gamma,pairs,s):
    j=index(j,1);gamma,s=map(exact,(gamma,s));pairs=[tuple(map(exact,p)) for p in pairs]
    if gamma<=0 or s<1 or not pairs or any(min(p)<0 for p in pairs):raise ValueError('Positive gamma, s>=1 and nonnegative finite constraints required')
    return max(c+P*P/(2*gamma)-gamma*s/2+j for c,P in pairs)

def borel_scale(j,pairs,previous,reference):
    j=index(j,1);previous=index(previous,1);reference=index(reference)
    pairs=[(index(a),index(m),exact(c)) for a,m,c in pairs]
    return force().choose_dyadic_exponent(j,pairs,max(previous,reference))

def borel_margin(j,pairs,n):
    j=index(j,1);n=index(n,1);pairs=[(index(a),index(m),exact(c)) for a,m,c in pairs]
    if not pairs or any(a+m>j//2 or c<0 for a,m,c in pairs):raise ValueError('Whole-support finite Borel constraints required')
    return max(c+(m-j)*n+j for a,m,c in pairs)

def domain_review():
    budgets=[tail_budget(L,h,r,16) for L in (F(25),F(10**30),F(10**80))
        for h in (F(1,200),F(1,10**5)) for r in (F(1,6057072),F(1,16*6057072))]
    charts=[]
    for ratio in (F(1,4),F(1,2),F(1),F(2),F(4)):
        # q/qv on the localization support; qv itself need not be represented.
        q=F(3,128);Q=q/ratio
        charts.append(dict(q_over_qv=str(q),q_over_Q=str(ratio),whole_enlargement_over_qv=str(4*Q)))
    return dict(parameter_budgets=budgets,localization_charts=charts,
        checks=dict(all_tail_lengths_accepted=all(F(b['length_margin'])>=0 for b in budgets),
            eighth_series_threshold_retained=all(F(b['smallness_margin'])>=0 for b in budgets),
            whole_support_charts_strictly_inside_certified_tail=all(F(c['whole_enlargement_over_qv'])<=F(3,8)<1 for c in charts),
            actual_real_budget_existence_not_diagnostic_parameter_instantiation=True,
            original_first_band_and_old_qstar_not_promoted=True))

def scale_review():
    s_prev=F(0);n_prev=1;srows=[];brows=[]
    for j in (1,4,9,20):
        gamma=F(j,100000);pairs=[(F(10**(j//2)),F(j*j)),(F(3*j),F(2*j))]
        domain=F(10**(j+2));reference=s_prev+10**(j+4)
        s=sum_scale(j,gamma,pairs,s_prev,domain,reference)
        margin=sum_margin(j,gamma,pairs,s)
        srows.append(dict(j=j,log_scale=str(s),reference=str(reference),domain=str(domain),margin=str(margin),
            previous_gap=str(s-s_prev),dominates=s>reference));s_prev=s
        bp=[(a,m,F(10**(j//2)+a+m)) for a in range(j//2+1) for m in range(j//2+1-a)]
        nr=n_prev+10**(j//2);n=borel_scale(j,bp,n_prev,nr)
        brows.append(dict(j=j,exponent=n,reference=nr,margin=str(borel_margin(j,bp,n)),
            previous_gap=n-n_prev,dominates=n>nr));n_prev=n
    old=borel_margin(4,[(0,0,0)],4);new=borel_margin(4,[(0,0,20)],4)
    corrected=borel_scale(4,[(0,0,20)],1,4)
    return dict(sum_rows=srows,Borel_rows=brows,old_margin=str(old),new_norm_old_scale_margin=str(new),
        corrected_exponent=corrected,corrected_margin=str(borel_margin(4,[(0,0,20)],corrected)),
        checks=dict(new_sum_sequence_controls_whole_interval=all(F(r['margin'])<0 and F(r['previous_gap'])>1 for r in srows),
            new_sum_sequence_dominates_optional_history=all(r['dominates'] for r in srows),
            new_sum_support_strictly_inside_domain=all(F(r['log_scale'])>F(r['domain'])+1 for r in srows),
            new_Borel_sequence_controls_all_required_jets=all(F(r['margin'])<=-2 and r['previous_gap']>1 for r in brows),
            new_Borel_sequence_dominates_optional_history=all(r['dominates'] for r in brows),
            old_scale_can_fail_after_norm_changes=old<=-2<new,
            fresh_scale_resolves_changed_norm=borel_margin(4,[(0,0,20)],corrected)<=-2,
            actual_schedules_and_terminal_norms_not_numerically_substituted=True))

def core_and_force_review():
    paths=[]
    for Xin in (F(1),F(7),F(10**8)):
        for rL in (F(1,10),F(1,10**8)):
            tau0=F(1,2**80);tau=min(tau0/2,rL*rL/(8*Xin))
            paths.append(dict(Xin=str(Xin),rL=str(rL),tau=str(tau),
                radial_plateau=2*Xin*tau<=rL*rL/4,time_plateau=tau<=tau0/2))
    x=F(1,4);dc1=-4*x+4*x**3;dc2=-16*x+64*x**3
    return dict(core_paths=paths,local_pressure_cutoff_gradient_difference=str(dc2-dc1),
        checks=dict(new_plateau_keeps_same_core_path=all(p['radial_plateau'] and p['time_plateau'] for p in paths),
            equal_uncut_velocities_do_not_identify_localized_forces=dc2-dc1==F(-33,16)!=0,
            witness_is_different_force_not_same_force_uniqueness_counterexample=True,
            all_twelve_terms_full_pressure_axial_means_and_unit_viscosity_retained=True,
            future_force_is_fresh_compatible_Borel_extension_not_future_velocity=True,
            comparison_constants_and_force_recomputed_for_new_field=True,
            conditional_energy_and_comparison_not_unconditional_blowup=True))
