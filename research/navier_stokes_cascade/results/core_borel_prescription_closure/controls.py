"""Universal polynomial/source checks, with explicitly manufactured log data."""
from fractions import Fraction as Q
from functools import lru_cache
from proof import HERE,load,poly,add,certificates,verify_certificate,source_rule_checks,real_stage_residuals,floor_inputs_from_rule,tangent


def algebra_review():
    checks={};rows={}
    for name,c in certificates().items():
        checks[name+'_universal_polynomial_identity_and_nonnegative_remainder']=verify_certificate(c)
        wrong=dict(c,target=add(c['target'],poly((1,()))))
        checks[name+'_wrong_target_is_rejected']=not verify_certificate(wrong)
        rows[name]={label:{'*'.join(k) or '1':str(v) for k,v in p.items()} for label,p in c.items()}
    return dict(checks=checks,certificates=rows,scope='arbitrary symbolic index and positive gains; universal algebra, not inference from finitely many stages')


def source_review():
    text=(HERE.parent/'background_jet_foundation/foundation.py').read_text()
    checks=source_rule_checks(text)
    replacements={
        'missing_all_order_row_guard':('len(log_constants) != n+1','len(log_constants) != n'),
        'wrong_gain':('mp.mpf(h.numerator)*n/h.denominator','mp.mpf(h.numerator)/h.denominator'),
        'dropped_index_budget':('v+n*mp.log(2)','v+mp.log(2)'),
        'dropped_unit_padding':('return 1+max(previous','return max(previous'),
        'dropped_domain_margin':('boundary_log+mp.log(2)','boundary_log'),
        'dropped_separation':('previous+mp.log(2)','previous'),
    }
    for name,(before,after) in replacements.items():
        assert before in text
        checks[name+'_changes_the_retained_rule']=not all(source_rule_checks(text.replace(before,after)).values())
    return dict(checks=checks,scope='AST identity of the historical diagnostic real formula; no certification of mpmath rounding or actual numerical inputs')


def stage_review():
    checks={};rows=[]
    # ell=7/10 is a rational ALGEBRA stand-in, not asserted to equal log(2).
    ell=Q(7,10);previous=Q(0);domain=Q(3)
    for n in range(1,9):
        g=Q(n,200);constants=[Q(n**3+(m+1)**2) for m in range(n+1)]
        b=1+max(previous+ell,domain+ell,max((v+n*ell)/g for v in constants))
        results=[real_stage_residuals(g,n,ell,v,previous,domain,b) for v in constants]
        checks[f'n{n}_every_finite_row_has_the_strict_gain_margin']=all(r['all_supplied_constraints_pass'] and r['normalized_log_margin']>=g for r in results)
        checks[f'n{n}_no_actual_scale_selected_or_infinite_prefix_promotion']=all(not r['scales_selected_or_changed'] and not r['infinite_schedule_verified_by_finite_inputs'] for r in results)
        rows.append(dict(n=n,gain=str(g),manufactured_log_scale=str(b),maximum_manufactured_log_constant=str(max(constants))))
        previous=b
    n=9;g=Q(n,200);bad_v=g*previous+1
    bad=real_stage_residuals(g,n,ell,bad_v,previous,domain,previous)
    checks['finite_prefix_does_not_certify_arbitrary_future_stage']=not bad['all_supplied_constraints_pass'] and bad['normalized_log_margin']<0
    old=real_stage_residuals(1,1,ell,1,0,0,3)
    enlarged=real_stage_residuals(1,1,ell,100,0,0,3)
    checks['larger_later_bound_is_not_automatically_satisfied_by_fixed_scale']=old['all_supplied_constraints_pass'] and not enlarged['all_supplied_constraints_pass']
    return dict(checks=checks,rows=rows,scope='rational parameter-space controls only; universal theorem is the arbitrary-index derivation, not this finite list')


def constant_review():
    foundation=load('prescription_foundation',HERE.parent/'background_jet_foundation/foundation.py')
    checks={}
    for n in (1,2,7):
        for m in range(n+1):
            for K in (Q(1),Q(4928),Q(10**30)):
                M=Q(n+m+1);h=Q(1,200)
                a=foundation.normalized_jet_constant(h,n,m,M,K)
                b=tangent.head.increment_constant(h,n,m,M,K)
                checks[f'n{n}_m{m}_K{K}_same_complete_tuple_constant']=a==b==max(Q(1),2*K*(1+2*n*h)**m*M)
    checks['second_Euler_order_needs_Bell_conversion']=foundation.conversion_cost(2)['frozen_to_euler']==2
    checks['first_low_order_curl_norm_still_has_extra_cutoff_derivative']=foundation.normalized_jet_constant(Q(1,200),1,1,1,100)>foundation.normalized_jet_constant(Q(1,200),1,1,1,1)
    return dict(checks=checks,scope='actual historical constant formulas on manufactured seminorms; arbitrary fixed smooth cutoff norms are retained, not bounded uniformly over all cutoffs')


def transfer_review():
    args=dict(Lead1=3,Rmin=1,Rmax=2,U=2,lambda_max=3,gmin=1,Llow=1,Lhigh=2,Dbox=6)
    r=floor_inputs_from_rule(**args);old=tangent.comparison(**args)
    checks={
        'same_Ctan_and_complete_transport_without_independent_diagonal_input':r['floor_inputs']==old,
        'independent_diagonal_is_discharged_only_for_rule_defined_construction':not r['independent_diagonal_assumption_needed_for_rule_defined_construction'] and r['all_order_diagonal_proved_relative_to_original_real_rule_and_full_coefficient_regularity'],
        'diagnostic_cutoff_identity_is_not_needed_for_this_floor':not r['diagnostic_cutoff_identity_needed_for_this_floor'],
        'actual_independent_numeric_sequence_and_band_are_not_certified':not r['actual_independently_fixed_numeric_sequence_identified'] and not r['actual_band_or_full_residual_verified'],
    }
    for t in (Q(1,100),Q(1,2),Q(1)):
        checks[f't{t}_closed_diagonal_feeds_same_head_free_first_order_sum']=tangent.diagonal_tail(t)<=t
    return dict(checks=checks,Ctan=str(old['Ctan']),retained_premises=['leading and full coefficient regularity certificate','same real max prescription with own cutoff norms and original finite constraints','actual physical chart and leading gaps','unchanged-band normal and slow inequalities'])


@lru_cache(None)
def review():
    parts=dict(algebra=algebra_review(),source=source_review(),stages=stage_review(),constants=constant_review(),transfer=transfer_review())
    return dict(checks={name+':'+k:v for name,p in parts.items() for k,v in p['checks'].items()},parts=parts)
