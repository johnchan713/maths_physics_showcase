"""Exact source identities, cutoff jets, supplied schedules and Cphi dependence."""
import ast
from copy import deepcopy
from fractions import Fraction as Q
from functools import lru_cache
from budgets import HERE,load,gaussian_polynomial,gaussian_raw,step_bounds,cutoff_bounds,stirling,verify_schedule,ceil_log2,comparison_dependence,comparison_value,selector_gate,fixture


@lru_cache(None)
def maps():return load('borel_certificate_maps',HERE.parent/'modulation_repair_all_order/map_checks.py')


def source_review():
    schedule=ast.parse((HERE.parent/'outer_pressure_pilot/schedule.py').read_text())
    checks_ast=ast.parse((HERE.parent/'background_summation/checks.py').read_text())
    step=next(n for n in schedule.body if isinstance(n,ast.FunctionDef) and n.name=='step')
    cutoff=next(n for n in checks_ast.body if isinstance(n,ast.FunctionDef) and n.name=='cutoff_derivatives')
    dump=lambda n:ast.dump(n,include_attributes=False)
    expected=ast.parse('1/(1+mp.exp(1/x**2-1/(1-x)**2))',mode='eval').body
    chi=ast.parse('lambda v: 1-step(2*v-1)',mode='eval').body
    return dict(checks=dict(
        implemented_step_body_is_the_squared_flat_logistic=dump(step.body[-1].value)==dump(expected),
        implemented_Borel_cutoff_is_exact_same_affine_step=any(dump(n)==dump(chi) for n in ast.walk(cutoff)),
        prior_numeric_cutoff_helper_is_not_used_as_an_interval_certificate=any(isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='float' for n in ast.walk(cutoff))),
        scope='AST identity of the frozen implemented functions; analytic real bounds, not their floating outputs, are used')


def jet_review():
    module=maps();J=module.JetClass();checks={};rows=[]
    for x0 in (Q(1),Q(1,2),Q(1,3)):
        x=J.variable(0,x0,n=1,order=11)
        normalized=module.exponential(1/(x0*x0)-x**-2)
        for j in range(8):
            value=sum(Q(c)*(1/x0)**d for d,c in gaussian_polynomial(j).items())
            checks[f'x{x0}_j{j}_Gaussian_raw_recurrence_identity']=normalized.partial((j,))==value
    s=J.variable(0,Q(3,4),n=1,order=11);x=2*s-1
    chi=1-(1+module.exponential(x**-2-(1-x)**-2))**-1
    all_bounds=cutoff_bounds(7);z=chi
    for j in range(8):
        val=z.value();checks[f'j{j}_exact_Euler_midpoint_jet_bounded']=abs(val)<=all_bounds['K'][j]
        rows.append(dict(order=j,midpoint_Euler_jet=str(val),K=str(all_bounds['K'][j])))
        if j<7:z=s*z.derivative(0)
    checks.update(midpoint_first_and_second_Euler_jets_are_minus_twelve=rows[1]['midpoint_Euler_jet']=='-12' and rows[2]['midpoint_Euler_jet']=='-12',
                  ordinary_second_jet_cannot_replace_Euler_second=chi.partial((2,))==0 and rows[2]['midpoint_Euler_jet']!='0',
                  padded_K2_three_cannot_bound_this_implemented_cutoff=Q(3)<12,
                  certified_low_Euler_bounds_are_1_64_4928=cutoff_bounds(2)['K']==[1,64,4928])
    for N in (1,4,16,64):
        width=Q(1,4*N);left=Q(3,4)-width/2
        checks[f'N{N}_qualitative_plateaus_do_not_bound_cutoff_norms']=Q(1,2)<=left<left+width<=1 and Q(6)/width==24*N
    return dict(checks=checks,rows=rows,scope='exact rational normalized exponential jets and actual implemented-cutoff midpoint jets; global bounds follow the analytic Gaussian quotient proof')


def recurrence_review():
    checks={};records=[]
    for m in range(8):
        s=step_bounds(m);b=cutoff_bounds(m)
        checks[f'm{m}_all_fixed_step_bounds_finite_positive']=all(v>0 for v in s['raw'])
        checks[f'm{m}_padded_Euler_norms_are_monotone']=all(a<=c for a,c in zip(b['K'],b['K'][1:]))
        checks[f'm{m}_sharp_low_bounds_do_not_exceed_independent_coarse_route']=all(v<=c for v,c in zip(s['raw'],s['coarse']))
        if m:
            expected=81*sum(Q(__import__('math').comb(m,j))*2*s['gaussian'][j]*s['reciprocal'][m-j] for j in range(1,m+1))
            checks[f'm{m}_raw_binomial_reciprocal_rule_retained']=s['reciprocal'][m]==expected
        records.append(dict(order=m,Gaussian_integer_bound=str(gaussian_raw(m)),step_raw=str(s['raw'][m]),Euler_K=str(b['K'][m])))
    checks['Euler_second_requires_first_and_second_ordinary_cutoff_rows']=cutoff_bounds(2)['Euler'][2]==2*32+4*1216
    return dict(checks=checks,rows=records,scope='integer/fraction polynomial recurrences; no finite order list is promoted to a uniform-in-order constant')


def supplied_stages():
    return [dict(n=n,gain_lower=Q(n,200),log2_lower=Q(10**9*n),log2_upper=Q(10**9*n)+Q(1,2),
                 constants=tuple(Q(2*(n+m+1))*cutoff_bounds(m+1)['K'][m+1] for m in range(n+1))) for n in range(1,5)]


def schedule_review():
    stages=supplied_stages();before=deepcopy(stages);good=verify_schedule(stages)
    bad=[dict(n=1,gain_lower=Q(1,200),log2_lower=Q(2),log2_upper=Q(2),constants=(Q(1),Q(1)))]
    wrong=verify_schedule(bad)
    sep=deepcopy(stages);sep[1]['log2_lower']=sep[0]['log2_upper'];sep[1]['log2_upper']=sep[0]['log2_upper']
    checks=dict(all_supplied_normalized_prefix_rows_pass=good['all_supplied_prefix_gates'],
                no_existing_scale_is_selected_or_changed=before==stages and not good['scales_selected_or_changed'],
                prefix_is_not_infinite_schedule_certification=not good['infinite_schedule_certified'],
                inadequate_existing_scale_reports_normalized_row_failure=not wrong['rows'][0]['gates']['all_normalized_rows_on_entire_support'],
                insufficient_separation_is_not_silently_repaired=not verify_schedule(sep)['rows'][1]['gates']['separation_from_previous'],
                original_common_q_domain_gate_is_retained=not verify_schedule(stages,boundary_log2_upper=10**9+1)['rows'][0]['gates']['common_q_domain'])
    for value in (Q(1),Q(2),Q(3),Q(7,2),Q(2**100+1,3)):
        k=ceil_log2(value);checks[f'C{value}_integer_binary_log_upper']=value<=2**k and (k==0 or value>2**(k-1))
    return dict(checks=checks,stage_count=len(stages),scope='manufactured supplied binary-log intervals and constants; no actual retained sequence is numerically evaluated')


def comparison_review():
    f=fixture();r=f['radial'];o=r['result'];checks={};rows=[]
    for M in (Q(0),Q(1,1000),o['first']['through'][1]):
        dep=comparison_dependence(r['h'],o['leading']['through'][1],o['leading']['through'][2],M,r['geometry'])
        for K in (Q(1),Q(3),Q(12),Q(64),Q(4928)):
            g=dict(r['geometry']);lo=g.pop('Rmin');hi=g.pop('Rmax')
            b=load('borel_direct_head_'+str(len(rows)),HERE.parent/'core_background_finite_head_budget/budgets.py').BackgroundBudget(
                r['h'],lambda j:o['leading']['through'][j],lambda n,j:M,lambda j:K,lo,hi,'euler')
            actual=b.comparison(**g)['Cphi'];expected=comparison_value(dep,K)
            checks[f'M{M}_K{K}_all_seven_groups_equal_cutoff_dependence']=actual==expected
            rows.append(dict(M=str(M),K=str(K),Cphi=str(expected)))
    checks['certified_implemented_cutoff_bound_joins_same_radial_repair']=f['joined']['comparison']['Cphi']==comparison_value(f['dependence'],4928)
    checks['earlier_K_three_is_only_a_manufactured_bound']=f['joined']['comparison']['Cphi']>r['result']['comparison']['Cphi']
    checks['midpoint_lower_recipe_can_block_all_admissible_cutoff_upper_norms']=not f['gate']['necessary_recipe_gate']
    checks['large_supplied_length_can_pass_both_recipe_gates']=selector_gate(f['dependence'],1,2,Q(1,6057072),10**60)['sufficient_recipe_gate']
    return dict(checks=checks,rows=rows,scope='exact same-recipe substitution at manufactured full native/pressure inputs; failed sufficient certificates do not prove actual normal/cone failure')


@lru_cache(None)
def review():
    parts=dict(source=source_review(),jets=jet_review(),recurrences=recurrence_review(),schedule=schedule_review(),comparison=comparison_review())
    return dict(checks={name+':'+k:v for name,p in parts.items() for k,v in p['checks'].items()},parts=parts)
