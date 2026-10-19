"""Sharp full affine cutoff cost and the mixed current-frame scale certificate."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
from pathlib import Path
import importlib.util
import sys

HERE=Path(__file__).resolve().parent
def rat(x):
    if isinstance(x,(bool,float)):raise TypeError('Exact rational inputs required')
    return F(x)
def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [serial(v) for v in x]
    return x
@lru_cache(None)
def module(name,file):
    p=HERE.parent/name/file;s=importlib.util.spec_from_file_location('sharp_current_'+name,p)
    m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
def shape():return module('cutoff_shape_feasibility','shape_bounds.py')
def current():return module('current_frame_freezing','mixed_bounds.py')
def sign(x):return int(x>0)-int(x<0)
def quad(a,b=0,q=0):
    """a+b sqrt(q), with exact rational coefficients; normalize rational roots."""
    a,b,q=map(rat,(a,b,q))
    if q<0:raise ValueError('Nonnegative rational radicand required')
    if not b or not q:return dict(a=a,b=F(0),q=F(0))
    n,d=isqrt(q.numerator),isqrt(q.denominator)
    if n*n==q.numerator and d*d==q.denominator:return dict(a=a+b*F(n,d),b=F(0),q=F(0))
    return dict(a=a,b=b,q=q)
def qsign(x):
    a,b,q=(x[k] for k in ('a','b','q'));sa,sb=sign(a),sign(b)
    if not b or not q:return sa
    if not a:return sb
    if sa==sb:return sa
    return sa*sign(a*a-b*b*q)
def qcmp(x,y):
    if x['b'] and y['b'] and x['q']!=y['q']:raise ValueError('Comparison requires one shared quadratic field or a rational value')
    q=x['q'] if x['b'] else y['q']
    return qsign(quad(x['a']-y['a'],x['b']-y['b'],q))
def linear(x,a,b):
    a,b=map(rat,(a,b));return quad(a+b*x['a'],b*x['b'],x['q'])
def square(x):return quad(x['a']**2+x['b']**2*x['q'],2*x['a']*x['b'],x['q'])
def qmax(x,y):return x if qcmp(x,y)>=0 else y

def full_floor(rows):
    """Infimum over ALL smooth-shape enclosure pairs, not the unchanged cutoff."""
    rows=shape().normalized_rows(rows);D0,D1,_=rows['derivative'];A0,A1,A2=rows['stress']
    C0=max(rows['unit'][0],rows['value'][0]);P,Q=A1+A2,A2;S0=A0+A1+2*A2;G0=D0+D1
    ts=quad(0,1,Q/P);comparison=qsign(quad(G0-S0,D1-2*P,Q/P))
    if comparison<=0:
        branch='stress_stationary';t=ts;variable=linear(ts,S0,2*P)
    else:
        branch='derivative_stress_crossing';a=P-D1;b=G0-S0;disc=b*b-4*a*Q
        if a==0:
            if b<=0:raise ValueError('Invalid linear crossing')
            t=quad(Q/b)
        else:
            if disc<=0:raise ValueError('Crossing must have positive discriminant')
            t=quad(b/(2*a),-1/(2*a),disc)
        variable=linear(t,G0,D1)
    dominated=qcmp(quad(C0),variable)>0
    infimum=quad(C0) if dominated else variable
    return dict(rows=rows,P=P,Q=Q,S0=S0,G0=G0,constant_floor=C0,
        stationary_branch_comparison=comparison,active_variable_branch=branch,
        closure_t=t,variable_floor=variable,infimum=infimum,
        infimum_attained_by_some_smooth_shape=dominated,
        strict_floor_for_every_smooth_shape=not dominated)

def decision(rows,B):
    B=rat(B);floor=full_floor(rows);cmp=qcmp(quad(B),floor['infimum'])
    possible=cmp>0 or (cmp==0 and floor['infimum_attained_by_some_smooth_shape'])
    return dict(possible=possible,budget=B,budget_minus_infimum_sign=cmp,
        infimum=floor['infimum'],infimum_attained=floor['infimum_attained_by_some_smooth_shape'])

def checked_cost(cost):
    out={k:rat(cost[k]) for k in ('Cmesh','Cpert','chart_distance_threshold','profile_normalization')}
    if out['profile_normalization']!=8 or out['Cpert']<=0 or min(out['Cmesh'],out['chart_distance_threshold'])<0:
        raise ValueError('Retained nonnegative mesh/geometry and positive mixed value costs required')
    return out

def mixed_budget(data,cost,slack,S,ell):
    """S is the supplied error-decay factor E<=Cstar/S, not necessarily an exact physical power."""
    S=rat(S);slack=rat(slack)
    if S<=0 or type(ell) is not int or ell<2:raise ValueError('Positive decay factor and integer ell>=2 required')
    cost=checked_cost(cost)
    cert=current().deficit().certificate(data,slack);gamma=cert['floors']['unfrozen_cone_margin'];w=data['inputs']['omega']
    headroom=gamma/8-cost['Cmesh']/F(ell**6)
    caps=dict(base=S*data['base_radius']/8,cone=S*w*w*slack/(16*data['Lpsi']),
        freezing=S*headroom/(8*cost['Cpert']))
    return dict(decay_factor=S,slack=slack,unfrozen_margin=gamma,mesh_headroom=headroom,
        caps=caps,effective_budget=min(caps.values()),
        chart_gate=F(ell**6)>=cost['chart_distance_threshold'])
def joint_test(data,cost,rows,slack,S,ell):
    budget=mixed_budget(data,cost,slack,S,ell);r=decision(rows,budget['effective_budget'])
    return dict(budget=budget,shape_decision=r,
        possible=budget['chart_gate'] and budget['mesh_headroom']>0 and r['possible'])
def scale_floor(data,cost,rows,S,ell):
    S=rat(S)
    if S<=0 or type(ell) is not int or ell<2:raise ValueError('Positive decay factor and integer ell>=2 required')
    cost=checked_cost(cost)
    opt=full_floor(rows);C=opt['infimum'];w=data['inputs']['omega'];k=data['kappa']
    background=linear(C,0,16*data['Lpsi']/(w*w*S))
    freezing=linear(C,72*k*cost['Cmesh']/(w*w*ell**6),576*k*cost['Cpert']/(w*w*S))
    return dict(cost_floor=C,profile_error_floor=linear(C,0,8/S),background_deficit_floor=background,
        freezing_deficit_floor=freezing,required_deficit_floor=qmax(background,freezing),
        attainable_by_some_smooth_shape=opt['infimum_attained_by_some_smooth_shape'],
        strict_for_every_smooth_shape=opt['strict_floor_for_every_smooth_shape'])

def toy(D0=0,D1=1,A0=0,A1=3,A2=1,C0=1):
    return dict(unit=(F(1),F(0),F(0)),value=(rat(C0),F(0),F(0)),
        derivative=(rat(D0),rat(D1),F(0)),stress=tuple(map(rat,(A0,A1,A2))))
def fixtures():
    return dict(stress_stationary=toy(),derivative_crossing=toy(D1=3,A1=0),
        linear_crossing=toy(D0=3,D1=1,A1=0),irrational_stress=toy(A1=1),
        irrational_crossing=toy(D0=1,D1=7),positive_quadratic_crossing=toy(D0=10),constant_attained=toy(C0=10),
        constant_tied_unattained=toy(C0=9),stationary_equal_derivative=toy(D0=F(9,2),D1=3),
        endpoint_stationary=toy(A1=0))

def as_symbolic(x):
    import sympy as s
    return s.Rational(x['a'].numerator,x['a'].denominator)+s.Rational(x['b'].numerator,x['b'].denominator)*s.sqrt(s.Rational(x['q'].numerator,x['q'].denominator))

def optimizer_review():
    import sympy as s
    out=[];budgets=[]
    for name,rows in fixtures().items():
        f=full_floor(rows);t=as_symbolic(f['closure_t']);v=as_symbolic(f['variable_floor'])
        fs=f['S0']+f['P']*t+f['Q']/t;gs=f['G0']+rows['derivative'][1]*t
        out.append(dict(name=name,optimization=serial(f),checks=dict(
            closure_parameter_positive=qsign(f['closure_t'])>0,
            closure_parameter_at_or_before_stress_stationary=qcmp(square(f['closure_t']),quad(f['Q']/f['P']))<=0,
            exact_stress_value_at_optimizer=s.simplify(fs-v)==0,
            derivative_value_at_or_below_variable_floor=qcmp(linear(f['closure_t'],f['G0'],rows['derivative'][1]),f['variable_floor'])<=0,
            correct_stationary_or_crossing_equation=s.simplify(f['P']-f['Q']/t**2)==0 if f['active_variable_branch']=='stress_stationary' else s.simplify(fs-gs)==0,
            equality_attainment_iff_constant_strictly_dominates=f['infimum_attained_by_some_smooth_shape']==(qcmp(quad(f['constant_floor']),f['variable_floor'])>0))))
        Bs={F(k,2) for k in range(1,41)}
        if f['infimum']['b']==0:Bs.update(f['infimum']['a']+e for e in (-F(1,1000),F(0),F(1,1000)))
        for B in sorted(Bs):
            new=decision(rows,B);old=shape().shape_gate(rows,B)
            budgets.append(dict(name=name,B=str(B),possible=new['possible'],checks=dict(
                exact_closed_form_agrees_with_independent_capped_quadratic=new['possible']==old['possible'],
                accepted_witness_capacity_and_original_rows=not old['possible'] or all(old['checks'].values()))))
    return dict(auxiliary_rows_not_actual_cutoff_or_profile=True,regimes=out,budget_controls=budgets,
        checks=dict(all_ten_exact_optimizer_regimes=all(all(r['checks'].values()) for r in out),
            all_closed_form_decisions_match_original_quadratic=all(all(r['checks'].values()) for r in budgets),
            derivative_row_can_raise_floor_above_stress_infimum=full_floor(toy(D1=3,A1=0))['infimum']==quad(F(9,2)),
            constant_row_minimum_can_be_attained=decision(toy(C0=10),10)['possible'],
            tied_constant_row_does_not_attain_smooth_limit=not decision(toy(C0=9),9)['possible']))

def symbolic_review():
    import sympy as s
    t,P,Q,S0,G0,D1,b,a,C,Cm,Cp,w,k,S,gamma=s.symbols('t P Q S0 G0 D1 b a C Cm Cp w k S gamma',positive=True)
    f=S0+P*t+Q/t;g=G0+D1*t
    return dict(checks=dict(
        exact_derivative_stress_crossing_polynomial=s.simplify(t*(f-g)-((P-D1)*t*t+(S0-G0)*t+Q))==0,
        crossing_root_rationalization=s.simplify(2*Q/(b+s.sqrt(b*b-4*a*Q))-(b-s.sqrt(b*b-4*a*Q))/(2*a))==0,
        stationary_stress_minimum=s.simplify(f.subs(t,s.sqrt(Q/P))-S0-2*s.sqrt(P*Q))==0,
        mixed_budget_headroom_solved_without_dropping_value_cost=s.simplify((gamma/8-Cm)/(8*Cp)*S-(S*(gamma-8*Cm)/(64*Cp)))==0,
        freezing_raw_deficit_keeps_mesh_and_value_terms=s.simplify(72*k/w**2*(Cm+8*Cp*C/S)-(72*k*Cm/w**2+576*k*Cp*C/(w**2*S)))==0,
        stress_cost_strictly_increases_with_second_enclosure=s.diff(S0+P*t+Q/t+Q*s.Symbol('eps'),s.Symbol('eps'))==Q))

def coupled_review():
    cur=current();data=cur.deficit().constants(F(1,2),2,2,6,4,F(1,2),F(3,4));slack=F(1,100)
    cost=cur.current_costs(data,0,0,6,F(1,20000));gamma=cur.deficit().certificate(data,slack)['floors']['unfrozen_cone_margin']
    rows=toy();ell=10002;S9=9*64*cost['Cpert']/gamma;S10=10*64*cost['Cpert']/gamma
    at=joint_test(data,cost,rows,slack,S9,ell);after=joint_test(data,cost,rows,slack,S10,ell)
    old_budget=S9*cur.deficit().certificate(data,slack)['radius']/8
    constant=joint_test(data,cost,toy(C0=10),slack,S10,ell)
    old=shape().shape_gate(rows,old_budget)
    best=scale_floor(data,cost,rows,S9,ell)
    return dict(auxiliary_scale_factors_not_selected_physical_bands=True,boundary=serial(at),strictly_after=serial(after),
        attainable_constant_boundary=serial(constant),best_case=serial(best),checks=dict(
            mixed_value_freezing_budget_is_exactly_nine=at['budget']['effective_budget']==9,
            old_background_shape_and_mesh_tests_pass_but_mixed_test_fails=old['possible'] and cost['Cmesh']==0 and not at['possible'],
            every_positive_full_cost_slack_can_pass_here=after['possible'],
            constant_dominating_boundary_can_be_attained=constant['possible'] and constant['budget']['effective_budget']==10,
            raw_deficit_floor_equals_supplied_slack_but_is_unattained=qcmp(best['required_deficit_floor'],quad(slack))==0 and best['strict_for_every_smooth_shape']))

def full_composition_review():
    cur=current();g=shape().gate_module();p=g.fixture();rows=g.affine_gate(**p)['rows']
    data=cur.deficit().constants(F(1,2),2,2,6,4,F(1,2),F(3,4));chart=cur.leading().chart_caps(1,4);jets=cur.leading().join_review()['jets']
    cost=cur.current_costs(data,F(jets['slow_F_Lipschitz']),F(jets['slow_g_Lipschitz']),6,chart['admissible_slow_distance'])
    witness=cur.later_witness(data,cost,rows,32,1216,F(1,100),p['h']);ell=witness['sufficient_ell'];scale=cur.deficit().dyadic_scale(p['h'],ell)
    joined=joint_test(data,cost,rows,F(1,100),scale['factor'],ell);B=joined['budget']['effective_budget'];exist=shape().shape_gate(rows,B)
    M,N=exist['M'],exist['N'];C=g.cstar(rows,M,N);direct=cur.mixed_gate(data,cost,C/scale['factor'],ell)
    erased=g.affine_gate(**dict(p,leading=(*p['leading'][:5],F(0),F(0),F(0))))['rows']
    floor=scale_floor(data,cost,rows,scale['factor'],ell)
    return dict(auxiliary_full_pressure_and_global_jet_inputs_not_actual_profile=True,rows=serial(rows),optimization=serial(full_floor(rows)),
        conservative_dyadic=serial(scale),joint=serial(joined),auxiliary_witness=serial(exist),best_case=serial(floor),checks=dict(
            complete_current_frame_and_shape_join_passes=joined['possible'] and exist['possible'],
            rational_smooth_witness_passes_original_current_frame_gate=all(exist['checks'].values()) and direct['background_base_gate'] and direct['chart_gate'] and direct['required_deficit']<=F(1,100),
            cutoff_shape_is_only_an_auxiliary_existence_witness=shape().capacity(M,N)>1,
            full_pressure_and_shifted_viscosity_costs_retained=rows['stress'][0]>erased['stress'][0] and rows['stress'][1:]==erased['stress'][1:],
            irreducible_deficit_below_supplied_slack_at_this_later_witness=qcmp(floor['required_deficit_floor'],quad(F(1,100)))<0,
            exact_integer_floor_whole_chart_shift_retained=scale['exact_exponent']==2*p['h']*(ell-2) and scale['factor']==F(2)**scale['integer_floor']))

def review():
    o,s,c,f=optimizer_review(),symbolic_review(),coupled_review(),full_composition_review()
    return dict(optimizer=o,symbolic=s,coupled=c,full_composition=f,checks=dict(
        sharp_full_maximum_floor_and_attainment_regimes=all(o['checks'].values()),
        exact_stationary_crossing_and_mixed_budget_identities=all(s['checks'].values()),
        coupled_fixed_scale_obstruction_and_attainable_constant_case=all(c['checks'].values()),
        complete_global_jet_pressure_and_shape_scale_composition=all(f['checks'].values())))
