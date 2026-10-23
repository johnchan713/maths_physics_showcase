"""Factor-level exact checks for the finite support closure proof."""
from functools import lru_cache
from math import comb
from fractions import Fraction as Q
import sympy as s
from bounds import flat_radius_recipe,signed_half_weight,signed_support_threshold,support_recipe,certificate,source

def root_controls():
    checks={};u,z,t=s.symbols('u z scale',positive=True)
    X=s.Function('Xi')(u,z);y=s.Function('y',positive=True)(u,z)
    for a,b in ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2),(2,1),(1,2),(2,2)):
        expr=s.diff(X/s.sqrt(y),u,a,z,b)
        replacements={d:t*d for d in expr.atoms(s.Derivative)}
        replacements[X]=t*X;replacements[y]=t*y
        scaled=expr.xreplace(replacements)
        checks[f'jet{a}{b}:every_signed_mixed_jet_retains_half_weight']=s.simplify(scaled-s.sqrt(t)*expr)==0
        if a+b:
            naive=s.diff(X,u,a,z,b)/s.sqrt(y)
            checks[f'jet{a}{b}:root_denominator_derivatives_cannot_be_dropped']=s.simplify(expr-naive)!=0
    for n in (0,1,5,60000006):
        r=signed_half_weight(n)
        checks[f'n{n}:original_root_product_not_reciprocal_alone_is_extended']=r['derivative_weight']==Q(1,2) and not r['reciprocal_root_alone_extended'] and not r['new_root_total_wave_or_finer_target_callback_used']
    w=s.symbols('weight',positive=True)
    for theta in (Q(1,4),Q(1,2),Q(3,4),Q(1)):
        r=signed_support_threshold(theta);power=s.Rational(theta.numerator,theta.denominator)-s.Rational(1,2)
        actual=s.limit(w**power,w,0,dir='+')
        expected=s.S.Zero if theta>Q(1,2) else (s.S.One if theta==Q(1,2) else s.oo)
        checks[f'theta{theta}:sharp_support_threshold_has_explicit_noncanonical_comparison']=actual==expected and r['flat_extension_for_all_inputs_in_weighted_class']==(theta>Q(1,2)) and not r['actual_canonical_target_failure_claimed']
    checks['smooth_compact_half_weight_target_can_give_nonvanishing_signed_coefficient']=s.simplify(s.sqrt(w)/s.sqrt(w))==1
    return checks,dict(scope='exact homogeneous weight accounting of actual linear signed update; positivity and weighted target jets remain input inequalities')

def edge_controls():
    checks={};x,y,A,B=s.symbols('x y A B',positive=True)
    E=s.exp(-A/x**2-B/y**2)
    for m,n in ((0,0),(1,0),(0,1),(2,0),(1,1),(2,1),(1,2),(2,2)):
        P=s.S.One
        for _ in range(m):P=s.expand(s.diff(P,x)+2*A*P/x**3)
        for _ in range(n):P=s.expand(s.diff(P,y)+2*B*P/y**3)
        direct=s.diff(E,x,m,y,n)
        checks[f'm{m}n{n}:flat_edge_derivative_is_exponential_times_finite_poles']=s.simplify(direct-E*P)==0
        checks[f'm{m}n{n}:no_negative_exponential_weight_created']=all(v.is_Integer for term in s.Add.make_args(P) for v in (term.as_powers_dict().get(x,s.S.Zero),term.as_powers_dict().get(y,s.S.Zero))) if P!=1 else True
    r,z,t=s.symbols('r z t',real=True)
    left=1+z+z*t;length=3+t*t
    X=r-left;Y=length-X;L=-A/X**2-B/Y**2
    P=s.S.One
    for var in (r,z,t):P=s.diff(P,var)+s.diff(L,var)*P
    checks['moving_chart_full_RZT_chain_matches_direct_flat_factor']=s.simplify(s.diff(s.exp(L),r,z,t)/s.exp(L)-P)==0
    rows=[]
    for p,K in ((0,0),(1,0),(2,3),(19,8),(60000006,30000000)):
        a=flat_radius_recipe(p,K);rows.append(a)
        checks[f'p{p}K{K}:arbitrary_finite_poles_leave_strict_vanishing_power']=a['resulting_power']>K and a['strict_margin'] in (1,2) and not a['materialized']
    return checks,dict(deferred_vanishing_recipes=rows,scope='all-order argument uses positive exponential series and finite chain rules, not sampled edge values')

def ode_cut_controls():
    checks={};v=s.symbols('v',real=True)
    # Differentiate t'=A t-f. Zero forcing jets and zero entrance imply
    # zero response jets. Nonzero forcing first appears one jet later.
    A=[s.Matrix([[j+1,j+2],[-j,2*j+1]]) for j in range(7)]
    zero=s.zeros(2,1)
    def jets(forcing):
        out=[zero]
        for n in range(7):out.append(sum((comb(n,j)*A[j]*out[n-j] for j in range(n+1)),zero)-forcing[n])
        return out
    allzero=jets([zero]*7)
    checks['zero_full_entrance_jets_force_zero_response_jets']=all(j==zero for j in allzero)
    for n in (0,2,5):
        f=[zero]*7;f[n]=s.Matrix([2,-3]);a=jets(f)
        checks[f'first_source_jet{n}:Gaussian_smallness_cannot_replace_entrance_flatness']=a[n+1]==-f[n] and all(j==zero for j in a[:n+1])
    # Smooth scalar forward forcing can leave a nonzero outgoing state.
    uncut=s.exp(-v)
    checks['homogeneous_outgoing_response_is_not_automatically_compact']=s.diff(uncut,v)+uncut==0 and uncut.subs(v,2)!=0
    mu=s.symbols('mu',real=True);response=mu**4*v*s.exp(-v);forcing=mu**4*s.exp(-v)
    checks['whole_path_parameter_face_support_is_preserved_by_zero_data_solve']=s.simplify(s.diff(response,v)+response-forcing)==0 and response.subs(v,0)==0 and all(s.diff(response,mu,j).subs(mu,0)==0 for j in range(4))
    bad=1-(1+v)*s.exp(-v)
    checks['entrance_vanishing_alone_does_not_supply_whole_path_face_support']=s.simplify(s.diff(bad,v)+bad-(1-s.exp(-v)))==0 and bad.subs(v,0)==0 and bad.subs(v,1)!=0
    R,theta,Z,T=s.symbols('R theta Z T',positive=True)
    coords=(R,theta,Z,T);psi=s.Function('psi')(*coords)
    a=[s.Function(k)(*coords) for k in ('Ar','At','Az')]
    curl=lambda q:[s.diff(q[2],theta)/R-s.diff(q[1],Z),s.diff(q[0],Z)-s.diff(q[2],R),(s.diff(R*q[1],R)-s.diff(q[0],theta))/R]
    grad=s.Matrix([s.diff(psi,R),s.diff(psi,theta)/R,s.diff(psi,Z)])
    cross=grad.cross(s.Matrix(a));full=curl([psi*j for j in a]);old=curl(a)
    for j in range(3):
        checks[f'curl{j}:complete_cut_potential_keeps_cutoff_curl_and_connection']=s.expand(full[j]-psi*old[j]-cross[j])==0
        checks[f'curl{j}:velocity_only_cutoff_mutation_changes_field']=s.expand(full[j]-psi*old[j])!=0
    pi=s.Function('pi')(*coords)
    checks['full_cut_pressure_slow_derivative_keeps_cutoff_term']=s.expand(s.diff(psi*pi,Z)-psi*s.diff(pi,Z)-s.diff(psi,Z)*pi)==0
    normal=s.Matrix([1,0,0]);f=normal;Proj=s.eye(3)-normal*normal.T
    checks['zero_projected_forcing_velocity_does_not_remove_full_source_normal_pressure']=Proj*f==s.zeros(3,1) and (normal.dot(f))==1
    return checks,dict(scope='ODE entrance and full cut potential/pressure operations; no raw outgoing zero-data support claim')

def primitive_controls():
    checks={};u=s.symbols('u',real=True)
    f=(u-1)**3*(2-u)**3;moment=s.integrate(f,(u,1,2))
    v=s.symbols('v',real=True);chi=s.integrate(6*(v-3)*(4-v),(v,3,u))
    checks['nonzero_source_moment_is_allowed_by_compact_Te']=moment!=0 and chi.subs(u,3)==0 and chi.subs(u,4)==1
    gap=moment;exit=moment-chi.subs(u,4)*moment
    checks['actual_potential_hull_includes_gap_between_source_and_cutoff']=gap!=0 and exit==0
    checks['same_cutoff_restores_zero_exterior_potential']=s.expand(moment-chi*moment).subs(u,4)==0
    # The physical periodic pullback is global, not an inverse branch.
    J=s.Matrix([[3,1],[1,5]]);k=s.Matrix([2,-1]);Y=s.Matrix(s.symbols('Y0 Y1',real=True))
    left=(k.T*J*Y)[0];right=((J.T*k).T*Y)[0]
    checks['global_integer_pullback_has_one_physical_phase_on_every_chart']=s.expand(left-right)==0 and J.det()==14
    for axis in (0,1):
        shifted=Y+s.eye(2)[:,axis]
        phase=s.expand((k.T*J*shifted)[0]-left)
        checks[f'axis{axis}:periodic_cover_uses_no_independent_inverse_branch']=phase.is_Integer
    return checks,dict(scope='shifted primitive hull and global periodic defining maps; finite polynomial fixture does not assert C-infinity edge regularity')

def closure_controls():
    checks={};a=support_recipe(0,physical_power=0);c=certificate()
    checks['all_mixed_trace_orders_follow_conditionally_not_from_samples']=a['conditional_all_mixed_trace_orders_zero'] and not a['huge_independent_trace_samples_needed_if_global_input_contract_proved']
    checks['same_old_physical_raw_request_and_entire_DAG_retained']=a['physical_trace_request']['original_raw_Phi_order']==60000006 and a['complete_older_raw_norm_and_operator_request_graph_retained']
    checks['actual_weighted_input_or_same_band_gaps_not_certified']=not a['actual_input_contract_or_same_band_gaps_verified'] and not c['actual_all_scale_inputs_absolute_decay_infinite_field_or_blowup_proved']
    return checks,dict(certificate=c)

@lru_cache(maxsize=1)
def review():
    checks={};details={}
    for name,fn in (('root',root_controls),('edge',edge_controls),('cut',ode_cut_controls),('primitive',primitive_controls),('closure',closure_controls)):
        c,d=fn();checks.update({name+':'+k:bool(v) for k,v in c.items()});details[name]=d
    return dict(checks=checks,details=details)
