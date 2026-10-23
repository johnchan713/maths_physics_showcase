"""Exact piecewise polynomial checks of trace and primitive defects."""
from functools import lru_cache
from fractions import Fraction as Q
import sympy as s
from bounds import boundary_recipe,seam_orders,A1_trace_request,certificate,transported_modes,geometry

def compact_controls():
    checks={};u,z,a=s.symbols('u z a',real=True)
    h0=(u-1)**2;h1=(3-u)**2+a*(3-u)
    chi=(u-1)/2+z*(u-1)*(3-u)/8
    jump=s.expand(h1.subs(u,2)-h0.subs(u,2))
    total=s.integrate(s.diff(h0,u),(u,1,2))+s.integrate(s.diff(h1,u),(u,2,3))
    checks['whole_classical_derivative_has_negative_internal_jump']=s.expand(total+jump)==0 and jump==a
    for side,h in ((0,h0),(1,h1)):
        primitive=(s.integrate(s.diff(h0,u),(u,1,u)) if not side else
            s.integrate(s.diff(h0,u),(u,1,2))+s.integrate(s.diff(h1,u),(u,2,u)))
        actual=s.expand(primitive-chi*total)
        repaired=s.expand(h-a*side+chi*a)
        checks[f'side{side}:compact_primitive_has_cumulative_jump_and_cutoff_repair']=s.expand(actual-repaired)==0
        checks[f'side{side}:local_pressure_not_equal_to_primitive_without_trace_matching']=s.expand(actual-h)!=0
        checks[f'side{side}:smooth_value_gluing_restores_exact_primitive']=s.expand((actual-h).subs(a,0))==0
        checks[f'side{side}:moving_cutoff_derivative_of_seam_pressure_is_retained']=s.expand(s.diff(actual-h,z)-a*s.diff(chi,z))==0
    checks['compact_projection_defect_has_original_chi_R_and_sign']=s.expand(s.diff(chi,u)*total+a*s.diff(chi,u))==0
    checks['one_more_cutoff_cannot_erase_nonzero_projection_defect']=s.diff(chi,u)*total!=0
    # A constant middle piece has opposite jumps and zero classical
    # derivative. Its moment is zero but its primitive is not h=a there.
    flat_pieces=[s.S.Zero,a,s.S.Zero]
    flat_integral=sum(s.integrate(s.diff(f,u),(u,i+1,i+2)) for i,f in enumerate(flat_pieces))
    flat_jumps=[flat_pieces[i]-flat_pieces[i-1] for i in (1,2)]
    checks['zero_total_jump_does_not_remove_cumulative_primitive_defect']=sum(flat_jumps)==flat_integral==0 and s.expand(flat_integral-chi*flat_integral-flat_pieces[1])==-a
    # Fixed characteristic pullback: r^e h=H(r), so D_e h=r^-e H'.
    r,M,eta=s.symbols('r M eta',positive=True);y=s.symbols('y',real=True)
    H=s.Function('H')(r)
    for e in (0,1,2):
        for d in (1,2):
            h=r**(-e)*H*s.exp(s.I*eta*(y-M*r**d))
            De=s.diff(h,r)+s.I*eta*M*d*r**(d-1)*h+e*h/r
            checks[f'e{e}d{d}:weighted_drift_pullback_retains_full_phase']=s.simplify(r**e*De*s.exp(-s.I*eta*(y-M*r**d))-s.diff(H,r))==0
    checks['zero_Haar_and_two_Haar_moments_do_not_detect_nonzero_mode_trace']=s.integrate(s.cos(2*s.pi*y),(y,0,1))==0 and jump!=0
    # This last comparison is discontinuous data, never a continuum field.
    return checks,dict(scope='classical patch derivatives and exact jump defects; no distributional terms erased; fixtures are not canonical smooth velocities')

def integration_controls():
    checks={};u,w=s.symbols('u omega',real=True,nonzero=True)
    pieces=[(s.S.One+u+u*u,1,2),(3-u+u**3,2,4)]
    for p in (1,2,3,4):
        boundary=0;bulk=0
        for F,left,right in pieces:
            for r in range(p):
                d=s.diff(F,u,r)
                boundary+=(-1)**r/(s.I*w)**(r+1)*(d.subs(u,right)*s.exp(s.I*w*right)-d.subs(u,left)*s.exp(s.I*w*left))
            bulk+=(-1)**p/(s.I*w)**p*s.integrate(s.diff(F,u,p)*s.exp(s.I*w*u),(u,left,right))
        direct=sum(s.integrate(F*s.exp(s.I*w*u),(u,left,right)) for F,left,right in pieces)
        checks[f'p{p}:all_seam_and_exterior_traces_reconstruct_exact_integral']=s.simplify(direct-boundary-bulk)==0
        checks[f'p{p}:bulk_only_mutation_detects_boundary_terms']=s.simplify(direct-bulk)!=0
    rows=[]
    # n unit boxes convolved, shifted by 1. Degree n-1, C^(n-2).
    for n in (2,3,4):
        for odd in (1,3):
            omega=s.pi*odd
            polys=[]
            for j in range(n):
                F=s.expand(sum((-1)**l*s.binomial(n,l)*(u-1-l)**(n-1) for l in range(j+1))/s.factorial(n-1))
                polys.append(F)
            direct=sum(s.integrate(F*s.exp(s.I*omega*u),(u,j+1,j+2)) for j,F in enumerate(polys))
            convolution=s.exp(s.I*omega)*((-2)/(s.I*omega))**n
            checks[f'n{n}odd{odd}:finite_regular_compact_fixture_has_exact_M_minus_n_leak']=s.simplify(direct-convolution)==0
            checks[f'n{n}odd{odd}:nonzero_leak_has_exact_magnitude']=s.simplify(s.expand_complex(direct*s.conjugate(direct))-(2/omega)**(2*n))==0
            for r in range(n-1):
                traces=[-s.diff(polys[0],u,r).subs(u,1),s.diff(polys[-1],u,r).subs(u,n+1)]
                traces += [s.diff(polys[j-1]-polys[j],u,r).subs(u,j+1) for j in range(1,n)]
                checks[f'n{n}odd{odd}r{r}:all_lower_derivative_traces_match']=all(v==0 for v in traces)
            rows.append(dict(convolutions=n,smoothness=n-2,first_jump=n-1,decay_power=n,half_integer_M=Q(odd,2)))
    return checks,dict(fixed_mode_comparisons=rows,scope='compact finite-regularity comparisons, not canonical C-infinity sources; no asserted lower bound on actual residual')

def moving_controls():
    checks={};u,z,t=s.symbols('u z t',real=True)
    seam=1+z+t
    left=1+u+z;right=2*u+t
    whole=s.integrate(left,(u,0,seam))+s.integrate(right,(u,seam,4))
    for v in (z,t):
        bulk=s.integrate(s.diff(left,v),(u,0,seam))+s.integrate(s.diff(right,v),(u,seam,4))
        trace=(left-right).subs(u,seam)*s.diff(seam,v)
        checks[f'{v}:moving_integral_keeps_seam_speed_times_value_jump']=s.expand(s.diff(whole,v)-bulk-trace)==0
        checks[f'{v}:dropping_moving_seam_term_changes_classical_derivative']=s.expand(s.diff(whole,v)-bulk)!=0
    # Continuous hinge, yet its slow derivative has a value jump.
    left=seam-u;right=u-seam
    for r in (0,1,2):
        Delta=s.diff(left-right,u,r).subs(u,seam)
        slow=s.diff(left-right,z,u,r).subs(u,seam)
        radial=s.diff(left-right,u,r+1).subs(u,seam)
        checks[f'r{r}:mixed_trace_chain_rule_keeps_next_radial_trace']=s.expand(s.diff(Delta,z)-slow-s.diff(seam,z)*radial)==0
    checks['continuous_Phi_alone_does_not_make_Dz_Phi_value_trace_zero']=(left-right).subs(u,seam)==0 and s.diff(left-right,z).subs(u,seam)==2
    return checks,dict(scope='moving trace identities in one fixed physical U chart, before integer-cover transport; not permission to commute moving A1 or T1 through Dz')

def budget_controls():
    checks={};rows=[]
    for p in (1,2,5):
        a=Q(2);k=Q(1,10);target=a+(p-Q(1,2))*k
        zero=seam_orders(a,p,{r:None for r in range(p)},kappa=k,target=target)
        last=seam_orders(a,p,{r:(Q(0) if r==p-1 else None) for r in range(p)},kappa=k,target=target)
        checks[f'p{p}:all_exact_matching_recovers_bulk_power']=zero['guaranteed_power']==a+p*k and zero['strict_target_pass']
        checks[f'p{p}:last_trace_is_allowed_at_bulk_order']=last['guaranteed_power']==a+p*k
        if p>1:
            leak=seam_orders(a,p,{r:(Q(0) if r==0 else None) for r in range(p)},kappa=k,target=target)
            checks[f'p{p}:smallest_unmatched_trace_blocks_strict_target']=not leak['strict_target_pass'] and leak['guaranteed_power']==a+k
            improved=seam_orders(a,p,{r:(Q(1) if r==0 else None) for r in range(p)},kappa=k,target=target)
            checks[f'p{p}:quantified_extra_epsilon_smallness_can_replace_exact_matching']=improved['strict_target_pass']
        rows.append(last)
    a=A1_trace_request(physical_power=0);d=A1_trace_request()
    checks['physical_p30000001_and_original_Phi_raw60000006_retained']=a['integrations']==30000001 and a['original_raw_Phi_order']==60000006
    checks['physical_Dz_Phi_match_through29999999_Phi_through30000000']=a['matching_Dz_weighted_Phi_through']==29999999 and a['sufficient_matching_weighted_Phi_through_for_moving_seam']==30000000
    checks['normalized_direct_A1_has_no_added_nonstationary_trace_request']=d['integrations']==0 and d['boundary_recipe'] is None
    huge=boundary_recipe(30000001)
    checks['giant_trace_request_is_deferred_not_materialized']=not huge['materialized'] and huge['trace_last']==30000000
    c=certificate()
    checks['no_actual_canonical_trace_or_complete_decay_claim']=not c['actual_input_validation_complete_regular_decay_or_blowup_proved'] and not a['actual_closed_cover_trace_identities_verified']
    return checks,dict(certificate=c,small_exact_power_budgets=rows)

def cover_controls():
    checks={};seed={(1,0):Q(2),(0,1):Q(3),(-1,0):Q(2),(0,-1):Q(3)}
    checks['original_J_and_nonunimodular_cover_retained']=geometry.J==((3,1),(1,5)) and geometry.determinant(geometry.J)==14
    for delta in (1,2,3):
        transported=transported_modes(seed,delta)
        checks[f'delta{delta}:matching_physical_trace_uses_transposed_relative_cover']=transported_modes(transported,0)==transported_modes(seed,delta)
        checks[f'delta{delta}:equal_raw_mode_lists_are_not_equal_physical_traces']=transported_modes(seed,0)!=transported_modes(seed,delta)
        checks[f'delta{delta}:zero_Haar_is_preserved_but_does_not_make_traces_match']=(0,0) not in transported and transported!=seed
        checks[f'delta{delta}:constant_mode_Haar_is_unchanged']=transported_modes({(0,0):Q(7)},delta)=={(0,0):Q(7)}
    checks['relative_cover_composition_uses_same_original_matrix']=transported_modes(transported_modes(seed,1),2)==transported_modes(seed,3)
    return checks,dict(scope='finite Fourier trace comparisons on ORIGINAL relative covers; canonical infinite Fourier/mixed traces remain unverified')

@lru_cache(maxsize=1)
def review():
    checks={};details={}
    for name,fn in (('compact',compact_controls),('integration',integration_controls),('moving',moving_controls),('budget',budget_controls),('cover',cover_controls)):
        c,d=fn();checks.update({name+':'+k:v for k,v in c.items()});details[name]=d
    return dict(checks=checks,details=details)
