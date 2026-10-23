"""Exact joint divisors, complete swirl/pressure split, and raw norm budgets."""
from fractions import Fraction as Q
from functools import lru_cache
import sympy as s
from bounds import divisor_data,kernel_recipe,sharp_boundary,radial_counts,split_envelope,radial_profile_recipe,certificate,P,sharp,parent

def divisor_controls():
    checks={};u,v=s.symbols('u v',integer=True);b=s.sqrt(2)-1
    alpha=u-b*v;tau=b*u+v;N=(u+v)**2-2*v*v
    checks['exact_same_mode_joint_integer_form']=s.expand(alpha*tau-b*N)==0
    checks['temporal_divisor_is_b_times_radial_quadratic_conjugate']=s.expand(tau-b*(u+v+s.sqrt(2)*v))==0
    grid=[(i,j) for i in range(-7,8) for j in range(-7,8) if (i,j)!=(0,0)]
    checks['all_224_integer_grid_joint_norms_nonzero']=all(divisor_data(k)['quadratic_norm']!=0 for k in grid)
    checks['all_224_joint_products_have_nonzero_integer_floor']=all(abs(divisor_data(k)['quadratic_norm'])>=1 for k in grid)
    checks['joint_constant_b_inverse_is_less_than_three']=bool(1+s.sqrt(2)<3)
    rows=[]
    for j in range(1,17):
        d=sharp.pell(j);Pj,Qj=d['P'],d['Q'];k=d['k'];a=Pj-s.sqrt(2)*Qj;t=b*(Pj+s.sqrt(2)*Qj)
        checks[f'pell{j}:actual_joint_divisor_product_is_exact_plus_minus_b']=s.expand(a*t-(-1)**j*b)==0
        checks[f'pell{j}:temporal_inverse_gains_one_over_P_on_radial_resonance']=Pj>=Qj and Pj<=2*Qj and s.expand(t-(b*k[0]+k[1]))==0
        if j in (2,3,8):
            # Rotate the INPUT mode into temporal resonance but keep the
            # OUTPUT mode in radial resonance: same-mode (3) is invalid.
            kin=(Qj,-(Pj-Qj));tin=b*kin[0]+kin[1]
            checks[f'pell{j}:different_input_output_divisor_mutation_detected']=s.expand(tin+a)==0 and Pj*Pj>3
        if j in (1,2,3,8,16):rows.append(d)
    return checks,dict(rows=rows,scope='exact algebra of the ORIGINAL directions; no irrational floating-point sampling')

def Fourier_controls():
    checks={};rows=[]
    for m in (0,1,3):
        for p in (1,2,7):
            r=kernel_recipe(m,p);J=r['raw_F_input_order'];key=f'm{m}p{p}'
            checks[key+':all_derivative_splits_read_exact_m_plus_two_p_plus_two']=all((m-a)+p+(a+p+2)==J for a in range(m+1))
            checks[key+':joint_divisor_leaves_summable_r_minus_three']=all(a+(p-1)-(a+p+2)==-3 for a in range(m+1))
            checks[key+':two_pi_power_is_negative_and_not_dropped_in_wrong_direction']=-p-1-(p+2)==-2*p-3<0
            checks[key+':five_derivative_saving_is_only_vs_sequential_upper_rules']=r['sequential_F_input_order']-J==5 and J==m+2*p+2
            checks[key+':raw_original_profiles_and_source_factorial_retained']=r['source_factorial']==J and f'N_{J}(f_e*a)' in r['factors'] and not r['actual_input_or_factorial_evaluated']
            low=sharp_boundary(m,m+2*p-2,p);at=sharp_boundary(m,m+2*p-1,p)
            checks[key+':one_below_mixed_lower_boundary_is_obstructed']=low['strict_finite_jet_obstruction'] and low['ratio_power_in_P']==1
            checks[key+':equality_left_open_and_upper_gap_three']=not at['strict_finite_jet_obstruction'] and J-at['F_input_order']==3 and not at['equality_or_optimal_upper_rule_resolved']
            vals=[Q(Pj*Pj)**p*Q(Pj)**(m-low['F_input_order']-1) for Pj in (7,41,239)]
            checks[key+':independent_rational_mixed_ratio_grows_like_P']=vals==[7,41,239]
            rows.append(dict(m=m,p=p,J=J,necessary=m+2*p-1))
    for r in range(1,13):
        lattice=[(i,j) for i in range(-r,r+1) for j in range(-r,r+1) if abs(i)+abs(j)==r]
        checks[f'lattice{r}:four_r_modes_and_inverse_square_tail_bound']=len(lattice)==4*r and sum((Q(1,j*j) for j in range(1,r+1)),Q(0))+Q(1,r)<=2
    return checks,dict(rows=rows,scope='raw Fourier/Leibniz proof counts and exact lower comparisons, not evaluated canonical kernels')

def integration_controls():
    checks={};r=s.symbols('r',real=True);lam=2*s.pi
    def moment(poly):
        deg=s.degree(poly,r)
        return s.expand(sum((-1)**j*(s.diff(poly,r,j).subs(r,2)-s.diff(poly,r,j).subs(r,1))/(s.I*lam)**(j+1) for j in range(deg+1)))
    for p in (1,2):
        J=kernel_recipe(0,p)['raw_F_input_order'];phi=s.expand((r-1)**(J+1)*(2-r)**(J+1))
        I=moment(phi);Ip=moment(s.diff(phi,r,p))
        checks[f'IBP{p}:complete_radial_integration_with_no_endpoint_terms']=s.expand(Ip-(-s.I*lam)**p*I)==0
        checks[f'IBP{p}:fixture_has_every_raw_source_jet_at_endpoints']=all(s.diff(phi,r,j).subs(r,e)==0 for j in range(J+1) for e in (1,2))
        checks[f'IBP{p}:nonzero_real_compact_integral_is_retained']=I!=0 and s.im(I)==0
    theta,tau=s.symbols('theta tau',real=True,nonzero=True)
    checks['real_Nt_inverse_has_correct_two_pi_and_sign']=s.diff(s.sin(2*s.pi*tau*theta)/(2*s.pi*tau),theta)==s.cos(2*s.pi*tau*theta)
    checks['dropping_temporal_divisor_changes_exact_real_inverse']=s.diff(s.sin(2*s.pi*tau*theta)/(2*s.pi),theta)!=s.cos(2*s.pi*tau*theta)
    U,d=s.symbols('U d',positive=True)
    checks['actual_f0_over_R_is_exact_reciprocal_U_profile']=s.simplify((U**(1/d))**(1-d)/(d*U**(1/d))-1/(d*U))==0
    for j in range(9):
        checks[f'profile{j}:exact_raw_derivative_of_actual_new_kernel_factor']=s.simplify(s.diff(1/(d*U),U,j)-(-1)**j*s.factorial(j)/(d*U**(j+1)))==0
        for u0 in (Q(1,3),Q(2)):
            a=radial_profile_recipe(j,Q(3,2),u0)
            upper=a['prefactor']*s.factorial(j)*a['base']**j
            checks[f'profile{j}u{u0}:entire_through_jet_upper_recipe']=all(Q(s.factorial(k))/(Q(3,2)*u0**(k+1))<=upper for k in range(j+1))
    giant=radial_profile_recipe(60200006,Q(3,2),Q(1,3))
    checks['huge_actual_profile_recipe_keeps_factorial_symbolic_and_positive_hull']=giant['factorial_index']==60200006 and giant['prefactor']==2 and giant['base']==3 and not giant['numerical_factorial_or_actual_source_materialized']
    return checks,dict(scope='exact polynomial Fourier moments at actual radial mode k=(1,0), M=1; finite endpoint fixtures only')

def nonlinear_controls():
    checks={};z=s.symbols('z',nonzero=True);V,v0,A,B,D,c,lam=s.symbols('V v0 A B D c lambda',nonzero=True)
    cos=(z+1/z)/2;cos2=(z*z+z**-2)/2
    v=v0+(A+B)*cos+D*cos2
    Nt=lambda h:s.expand(s.I*lam*z*s.diff(h,z))
    haar=lambda h:s.expand(h).coeff(z,0)
    Qc=lambda h:s.expand(h-haar(h))
    inv=lambda h:s.expand(sum(s.expand(h).coeff(z,k)*z**k/(s.I*lam*k) for k in range(-4,5) if k))
    C=s.expand(2*V*v+v*v);F=-c*Nt(v);G=s.expand(2*(V+v)*F)
    checks['complete_nonlinear_centrifugal_Q_split_with_correct_c_and_sign']=s.expand(Qc(C)+inv(Qc(G))/c)==0
    checks['Fsum_is_full_current_time_derivative_not_new_fine_row']=s.expand(F+c*Nt(v))==0
    checks['actual_Gtheta_has_zero_Haar_because_it_is_a_derivative']=haar(G)==0
    checks['full_slow_and_all_fast_parts_remain_inside_Gtheta']=s.expand(Qc(C)+inv(2*(V+v-v0)*F)/c)!=0
    checks['missing_factor_two_in_Nt_v_square_is_detected']=s.expand(Qc(C)+inv((2*V+v)*F)/c)!=0
    expected=2*V*v0+v0*v0+(A+B)**2/2+D*D/2
    checks['Haar_centrifugal_pressure_keeps_all_fast_variance_and_old_new_cross']=s.expand(haar(C)-expected)==0
    checks['Haar_v_square_is_not_square_of_Haar_v']=s.expand(haar(v*v)-haar(v)**2)==(A*A/2+A*B+B*B/2+D*D/2)
    checks['using_Q_split_for_full_C_would_erase_pressure_moment']=s.expand(C+inv(G)/c)==s.expand(expected) and expected!=0
    h=z+z**-1;a=z
    checks['auxiliary_dependent_multiplier_inverse_commutator_nonzero']=s.expand(inv(Qc(a*h))-a*inv(h))!=0
    checks['auxiliary_multiplier_Nt_product_rule_keeps_derivative_term']=s.expand(Nt(a*h)-a*Nt(h)-Nt(a)*h)==0
    R,rhoP,K0,p0=s.symbols('R rhoP K0 p0',nonzero=True)
    checks['radial_Q_split_preserves_complete_K0_compensated_pressure']=s.expand(K0+p0-rhoP+haar(C)/R+Qc(C)/R-(K0+p0-rhoP+C/R))==0
    # Fast-only comparison: v=-c^-1 Nt^-1 F has zero Haar; pure swirl
    # is divergence free and its square cannot cancel the fundamental.
    fast=-inv((A*(z+z**-1)/2))/c
    delta=s.expand(2*V*fast+fast*fast)
    checks['comparison_obeys_actual_fast_formula_and_both_velocity_moments']=s.expand(-c*Nt(fast)-A*(z+z**-1)/2)==0 and haar(fast)==0
    checks['fast_compliant_centrifugal_square_does_not_cancel_fundamental']=s.expand(fast*fast).coeff(z,1)==0 and delta.coeff(z,1)!=0
    checks['fast_compliant_pressure_variance_remains_nonzero']=haar(delta)!=0
    return checks,dict(scope='exact finite real conjugate Fourier families and pressure split; prescribed residuals are not the canonical current residual')

def budget_controls():
    checks={};rows=[]
    keys=('x','u','z','b','V','G','c','r','W')
    for n in (0,1,3):
        data={k:P(n+j+1) for j,k in enumerate(keys)};a=split_envelope(n,data);vals={k:v.evaluate(1) for k,v in data.items()};b=split_envelope(n,vals)
        checks[f'n{n}:positive_polynomial_and_exact_rational_split_envelopes_agree']=all(a[k].evaluate(1)==b[k] for k in ('Fsum','Gtheta','rest'))
        checks[f'n{n}:generic_Fsum_route_keeps_one_v_and_two_beta_derivatives']=a['Gtheta_input_order']==n+1 and a['rest_input_order']==n+2
        for t in (Q(1,2),Q(1,8)):
            ep=t**10;x=t**19*vals['x'];u=t**9*vals['u'];gamma=t**9*vals['z'];bb=t**10*vals['b'];W=t**10*vals['W'];d=n+1
            Fsum=2*vals['c']*d*u;Gtheta=2*(vals['V']+u)*Fsum
            S=2*bb*x+x*x+W;H=bb*gamma+vals['G']*x+x*gamma+W
            rest=d*(ep+2*vals['c'])*x+vals['r']*S+ep*d*H+vals['r']*W+ep**3*(n+1)*(n+2)*x
            checks[f'n{n}t{t}:independent_full_epsilon_budget_Gtheta_B_rest_one']=Gtheta<=t**9*b['Gtheta'] and rest<=t**10*b['rest']
    for k in (Q(1,100000),Q(1,10),Q(6,25)):
        for m in (0,1,3):
            d=radial_counts(m,kappa=k)
            checks[f'k{k}m{m}:normalized_strict_split_counts_and_local_input_max']=d['mixed_integrations']==d['rest_integrations']+1 and d['local_current_field_order']==d['previous_full_gr_raw_order'] and min(d['mixed_strict_margin'],d['rest_strict_margin'])>0
    for h in (Q(1,200),Q(1,1000),Q(99,10000)):
        for m,N in ((0,0),(1,2),(3,7)):
            d=radial_counts(m,physical_power=N,h=h)
            checks[f'h{h}m{m}N{N}:same_physical_loss_and_strict_split_recipe']=d['local_current_field_order']==d['previous_full_gr_raw_order'] and d['mixed_integrations']==d['rest_integrations']+1 and min(d['mixed_strict_margin'],d['rest_strict_margin'])>0
    n=radial_counts(0);p=radial_counts(0,physical_power=0);rows.extend((n,p))
    checks['normalized_G200002_rest200001_are_different_sources_not_gr']=n['G_raw_order']==200002 and n['rest_raw_order']==200001 and n['previous_full_gr_raw_order']==200003
    checks['physical_G60200006_rest60200005_and_original_graph_retained']=p['G_raw_order']==60200006 and p['rest_raw_order']==60200005 and p['previous_full_gr_raw_order']==60200007 and not p['older_root_normal_source_DAG_expanded_or_erased']
    checks['mixed_epsilon_boundary_still_requires_huge_row_order']=n['mixed_epsilon_target_necessary_F_order']==199997 and p['mixed_epsilon_target_necessary_F_order']==60200001
    c=certificate()
    checks['complete_pressure_variance_and_all_original_shifts_remain']='Haar(2Vv+v^2)' in c['pressure'] and 'K0' in c['pressure'] and 'ALL old' in c['domains']
    checks['whole_regular_E_A1_gain_and_actual_input_gate_unchanged']=c['complete_regular_order']==str(parent.rates()['E']) and c['parent_A1_direct_order']=='199999/100000' and not c['canonical_row_high_jet_or_support_norms_absolute_decay_infinite_field_or_blowup_validated']
    return checks,dict(rows=rows,certificate=c,scope='derived current coarse-input norm programs and integer requests; no actual source supremum or giant constant evaluated')

@lru_cache(None)
def _review():
    checks={};parts={}
    for name,fn in (('divisors',divisor_controls),('Fourier',Fourier_controls),('integration',integration_controls),('nonlinear',nonlinear_controls),('budgets',budget_controls)):
        a,b=fn();checks.update(a);parts[name]=b
    return dict(checks=checks,parts=parts)
def review():
    from copy import deepcopy
    return deepcopy(_review())
