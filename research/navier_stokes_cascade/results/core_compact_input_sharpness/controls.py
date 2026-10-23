"""Exact resonance, raw-jet and inherited request controls; no pulse solve."""
from fractions import Fraction as Q
import sympy as s
from bounds import pell,select_mode,finite_jet_obstruction,epsilon_target,compare_request,actual_frontier,parent

def pell_controls():
    checks={};rows=[]
    for j in range(1,25):
        d=pell(j);P,Qv=d['P'],d['Q'];omega=P-s.sqrt(2)*Qv
        checks[f'pell{j}:nonzero_integer_norm_and_actual_vr_divisor']=d['norm']==(-1)**j and s.expand(omega*(P+s.sqrt(2)*Qv))==d['norm']
        checks[f'pell{j}:real_pair_nonzero_k_and_l1_size']=d['k']!=(0,0) and sum(abs(v) for v in d['k'])==P and P<=2*Qv and Qv<=P
        if j in (1,2,3,8,24):rows.append(d)
    for M in (Q(1),Q(2),Q(9,4),Q(10),Q(1000),Q(10)**12,Q(10)**40):
        d=select_mode(M);P,Qv=d['P'],d['Q']
        checks[f'M{M}:same_mode_selection_without_retuning_M']=P*P>=M and P*P<9*M and P<=2*Qv<=2*P and abs(d['norm'])==1
        # |omega|<=1/P, hence M|omega|<=M/P<=P; no float sqrt.
        checks[f'M{M}:full_radial_phase_frequency_at_most_P']=M<=P*P
    return checks,dict(pell_rows=rows,scope='exact algebra in Q(sqrt2), any original M>=1, no carrier changed')

def phase_controls():
    checks={}
    U,a,y1,y2,M=s.symbols('U a y1 y2 M',real=True)
    for j in (2,3,5):
        d=pell(j);P,Qv=d['P'],d['Q'];k1,k2=d['k'];b=s.sqrt(2)-1;omega=P-s.sqrt(2)*Qv
        phase=k1*y1+k2*y2-M*omega*U
        shifted=k1*(y1+M*a)+k2*(y2-M*a*b)-M*omega*(U+a)
        checks[f'phase{j}:complete_shift_cancels_radial_source_oscillation']=s.expand(shifted-phase)==0
        wrong=k1*y1+k2*y2-M*omega*(U+a)
        checks[f'phase{j}:unshifted_primitive_does_not_give_same_integral']=s.expand(wrong-phase)!=0
        reversed=k1*(y1-M*a)+k2*(y2+M*a*b)-M*omega*(U+a)
        checks[f'phase{j}:wrong_drift_sign_is_detected']=s.expand(reversed-phase+2*M*omega*a)==0 and s.expand(reversed-phase)!=0
        # Integrate the two conjugate Fourier modes before taking their
        # real part; do not rely on symbolic trigonometric simplification
        # of a large irrational but y2-independent phase offset.
        checks[f'phase{j}:real_cosine_pair_has_zero_Haar_and_moment']=k2>0 and all(s.integrate(s.exp(sign*2*s.pi*s.I*k2*y2),(y2,0,1))==0 for sign in (-1,1))
    r=s.symbols('r',positive=True)
    # Independent exact finite-regularity bump fixture on [1,2]. Smooth
    # flat bumps are used by the proof, not replaced by this diagnostic.
    for J in (0,1,2,4):
        phi=(r-1)**(J+2)*(2-r)**(J+2)
        mass=s.integrate(phi,(r,1,2));phi=s.expand(phi/mass)
        checks[f'bump{J}:unit_mass_and_compact_finite_jet_extension']=s.integrate(phi,(r,1,2))==1 and all(s.diff(phi,r,z).subs(r,edge)==0 for z in range(J+2) for edge in (1,2))
        # Actual specialization e=1,d_r=1 has f_e(U)=U, k_e(U)=1/U.
        checks[f'bump{J}:actual_kernel_f_e_inverse_factor_retained']=s.integrate(s.cancel(r*(phi/r)),(r,1,2))==1
        checks[f'bump{J}:dropping_actual_f_e_inverse_changes_kernel_mass']=s.integrate(r*phi,(r,1,2))==s.Rational(3,2)
        # A dropped shift already changes the M derivative at M=0:
        # integral (r-U)phi = 3/2-U, nonzero at output U*=3.
        checks[f'bump{J}:omitting_shift_changes_integral_not_just_formula']=s.integrate((r-3)*phi,(r,1,2))==-s.Rational(3,2)
        # A zero Haar source does not force A_e itself to be zero.
        checks[f'bump{J}:zero_Haar_compensation_leaves_nonzero_compact_output']=s.integrate(phi,(r,1,2))!=0
    return checks,dict(exact_output='A*k_e*chi_U*B*cos(2pi(k.y-M omega U))',
        continuum_bump='arbitrary fixed nonnegative C-infinity unit-mass phi on positive input hull',
        fixture_scope='polynomial bumps have finite endpoint regularity only; no actual cutoff or current field instantiated')

def norm_controls():
    checks={};rows=[]
    # Exact product-binomial proof of 16^J with 2pi bounded by 8.
    for J in range(9):
        checks[f'raw{J}:full_Leibniz_cost_dominates_all_lower_jets']=all(sum(s.binomial(r,z)*8**(r-z) for z in range(r+1))<=16**J for r in range(J+1))
        # With amplitude P^-J, raw order r is P^(r-J). The all-jet
        # class is not certified by the uniformly controlled first J jets.
        checks[f'raw{J}:finite_jet_fixture_not_uniform_all_jet_source']=(J+1-J)>0
    for m in (0,1,3):
        for p in (1,2,7):
            boundary=m+2*p
            below=finite_jet_obstruction(m,boundary-1,p)
            at=finite_jet_obstruction(m,boundary,p)
            upper=finite_jet_obstruction(m,boundary+3,p)
            key=f'm{m}p{p}'
            checks[key+':one_fewer_than_two_p_is_obstructed']=below['uniform_estimate_obstructed'] and below['ratio_growth_in_P']==1
            checks[key+':equality_not_disproved_or_promoted_to_upper_bound']=not at['uniform_estimate_obstructed'] and at['ratio_growth_in_P']==0
            checks[key+':inherited_upper_rule_is_three_derivatives_above_boundary']=upper['source_jet']-at['source_jet']==3 and upper['sufficient_source_order']==upper['source_jet']
            # Independent rational lower ratio M^p P^(m-J), using M=P^2.
            J=boundary-1
            vals=[Q(P*P)**p*Q(P)**(m-J) for P in (7,41,239)]
            checks[key+':Pell_subsequence_ratio_grows_exactly_like_P']=vals==[7,41,239]
            rows.append(dict(m=m,p=p,lower_J=boundary,upper_J=boundary+3))
    frac=finite_jet_obstruction(2,4,Q(3,2))
    checks['noninteger_gain_also_has_sharp_power_boundary']=frac['necessary_source_order']==5 and frac['uniform_estimate_obstructed']
    return checks,dict(rows=rows,raw_upper='abs(A)*16^J*A_J*P^J',raw_lower='abs(A)*c_star*P^m',
        all_fixed_profiles_retained='A_J=N_J(f_e^-1 phi B); c_star=abs(k_e(U*)chi_U(U*))>0',
        norms_or_huge_constants_materialized=False)

def route_controls():
    checks={};table=[]
    expected={(False,0):(199998,200003),(False,1):(2,7),
              (True,0):(60200002,60200007),(True,1):(60000006,60000011)}
    for phys in (False,True):
        for e in (0,1):
            r=compare_request(e,0,physical_power=0 if phys else None);table.append(r)
            lower,upper=expected[(phys,e)]
            key=f'{"physical" if phys else "normalized"}A{e}'
            checks[key+':exact_parent_source_jet_and_sharp_route_threshold']=r['necessary_finite_jet_source_order']==lower and r['retained_raw_source_order']==upper
            checks[key+':strict_rule_is_five_above_integer_epsilon_boundary']=r['excess_over_necessary_order']==5 and r['strict_epsilon_margin']>0
            a,target=r['source_order'],r['normalized_target']
            below=epsilon_target(0,lower-1,a,target);at=epsilon_target(0,lower,a,target)
            checks[key+':one_jet_below_boundary_has_strict_epsilon_deficit']=below['finite_jet_route_obstructed'] and below['epsilon_gap']>0
            checks[key+':equality_is_not_resolved_by_comparison']=not at['finite_jet_route_obstructed'] and not at['equality_case_decided']
    for h in (Q(1,200),Q(1,1000),Q(99,10000)):
        for m,N in ((0,0),(1,2),(3,7)):
            for e in (0,1):
                r=compare_request(e,m,physical_power=N,h=h)
                d=parent.compact_physical_orders(e,m,N,h)
                checks[f'h{h}m{m}N{N}e{e}:same_physical_route_not_new_loss_lower_bound']=r['retained_raw_source_order']==d['raw_source_order'] and 3<=r['excess_over_necessary_order']<=5
    # Independent fractional target gains exercise all three possible gaps.
    for x,gap in ((Q(2),5),(Q(9,4),4),(Q(11,4),3)):
        p=x.numerator//x.denominator+1;J=2*p+3
        d=epsilon_target(0,J,0,x,1)
        checks[f'fractional_gain{x}:strict_rule_excess{gap}']=J-d['necessary_source_order']==gap
    frontier=actual_frontier()
    checks['actual_full_current_pressure_source_inequalities_not_evaluated']=not frontier['actual_current_source_high_jet_norms_evaluated'] and 'all six' in frontier['source_identification'] and 'full forced' in frontier['source_identification']
    checks['scalar_obstruction_not_canonical_or_blowup_claim']=not frontier['actual_input_validity_residual_decay_or_blowup_certified'] and 'not in' in frontier['all_jet_class_warning']
    return checks,dict(table=table,actual_input_frontier=frontier,
        remaining_regular_exponent=parent.regular_physical(0,Q(1,200)),
        physical_loss_is_retained_sufficient_route_not_proved_necessary=True)

def review():
    checks={};parts={}
    for name,fn in (('pell',pell_controls),('phase',phase_controls),('norm',norm_controls),('routes',route_controls)):
        c,data=fn();checks.update(c);parts[name]=data
    return dict(checks=checks,parts=parts)
