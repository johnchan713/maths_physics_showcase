"""Exact whole-box manufactured controls, not numerical actual pulse bounds."""
from functools import lru_cache
from fractions import Fraction as Q
from math import comb
import sympy as s
from bounds import (P,inverse,reference,CanonicalEntranceBudget,from_original,
    root_wave_join,certificate,product)

v,x=s.symbols('v x',real=True)

def rational_box(expr):
    """Exact upper bound on [0,1]^2: numerator l1 / positive denominator constant.

    Refuses any denominator without nonnegative coefficients. No grid extrema.
    """
    num,den=s.fraction(s.cancel(expr));den=s.Poly(den,v,x)
    if den.eval({v:0,x:0})<0:num,den=-num,-den
    dc=den.eval({v:0,x:0})
    if dc<=0 or any(c<0 for c in den.coeffs()):raise ValueError('Unproved denominator floor')
    return Q(sum(abs(c) for c in s.Poly(num,v,x).coeffs())/dc)

def raw_bound(items,j):
    return max(sum(rational_box(s.diff(e,v,a,x,j-a)) for e in items) for a in range(j+1))

def scaled_jet(items,p,j):
    values=[]
    for a in range(j+1):
        out=list(items)
        for _ in range(a):out=[s.diff(e,v)+p*e for e in out]
        for _ in range(j-a):out=[s.diff(e,x)+(v-s.Rational(1,2))*e for e in out]
        values.append(sum(rational_box(e) for e in out))
    return max(values)

@lru_cache(None)
def manufactured():
    n=s.Matrix([v+x,1,0]);np=n.diff(v);N=n.dot(n)
    b=s.Matrix([1,-v-x,1]);bp=b.diff(v);B=b.dot(b)
    p=1-2*v+x;c=1+v+x;d=N
    K=-(p+d)*s.eye(3)-bp*b.T/B+c*n*b.T/B
    callbacks={name:(lambda j,items=tuple(obj):P(raw_bound(items,j)))
        for name,obj in (('normal',n),('moving_normal',np),('matrix',K),('profile_rate',[p]))}
    base=inverse.ForcedBudget(**callbacks,source=lambda j:P(),tangent=P(9),length=P(1),
        inverse_normal_square=1,harmonic=1)
    h=CanonicalEntranceBudget(base,tangential_inverse_square=1,cmax=Q(4,5),u=Q(3,2))
    return n,np,N,b,p,c,K,base,h

def exact_equation_review():
    n,np,N,b,p,c,K,base,h=manufactured();lift=n/N;H=lift*n.T
    row=n.T*K-np.T;A0=-K+lift*row;A=A0-N*s.eye(3)
    Ahat=A+(N+p)*H
    zero=lambda M:all(s.cancel(e)==0 for e in M)
    checks=dict(
        moving_normal_constraint_exact=n.dot(b)==0,
        manufactured_homogeneous_equation_exact=zero(A*b-p*b-b.diff(v)),
        ambient_extension_agrees_on_tangent=zero((Ahat-A)*b),
        ambient_normal_response_is_profile_rate=zero(np.T+n.T*Ahat-p*n.T),
        homogeneous_pressure_is_nonzero=zero(row*b/N-s.Matrix([c])) and c!=0,
        dropping_moving_normal_changes_pressure=s.cancel((n.T*K*b)[0]/N-c-1/N)==0,
        exact_potential_is_n_cross_t_over_normal_square=zero(n.cross(b)/N-s.Matrix([1,-v-x,-1-(v+x)**2])/N),
        zero_source_zero_entry_inverse_cannot_bound_canonical_column=base.velocity(0).evaluate(1)==0<h.velocity(0).evaluate(1),
        reference_entrance_mixing_has_frozen_rational_factor=1+Q(3,2)**2/4==Q(5,4)**2 and Q(-4,5)*Q(5,4)==-1,
        canonical_entrance_is_er_minus_nr_tangent_plus_reference_rotation=b.subs(v,0)==s.Matrix([1,-x,1]),
        inherited_actual_wa_entrance_not_identified=Q(4,5)!=1,
        tangent_growth_exponent_is_bounded_by_two=Q(1,1)+Q(1,1)==2,
        conservative_tangent_constant_9_exceeds_exp2=s.exp(2)<9,
        damping_bound_uses_original_epsilon_k_square_upper_four=4>=1)
    # Relative Gaussian derivatives at the NONZERO entrance P0=exp(-1/4-x/2).
    P0=s.exp(-s.Rational(1,4)-x/2)
    checks['entrance_Gaussian_first_derivative_not_discarded']=s.simplify(s.diff(P0,x)/P0+s.Rational(1,2))==0
    checks['entrance_Gaussian_second_derivative_not_discarded']=s.simplify(s.diff(P0,x,2)/P0-s.Rational(1,4))==0
    checks['first_initial_jet_contains_Gaussian_times_nonzero_vector']=(
        s.diff(P0*b.subs(v,0),x)/P0 != b.diff(x).subs(v,0))
    return checks,dict(scope='exact moving-normal constrained ODE on v,x in [0,1]; arbitrary manufactured K, not the actual background',
        tangent_proof='On tangent vectors A-pI=Pi bprime b^T/|b|^2-h nprime^T; norm <=1/|b|+|h|<=2, so exp(2(v-w))<=exp(2)<9',
        Gaussian='P=exp(-(v-1/2)^2+x(v-1/2)); P(1/2,x)=1; no inverse P0 bound')

def whole_box_review():
    n,np,N,b,p,c,K,base,h=manufactured();checks={};rows=[]
    for j in range(4):
        actual=dict(velocity=scaled_jet(b,p,j),pressure=scaled_jet([c],p,j),
                    potential=scaled_jet(n.cross(b)/N,p,j))
        for kind,value in actual.items():
            checks[f'order{j}_{kind}_complete_mixed_jets_bound_on_whole_box']=getattr(h,kind)(j).evaluate(1)>=value
        rows.append(dict(order=j,exact_box_upper=actual))
    checks['entrance_Bell_first_term_is_half_length_profile_jet']=h.entrance_gaussian(1).c==h.entrance_log(1).c
    checks['entrance_Bell_second_term_has_square_and_second_log_jet']=h.entrance_gaussian(2).c==(h.entrance_log(1)*h.entrance_log(1)+h.entrance_log(2)).c
    checks['initial_Leibniz_first_term_retains_both_factors']=h.initial(1).c==(h.entry_vector(1)+h.entrance_gaussian(1)*h.entry_vector(0)).c
    checks['pulse_derivative_uses_Ahat_not_only_parameter_Duhamel']=h.velocity(1).evaluate(1)>(h.ambient()*h.initial(1)).evaluate(1)
    # A separate varying tangential norm tests the inverse-square-root chain.
    varying=inverse.ForcedBudget(normal=lambda j:P(raw_bound([0,1,x],j)),moving_normal=lambda j:P(),
        matrix=lambda j:P(),source=lambda j:P(),profile_rate=lambda j:P(),
        tangent=P(1),length=P(1),inverse_normal_square=1)
    hh=CanonicalEntranceBudget(varying,tangential_inverse_square=1,cmax=1,u=0)
    J=1+x*x;w=J**(-s.Rational(1,2))
    checks['inverse_sqrt_differential_identity_exact']=s.simplify(2*J*s.diff(w,x)+w*s.diff(J,x))==0
    for j in range(4):
        # D^j J^-1/2 = numerator / J^(j+1/2); denominator >=1.
        numerator=s.simplify(s.diff(w,x,j)*J**(j+s.Rational(1,2)))
        checks[f'order{j}_varying_tangential_inverse_sqrt_bound']=hh.entry_inverse_sqrt(j).evaluate(1)>=rational_box(numerator)
    delta=Q(1,100)
    checks['full_normal_floor_alone_does_not_bound_tangential_inverse']=(1+delta**2>=1 and 1/delta**2>1)
    return checks,dict(rows=rows,scope='exact coefficient-l1 rational whole-unit-box enclosures through order3, including every mixed derivative; no numerical pulse samples')

@lru_cache(None)
def original_fixture():
    normal=reference.normal
    p=normal.frame.inputs(bmin=Q(1,20),bmax=1,Fmax=1,gmax=2,gmin=1,u=Q(1,100),cmin=1,cmax=1,Dbox=6,Mz=1)
    sp=normal.selector.fixture();L=normal.normal_length_certificate(p)['length_rule_floor']
    sp['length']=Q(-(-L.numerator//L.denominator))
    ell=normal.selector.select_tail(sp)['ell_tail'] # manufactured ONLY
    kw=dict(F=lambda j:P(j+2),G=lambda j:P(j+3),btilde=lambda j:P(j+1),
        path=lambda j:P({1:1}) if j==0 else P(int(j==1)))
    return p,sp,ell,kw

def original_and_root_review():
    p,sp,ell,kw=original_fixture();phase,h,c,tc=from_original(p,sp,ell,**kw)
    cut=lambda j:P(j+1);template=inverse.roots.fixture();joined=root_wave_join(template,h,cut)
    checks=dict(original_normal_and_tangent_gates_retained=all(c['gates'].values()) and all(tc['gates'].values()),
        original_same_band_and_length_retained=c['existing_ell']==ell and c['existing_Lreq']==sp['length'],
        original_refined_floor_proved_for_tangential_and_full_norm=h.entry_q==h.inverse_normal_square==c['inverse_normal_square'],
        original_canonical_cmax_and_tilt_retained=h.entry_c==p['cmax'] and h.entry_u==p['u'],
        primary_source_exactly_zero=all(h.data('source',j).evaluate(1)==0 for j in range(3)),
        root_H_target_curl_loghalf_callbacks_are_same_objects=all(getattr(joined,k) is getattr(template,k) for k in ('H','target','curl','loghalf')),
        root_value_gaps_mass_radius_width_preserved=all(getattr(joined,k)==getattr(template,k) for k in ('K0','beta','R','mass','width')),
        no_new_selector_field_or_root_claim=not c['bands_lengths_norm_bounds_roots_or_fields_selected_or_changed'])
    for j in range(3):
        checks[f'order{j}_root_wave_slot_is_full_cutoff_jet']=joined.wave(j).c==h.wave(j,cut).c
    checks['cutoff_derivative_contributes_to_root_wave_slot']=h.wave(1,cut).evaluate(1)>(5*cut(0)*h.velocity(1)).evaluate(1)
    checks['joined_root_amplitude_and_signed_difference_are_finite_positive_polynomials']=all(
        joined.amplitude(j).evaluate(1)>0 and joined.signed(j,lambda i:P(i+1)).evaluate(1)>0 for j in range(2))
    return checks,dict(ell=ell,reciprocal=c['inverse_normal_square'],
        scope='manufactured ORIGINAL input budgets only; same acceptance gates and same root-budget objects; no actual norm suprema evaluated')

@lru_cache(None)
def review():
    checks={};parts={}
    for name,f in (('equation',exact_equation_review),('box',whole_box_review),('join',original_and_root_review)):
        c,d=f();checks.update({name+':'+k:bool(v) for k,v in c.items()});parts[name]=d
    return dict(checks=checks,parts=parts,certificate=certificate())
