"""Exact chart/Euler/cutoff identities and independent floor controls."""
from functools import lru_cache
from fractions import Fraction as Q
import sympy as s
from bounds import (directional_envelope,leading_direction,comparison,normal_gate,
                    uniform_slow_bound,certificate)

def chart_controls():
    checks={};v,R,Z,eta,q,X=s.symbols('s R Z eta q X',positive=True)
    Qlabel=s.Rational(1,4)
    for h in (s.Rational(1,200),s.Rational(1,1000),s.Rational(99,10000)):
        A=s.Rational(1,2)+h;D=s.Rational(1,2)-h
        L=1-2*h*Z**2*v**(2*h-1)
        sv=2*Z*v**(2*h)/L
        zeta=Z*v**(-D);xx=R**2/(2*v)
        checks[f'h{h}:implicit_s_Z_keeps_full_denominator']=s.simplify((1-2*h*Z**2*v**(2*h-1))*sv-2*Z*v**(2*h))==0
        checks[f'h{h}:eta_Z_has_exact_one_minus_eta_square']=s.simplify(v**(-D)-D*Z*v**(-D-1)*sv-v**(-D)*(1-zeta*zeta)/L)==0
        B=q*q*(1+X*eta)+X*X+3*eta**2+q*(1-q)*(1+X+eta)
        sub={q:Qlabel*v,X:xx,eta:zeta}
        profile=B.subs(sub,simultaneous=True)
        native=2*eta*(q*s.diff(B,q)-A*B-X*s.diff(B,X))+(1-eta*eta)*s.diff(B,eta)
        for divisor,name in ((R,'F'),(s.S.One,'G')):
            field=v**(-A)*profile/divisor
            actual=s.diff(field,Z)+s.diff(field,v)*sv
            expected=v**(-1)*native.subs(sub,simultaneous=True)/(divisor*L)
            checks[f'h{h}{name}:direct_implicit_chain_matches_native_Euler_direction']=s.simplify(s.powsimp(actual-expected,force=True))==0
        omitted=2*eta*(-A*B-X*s.diff(B,X))+(1-eta*eta)*s.diff(B,eta)
        checks[f'h{h}:freezing_q_in_Borel_coefficient_loses_Euler_transport']=s.expand(native-omitted-2*eta*q*s.diff(B,q))==0
        checks[f'h{h}:actual_q_dependent_cutoff_Euler_term_is_nonzero']=s.expand(native-omitted)!=0
    return checks,dict(scope='exact frozen-label chart differentiation including actual Euler q dependence; polynomial cutoff factors are diagnostic, not the selected field')

def heat_controls():
    checks={};X,eta,h=s.symbols('X eta h',positive=True)
    A=s.Rational(1,2)+h;H=s.Function('H');z=2*(1-eta**2)/X
    E=(2*X)**(-A)*H(z)
    direction=2*eta*(-A*E-X*s.diff(E,X))+(1-eta**2)*s.diff(E,eta)
    checks['actual_heat_exterior_direction_cancels_for_full_heat_function']=s.simplify(direction)==0
    checks['individual_heat_profile_first_jets_are_not_zero']=s.diff(E,X)!=0 and s.diff(E,eta)!=0
    checks['heat_cancellation_needs_both_native_radial_and_eta_terms']=s.simplify(2*eta*(-A*E-X*s.diff(E,X)))!=0
    return checks,dict(scope='exact retained heat exterior E0=C(2X)^-A H(2(1-eta^2)/X), U0=0; no pressure reset or interior heat substitution')

def directional_controls():
    checks={};t,a,b=s.symbols('t a b',real=True)
    poly=2*a*t+b*(1-t*t)
    checks['interior_vertex_bound_has_exact_square_remainder']=s.expand(b+a*a/b-poly-b*(t-a/b)**2)==0
    checks['endpoint_bound_has_exact_positive_factorization']=s.expand(2*a-poly-(1-t)*(2*a-b*(1+t)))==0
    for aa,bb in ((0,0),(1,0),(0,2),(1,2),(2,1),(3,3)):
        upper=directional_envelope(aa,bb)
        point=Q(1) if aa>=bb else Q(aa,bb)
        checks[f'a{aa}b{bb}:sharp_eta_envelope_attained']=upper==2*aa*point+bb*(1-point*point)
        checks[f'a{aa}b{bb}:sharp_eta_bound_no_larger_than_sum']=upper<=2*aa+bb
    return checks,dict(scope='sharp norm envelope for independent native component enclosures, not a claim their extrema occur together in the actual profile')

def fixture():
    return comparison(Q(1,200),dict(value=1,radial=0,angular=0),dict(value=0,radial=0,angular=0),
        Rmin=1,Rmax=2,tilt_upper=0,lambda_max=1,gmin=1,Llow=1,Lhigh=2,Dbox=6)

def floor_controls():
    checks={};rows=[];r=fixture();c=r['coefficients']
    for S in (Q(10),Q(100),Q(1000)):
        e=1/(100*S**4);ki=1/(10*S*S)
        out=normal_gate(r,Q(2,5),S,e,ki);rows.append(out)
        checks[f'S{S}:independent_full_error_is_bounded_with_distinct_powers']=out['gates']['original_slow_gate'] and out['raw_error']<=out['retained_power_error']
        checks[f'S{S}:retained_power_error_no_larger_than_same_collapsed_program']=out['retained_power_error']<=out['collapsed_same_program_error']
        # Independent full q/Euler diagonal-to-direction arithmetic.
        Xhi=Q(8);A=Q(101,200);K=Q(400,99)
        DE=2*A;HZ0=K*3*DE;HZ1=K*4*4*(1+A+Xhi)
        direct=ki+Q(12)/S**3+2*S*e*(HZ0+e*HZ1)
        checks[f'S{S}:native_profile_and_full_diagonal_formula_agree']=direct==out['raw_error']
        checks[f'S{S}:existing_scale_and_data_not_promoted_or_changed']=not out['new_scale_or_length_selected'] and not out['actual_profile_norms_band_or_field_certified']
    out=rows[-1]
    checks['retained_powers_can_pass_when_collapsed_same_program_fails']=out['gates']['tangential_margin'] and 4*out['collapsed_same_program_error']**2>out['beta']
    for S in (Q(1),Q(2),Q(10)):
        envelope=c['A']/S+c['B']/S**2+c['C']/S**3
        later=c['A']/(2*S)+c['B']/(2*S)**2+c['C']/(2*S)**3
        checks[f'S{S}:same_bounds_margin_improves_on_later_bands']=later<envelope
    checks['native_norm_domain_is_exact_original_enlarged_chart']=r['native_Xlo']==Q(1,8) and r['native_Xhi']==8
    return checks,dict(manufactured_norm_budget_rows=rows,scope='exact positive norm data; not actual profile/gap enclosures or accepted physical bands; no comparison of unrelated input norms')

def slow_controls():
    checks={};rows=[]
    for h in (Q(1,200),Q(1,1000),Q(99,10000)):
        for multiplier in (1,2,7):
            ell=multiplier*512/h**2;r=uniform_slow_bound(h,ell);rows.append(r)
            checks[f'h{h}m{multiplier}:original_lower_band_implies_whole_slow_gate']=r['original_band_gate'] and r['slow_upper']<=r['universal_upper_if_original_band_gate']<Q(1,10**24)
            checks[f'h{h}m{multiplier}:rounded_carrier_not_replaced_and_no_new_band']=not r['new_band_selected'] and not r['actual_supplied_band_identified']
    low=uniform_slow_bound(Q(1,200),16)
    checks['missing_original_band_gate_is_reported_not_repaired']=not low['original_band_gate'] and not low['uniform_gate_from_original_bound'] and not low['new_band_selected']
    # All constants in the exponential-series proof are fixed finite integers.
    x=s.symbols('x',positive=True)
    series=s.Add(*(x**j/s.factorial(j) for j in range(17)))
    checks['sixteenth_exponential_series_term_has_positive_remaining_coefficients']=all(v>0 for v in s.Poly(series-x**16/s.factorial(16),x).all_coeffs())
    c=certificate()
    checks['no_actual_input_full_residual_or_blowup_claim']=not c['actual_leading_input_diagonal_same_band_full_residual_or_blowup_certified']
    return checks,dict(rows=rows,certificate=c,scope='uniform analytical exponential majorant on original 0<h<1/100, ell>=512/h^2; examples are exact rational controls')

@lru_cache(maxsize=1)
def review():
    checks={};details={}
    for name,fn in (('chart',chart_controls),('heat',heat_controls),('directional',directional_controls),('floor',floor_controls),('slow',slow_controls)):
        c,d=fn();checks.update({name+':'+k:bool(v) for k,v in c.items()});details[name]=d
    return dict(checks=checks,details=details)
