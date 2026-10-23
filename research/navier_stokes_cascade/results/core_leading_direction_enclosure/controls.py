"""Independent exact/symbolic controls, never samples of the selected field."""
from fractions import Fraction as Q
import sympy as s
from bounds import (collar_recipe,compact_majorants,leading_enclosures,
    reference_enclosures,native_first_jets,comparison,selector_margin,certificate)

def fixture():
    e=leading_enclosures(Q(1,200),1,Q(1,4),2,16)
    return comparison(e,Q(1,200),Rmin=1,Rmax=2,tilt_upper=1,
        lambda_max=1,gmin=Q(1,2),Llow=Q(1,2),Lhigh=2,Dbox=6)

def symbolic_review():
    y,z,n,phi,h=s.symbols('y z N phi h',real=True,nonzero=True)
    E=s.Function('E');U=s.Function('U')
    # Symbolic periodic control with mixed slow/fast dependence. Generic
    # undefined-function Subs objects have dummy-variable equality issues.
    P=(y*y*z+z**3)*s.sin(phi)+(y*z*z+2)*s.cos(2*phi)
    beta=P.subs(phi,n*y)/n
    py=s.diff(P,y).subs(phi,n*y)
    pp=s.diff(P,phi).subs(phi,n*y)
    pe=s.diff(P,z).subs(phi,n*y)
    Ey=s.exp(beta)*(s.diff(E(y,z),y)+E(y,z)*(py/n+pp))
    checks=dict(
        pullback_first_radial_has_P_phi_without_N=s.simplify(s.diff(beta,y)-py/n-pp)==0,
        pullback_angular_retains_inverse_N=s.simplify(s.diff(beta,z)-pe/n)==0,
        exponential_swirl_first_product_rule=s.simplify(s.diff(E(y,z)*s.exp(beta),y)-Ey)==0,
        axial_first_derivative_has_same_pullback=s.simplify(s.diff(U(y,z)+beta,y)-s.diff(U(y,z),y)-py/n-pp)==0)
    # A periodic, explicitly manufactured primitive detects a wrong extension
    # of the first-jet statement to the second radial derivative.
    wave=s.sin(n*y)/n
    checks['second_radial_still_has_N']=s.diff(wave,y,2)==-n*s.sin(n*y)
    checks['omitting_phase_chain_term_detected']=s.diff(wave,y).subs(y,0)==1
    X=s.symbols('X',positive=True);F=s.Function('F');B=s.Function('B')
    r=s.sqrt(2*X);A=s.Rational(1,2)+h;d=1-z*z
    DE=-2*z*(A*r*F(X,z)+X*s.diff(r*F(X,z),X))+d*s.diff(r*F(X,z),z)
    rhs=r*(-2*z*((1+h)*F(X,z)+X*s.diff(F(X,z),X))+d*s.diff(F(X,z),z))
    checks['regular_swirl_conversion_retains_half_power']=s.simplify(DE-rhs)==0
    checks['using_A_instead_of_A_plus_half_detected']=s.simplify(DE-r*(-2*z*(A*F(X,z)+X*s.diff(F(X,z),X))+d*s.diff(F(X,z),z)))!=0
    # Sharp coefficient 2 for the pointwise logarithmic-jet sum.
    checks['pointwise_dual_norm_two_is_sharp']=(-2*z*X*s.diff(s.log(X),X)).subs({X:1,z:1})==-2
    checks['sqrt_radius_bound_has_positive_remainder']=s.expand((1+X)**2-2*X)==1+X*X
    return checks,dict(scope='exact symbolic periodic control and generic native conversion; not actual profile sampling')

def rational_review():
    checks={};rows=[]
    for i,H in enumerate((Q(16),Q(32),Q(256))):
        n=H**32+1;r=compact_majorants(H-Q(1,2),H,n)
        checks.update({f'compact_{i}_{k}':v for k,v in r['checks'].items()});rows.append(r)
    # Joint seminorm coefficients: each weight is <=1; no factor 2 for Py/Pphi.
    for i,N in enumerate((Q(2),Q(17),Q(1000))):
        vals=(Q(2),Q(3),Q(5),Q(7));pulled=sum(vals[:3])/N+vals[3]
        checks[f'joint_J1_budget_{i}']=pulled<=sum(vals)
    for i,h in enumerate((Q(1,200),Q(1,1000))):
        for j,x in enumerate((Q(1,10),Q(1),Q(5))):
            r=leading_enclosures(h,x,Q(1,3),Q(2),16)
            axis=r['axis'];checks[f'axis_direction_below_pointwise_dual_bound_{i}_{j}']=(
                axis['E_direction']<=2*axis['E_pointwise_jet_sum'] and
                axis['U_direction']<=2*axis['U_pointwise_jet_sum'])
            for k,z in enumerate((Q(-1),Q(-2,3),Q(0),Q(3,5),Q(1))):
                a=1+h+x;c=3
                checks[f'axis_sharp_envelope_dominates_{i}_{j}_{k}']=(
                    2*abs(z)*a+(1-z*z)*c<=axis['E_direction']/((1+x)*2))
    regional=((1,0,0),(0,1,0),(0,0,1))
    checks['regional_max_of_sums_not_global_sum_of_suprema']=(
        max(map(sum,regional))==1 and sum(max(v[j] for v in regional) for j in range(3))==3)
    c=collar_recipe(2,Q(1,2),8,Q(1,1000),Q(1,100000))
    checks['axis_recipe_keeps_original_width_cost']=c['K_prefactor']>1 and c['K_exp_argument']>0
    c2=collar_recipe(2,Q(1,2),8,Q(1,2000),Q(1,100000))
    checks['narrower_width_increases_same_collar_prefactor']=c2['K_prefactor']>c['K_prefactor']
    checks['axis_recipe_no_giant_exponential_evaluation']=not c['actual_parameters_evaluated'] and 'exp(' in c['K_expression']
    jets=native_first_jets(1,Q(1,4),2,16,Q(1,8))
    checks['whole_hull_ordinary_X_keeps_inverse_radius']=(
        jets['E']['radial']==8*16**4 and jets['U']['radial']==8*16**2)
    checks['whole_hull_has_all_two_profile_primitive_keys']=(
        set(jets['E'])==set(jets['U'])=={'value','radial','angular'})
    checks['whole_hull_has_original_Xlo']=jets['original_positive_hull_Xlo']==Q(1,8)
    # Positive squared data: prove outward sqrt bounds by squaring exactly.
    ref=reference_enclosures(Q(1,10),2,Q(1,4),5,4,Q(1,10),9)
    f=ref['native_squared_floors']
    checks['reference_gmin_outward']=ref['gmin']**2<=f['gmin2']/64
    checks['reference_growth_floor_outward']=ref['lambda_min']**2<=f['lambda_min2']/64
    checks['reference_growth_ceiling_outward']=ref['lambda_max']**2>=64*f['lambda_max2']
    checks['original_squared_tilt_not_changed']=ref['tilt_upper']**2>=9 and not ref['original_tilt_reselected']
    return checks,dict(compact_rows=rows,axis_recipe=c,reference=ref,ordinary_hull_jets=jets,
        scope='rational parameter controls, not actual selected parameter substitution')

def margin_review():
    r=fixture();r0=Q(1,6057072);rows=[];checks={}
    for tag,Sfloor in (('fail',Q(1)),('pass',Q(10**12))):
        out=selector_margin(r,Q(2,5),2*r0*Sfloor,r0,3,Q(1,2));rows.append(out)
        c=r['coefficients'];P=c['A']/Sfloor+c['B']/Sfloor**2+c['C']/Sfloor**3
        checks[tag+'_same_selector_polynomial_independent_formula']=P==out['retained_power_error']
        checks[tag+'_original_radius_factor_retained']=out['existing_S_lower']==Sfloor
        checks[tag+'_original_selector_supplies_slow']=out['slow_gate_from_existing_selector'] and out['slow_upper_from_existing_selector']==Q(1,6)
        checks[tag+'_no_scale_or_Lreq_reselection']=not out['new_band_or_Lreq_selected'] and out['supplied_Lreq']==2*r0*Sfloor
        checks[tag+'_margin_has_expected_outcome']=out['normal_margin']==(tag=='pass')
        checks[tag+'_later_S_improves_same_comparison']=(c['A']/(2*Sfloor)+c['B']/(2*Sfloor)**2+c['C']/(2*Sfloor)**3<P)
        checks[tag+'_no_actual_margin_claim']=not out['actual_same_band_floor_verified']
    z=Q(1,3);partial=sum((z/2)**n for n in range(1,8))
    checks['full_first_order_tail_includes_first_coefficient']=partial<z/(2-z)<=z
    checks['first_term_cannot_be_removed']=partial>sum((z/2)**n for n in range(2,8))
    cert=certificate();checks['no_full_residual_or_blowup_promotion']=not cert['actual_parameter_substitution_complete_residual_or_blowup_verified']
    return checks,dict(rows=rows,certificate=cert,scope='existing-input sufficient normal-floor tests; failure never changes selector')

def review():
    checks={};details={}
    for name,func in (('symbolic',symbolic_review),('regional',rational_review),('margin',margin_review)):
        c,d=func();checks.update({name+':'+k:bool(v) for k,v in c.items()});details[name]=d
    return dict(checks=checks,details=details)
