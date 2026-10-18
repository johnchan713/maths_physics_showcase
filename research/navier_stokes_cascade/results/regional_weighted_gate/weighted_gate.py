"""Finite-seminorm Borel bound and regional gate; diagnostics are not actual fields."""
from fractions import Fraction as F
from itertools import product
import importlib.util
from pathlib import Path

def rat(x):
    if isinstance(x,(bool,float)):raise TypeError('Exact rational data required')
    return F(x)

def tail(z,start=2):
    z=rat(z)
    if not isinstance(start,int) or isinstance(start,bool) or start<2 or not 0<z<=1:
        raise ValueError('0<q^h<=1 and integer start>=2 required')
    return (z/2)**start/(1-z/2)

def value_constant(E1,EX1,UX1,W1,Rlo,Rhi):
    E1,EX1,UX1,W1,Rlo,Rhi=map(rat,(E1,EX1,UX1,W1,Rlo,Rhi))
    if min(E1,EX1,UX1,W1)<0 or Rlo<=0 or Rhi<Rlo:raise ValueError('Nonnegative seminorms and positive ordered radii required')
    CE,CEX,CUX,CW=E1+F(1,2),EX1+F(1,2),UX1+F(1,2),W1+F(1,2)
    CF=5*CE/Rlo;Cgtheta=5*(Rhi*CEX+CE/Rlo);Cgz=5*Rhi*CUX
    return dict(E=CE,EX=CEX,UX=CUX,W=CW,F=CF,gtheta=Cgtheta,gz=Cgz,
        g_Euclidean_sufficient=Cgtheta+Cgz,Cval=max(F(1),CF,Cgtheta+Cgz,CW))

def schedule(n,h,Cvelocity,Cstress,previous_exp=0):
    h,Cvelocity,Cstress=map(rat,(h,Cvelocity,Cstress))
    if not isinstance(n,int) or isinstance(n,bool) or n<2 or not 0<h<=F(1,100) or min(Cvelocity,Cstress)<0:
        raise ValueError('Finite n>=2, retained positive small h and finite nonnegative bounds required')
    if not isinstance(previous_exp,int) or isinstance(previous_exp,bool) or previous_exp<0:raise ValueError('Nonnegative previous dyadic exponent required')
    C=max(F(1),Cvelocity,Cstress);p=0
    while 2**p<C:p+=1
    requirement=F(p+n)/(n*h);ceil_req=-(-requirement.numerator//requirement.denominator)
    exponent=max(previous_exp+1,ceil_req+1)
    return dict(n=n,h=h,C=C,ceil_log2_C=p,dyadic_exponent=exponent,
        logarithmic_margin=exponent*n*h-p-n)

def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x

def summation_review():
    rows=[]
    for z,j in product((F(1),F(1,2),F(1,8),F(1,10**12)),(2,3,7)):
        total=tail(z,j);last=j+12
        finite=sum(((z/2)**n for n in range(j,last+1)),F(0))
        remainder=(z/2)**(last+1)/(1-z/2)
        rows.append(dict(z=str(z),start=j,checks=dict(exact_geometric_identity=(1-z/2)*total==(z/2)**j,
            finite_plus_remainder=total==finite+remainder,whole_small_q_half_power=total<=z*z/2)))
    schedules=[];previous=0
    for n in range(2,12):
        r=schedule(n,F(1,200),2**(n*n),2**(n*n+3),previous)
        checks=dict(both_budgets_included=r['C']>=2**(n*n+3),strict_full_support_margin=r['logarithmic_margin']>0,
            geometric_growth=r['dyadic_exponent']>=previous+1)
        previous=r['dyadic_exponent'];schedules.append(dict(certificate=serial(r),checks=checks))
    constants=value_constant(46,19,6,3,1,3)
    # Velocity-only bounds permit an arbitrarily large independently supplied stress coefficient.
    z=F(1,4);W2=16;uncertified_relative=W2*z**4
    sigma=F(2**300,2**400)
    missed_stress=dict(q_h=str(z),q='2^-400',cutoff_scale='2^300',cutoff_argument=str(sigma),
        cutoff_equals_one=sigma<=F(1,2),padded_velocity_support_budget_passes=F(1,8)<=F(1,4),
        stress_support_budget_fails=F(W2,8)>F(1,4),velocity_tail='0',stress_order_two=str(W2),
        stress_value=str(uncertified_relative),claimed_stress_tail=str(z*z/2),
        violates_velocity_only_inference=uncertified_relative>z*z/2)
    return dict(geometric_rows=rows,shared_dyadic_schedules=schedules,finite_seminorm_constants=serial(constants),
        independent_stress_budget_control=missed_stress,
        checks=dict(all_exact_geometric_controls=all(all(r['checks'].values()) for r in rows),
            all_shared_support_budgets=all(all(r['checks'].values()) for r in schedules),
            stress_seminorm_cannot_be_inferred_from_velocity=missed_stress['violates_velocity_only_inference'],
            componentwise_vector_norm_pad_retained=F(2)<=F(2)**2 and F(2)>F(1)**2,
            uniform_chart_five_factor=5**8>4**9 and F(1,100)<F(1,8)))

def chart_review():
    rows=[];C=value_constant(46,19,6,3,1,3);z=F(1,4);power=z*z
    for rp,ss,B in product((F(1),F(2),F(3)),(F(1,2),F(1),F(2)),(F(1,2),F(1),F(2))):
        s=ss*ss;R=rp*ss;X=rp*rp/2
        e=power*(1+X+2*X*X);ex=power*(1+4*X);u=power*(-3*X+X*X);ux=power*(-3+2*X)
        fv=B*e/R;FR=B*(ex/s-e/(R*R));gz=B*R*ux/s;gtheta=R*FR
        prefactor=B/ss
        checks=dict(F_profile_identity=fv==prefactor*e/rp,
            theta_radial_chain_rule=gtheta==prefactor*(rp*ex-e/rp),
            axial_radial_chain_rule=gz==prefactor*rp*ux,
            nonunit_chart_factor_retained=prefactor<=5,
            F_constant=abs(fv)<=C['F']*power,
            theta_constant=abs(gtheta)<=C['gtheta']*power,
            axial_constant=abs(gz)<=C['gz']*power,
            Euclidean_constant_sufficient=gtheta*gtheta+gz*gz<=(C['g_Euclidean_sufficient']*power)**2,
            missing_inverse_radius_derivative_detected=gtheta!=B*R*ex/s)
        rows.append(dict(Rprof=str(rp),sqrt_s=str(ss),prefactor_B=str(B),checks=checks))
    return dict(model='rational fixed-s chain rules with an independent positive band prefactor; not fractional-power actual profiles',
        rows=rows,checks=dict(all_exact_chart_and_constant_controls=all(all(r['checks'].values()) for r in rows)))

def outer_review():
    rows=[]
    for angular,axial,first,tau,u in product((F(1),F(2)),(F(-1,5),F(0),F(1,5)),
        (F(0),F(1),F(1000)),(F(0),F(1),F(2**80)),(F(1,10),F(1))):
        c2=F(1,200);old=(angular,axial);new=(angular+tau*first,axial+F(0));newangular=new[0]
        old_ratio=(1+u*u)*c2*axial*axial/(u*u*angular*angular)
        new_ratio=(1+u*u)*c2*axial*axial/(u*u*newangular*newangular)
        old_orientation2=angular*angular/(angular*angular+axial*axial)
        new_orientation2=newangular*newangular/(newangular*newangular+axial*axial)
        checks=dict(leading_cone_accepted=old_ratio<1,angular_orientation_positive=newangular>0,
            same_tilt_squared_cone_ratio_decreases=new_ratio<=old_ratio,
            unit_orientation_does_not_decrease=new_orientation2>=old_orientation2,
            axial_stress_is_unchanged=new[1]==old[1])
        rows.append(dict(angular=str(angular),axial=str(axial),first=str(first),tau=str(tau),u=str(u),checks=checks))
    h=F(1,200);Xb=F(100);edge=2*Xb/(1-2*h)**2
    z=F(1,2);relative=edge*z*z;model_radius=F(1,68992)
    oldangular=F(1);newangular=oldangular+relative;axial=F(1,5);u=F(1);c2=h
    old_ratio=2*c2*axial*axial;new_ratio=2*c2*axial*axial/(newangular*newangular)
    # A signed cutoff can reverse the angular orientation and is not covered.
    signed_after=F(1)+F(-2)*F(1)
    return dict(model='algebraic use of reviewed exact outer support/sign identity, not a new heat-tail derivation',rows=rows,
        large_edge_control=dict(actual_edge_formula='2*Xb/(1-2h)^2 at eta=+/-1, with first cutoff one',
            diagnostic_edge_coefficient=str(edge),relative_size=str(relative),model_radius=str(model_radius),
            sufficient_global_small_error_test_fails=relative>model_radius,
            nevertheless_same_tilt_cone_improves=new_ratio<=old_ratio<1),
        signed_cutoff_orientation_control=str(signed_after),
        checks=dict(all_outer_same_tilt_controls=all(all(r['checks'].values()) for r in rows),
            large_error_not_cone_failure=relative>model_radius and new_ratio<old_ratio<1,
            nonnegative_cutoff_requirement_is_essential=signed_after<0))

def coupled_review():
    path=Path(__file__).resolve().parent.parent/'background_gap_certificate/gap_bounds.py'
    spec=importlib.util.spec_from_file_location('regional_parent_gap_formulas',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    parent=module.certificate(F(1,2),2,2,5,4,F(1,2),F(1,2),1)
    C=value_constant(46,19,6,3,1,3)['Cval'];j=0
    while C/F(2)**(2*j)>parent['radius']:j+=1
    z=F(1,2**j);power=z*z;delta=C*power
    # Compatible local jets at Rprof=s=1: F0=R^-4, G0=0,
    # e=q^(2h)(1+X+2X^2), u=q^(2h)(-3X+X^2), X=1/2.
    fv=1+2*power;gt=-4+power;gz=-2*power
    target=(F(24,25)-power,F(7,25)+power)
    g2=gt*gt+gz*gz;D=g2+2*fv*gt
    lam2=-2*fv*gt*D/g2;c2=D/(-2*fv*gt)
    Hnum=-gt*target[0]-gz*target[1];Jnum=-gz*target[0]+gt*target[1]
    target2=target[0]**2+target[1]**2;floor=parent['floors']
    checks=dict(combined_value_gate=delta<=parent['radius'],
        previous_grid_point_not_certified=C/F(2)**(2*(j-1))>parent['radius'],
        exact_input_norms_fit_gate=2*power<=delta and 5*power*power<=delta*delta and 2*power*power<=delta*delta,
        actual_spectral_floor=lam2>=floor['lambda2'],actual_mixing_floor=c2>=floor['c2'],
        positive_oriented_unit_floor=Hnum>0 and Hnum*Hnum>=floor['oriented_unit_stress']**2*g2*target2,
        same_tilt_margin=2*c2*Jnum*Jnum<=(1-floor['cone_margin'])**2*Hnum*Hnum)
    return dict(model='compatible manufactured local radial jets, not actual selected fields or a global supremum certificate',
        diagnostic_h='1/200',q_point_exponent=200*j,full_chart_band_exponent=200*j+2,
        Cval=str(C),radius=str(parent['radius']),delta=str(delta),lambda2=str(lam2),c2=str(c2),checks=checks)

def review():
    s=summation_review();c=chart_review();o=outer_review();g=coupled_review()
    return dict(shared_Borel_tail_and_finite_seminorms=s,exact_band_radial_pullback=c,outer_monotonic_cover=o,
        composed_local_input_to_parent_gap=g,
        checks=dict(shared_Borel_value_bound=all(s['checks'].values()),
            exact_radial_chart_interface=all(c['checks'].values()),outer_sign_cover=all(o['checks'].values()),
            composed_local_gap_controls=all(g['checks'].values())))
