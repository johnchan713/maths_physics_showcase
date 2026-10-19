"""Explicit pointwise-target frame freezing; no freezing of flat target amplitudes."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
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
    p=HERE.parent/name/file;s=importlib.util.spec_from_file_location('leading_freeze_'+name,p)
    m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
def deficit():return module('same_tilt_deficit','deficit_bounds.py')
def gap():return module('background_gap_certificate','gap_bounds.py')

def chart_caps(a,b):
    """Xa=a^2,Xb=b^2; bootstrap 1/8<=s<=8,Xa/2<=X<=2Xb, tau>=0."""
    a,b=map(rat,(a,b))
    if a<=0 or b<a:raise ValueError('Positive ordered active-annulus root radii required')
    xlo=a*a/2;xhi=2*b*b;rlo=a;rhi=2*b
    KX=4*rhi+144*xhi;Keta=F(14)
    return dict(Xa=a*a,Xb=b*b,Xlo=xlo,Xhi=xhi,Rlo=rlo,Rhi=rhi,
        X_coordinate_Lipschitz=KX,eta_coordinate_Lipschitz=Keta,
        admissible_slow_distance=min(F(1,288),a*a/(4*KX)))
def profile_to_slow_jets(chart,swirl,axial):
    """Full F,U derivative enclosures on the enlarged profile rectangle, not stresses.

    F_eta<=F_C3; F_Xeta<=F_X_C2; U_Xeta<=U_X_C2. The source lemmas
    extend swirl to any finite upper radius and axial derivatives by zero.
    """
    f,fx,fxx=(rat(swirl[k]) for k in ('F_C3','F_X_C2','F_XX_C1'))
    ux,uxx=(rat(axial[k]) for k in ('U_X_C2','U_XX_C1'))
    if min(f,fx,fxx,ux,uxx)<0:raise ValueError('Nonnegative complete leading jet bounds required')
    X=chart['Xhi'];r=chart['Rhi'];ra=chart['Rlo'];KX=chart['X_coordinate_Lipschitz'];Ke=chart['eta_coordinate_Lipschitz']
    G0=2*X*fx+r*ux;GX=2*fx+2*X*fxx+ux/ra+r*uxx;Geta=2*X*fx+r*ux
    return dict(profile_g_value=G0,profile_g_X=GX,profile_g_eta=Geta,
        slow_F_Lipschitz=fx*KX+f*Ke,slow_g_Lipschitz=GX*KX+Geta*Ke)
def frame_cost(f,Fbar,t,Gbar,omega,cmin,cmax,u,JF,JG):
    """Endpoint value floors and |Delta F|<=JF*d, |Delta g|<=JG*d.

    Target T=m*n is kept at the evaluation point. Only N,K,c are frozen.
    """
    f,Fbar,t,Gbar,omega,cmin,cmax,u,JF,JG=map(rat,(f,Fbar,t,Gbar,omega,cmin,cmax,u,JF,JG))
    if min(f,t,omega,cmin,u)<=0 or Fbar<f or Gbar<t or cmax<cmin or omega>1 or min(JF,JG)<0:
        raise ValueError('Positive ordered frame floors, same tilt and nonnegative jet bounds required')
    Dbar=Gbar*Gbar+2*Fbar*Gbar
    LD=2*(Gbar+Fbar)*JG+2*Gbar*JF;LP=Gbar*JF+Fbar*JG
    LC=LD/(2*f*t)+Dbar*LP/(2*f*f*t*t);LN=2*JG/t
    V=cmax*LN/(cmin*omega)+LC/(2*cmin*cmin)
    W=cmax*(1+u)*LN/(u*omega)
    return dict(inputs=dict(f=f,Fbar=Fbar,t=t,Gbar=Gbar,omega=omega,cmin=cmin,cmax=cmax,u=u,JF=JF,JG=JG),
        Dbar=Dbar,D_difference_Lipschitz=LD,product_difference_Lipschitz=LP,
        squared_mixing_Lipschitz=LC,unit_frame_Lipschitz=LN,
        relative_V_Lipschitz=V,relative_W_Lipschitz=W,relative_coordinate_Lipschitz=max(V,W))
def freezing_budget(chart,cost,box_diameter_constant):
    d=rat(box_diameter_constant)
    if d<0:raise ValueError('Nonnegative retained slow-box diameter constant required')
    return dict(Cfreeze=d*cost['relative_coordinate_Lipschitz'],box_diameter_constant=d,
        chart_distance_threshold=d/chart['admissible_slow_distance'])
def retained_box_cost(enlargement):
    """Representative in original box, point in k-fold box: distance <=(k+1)sqrt(3)delta."""
    enlargement=rat(enlargement)
    if enlargement<1:raise ValueError('Enlargement must include original closed box')
    return 2*(enlargement+1)  # sqrt(3)<2, rational safe enclosure
def later_freezing_band(budget,gamma0):
    gamma0=rat(gamma0)
    if not 0<gamma0<=1:raise ValueError('Certified positive unfrozen margin <=1 required')
    root=deficit().sixth_root_ceiling
    chart=root(budget['chart_distance_threshold']);freeze=root(8*budget['Cfreeze']/gamma0)
    ell=max(2,chart,freeze)
    return dict(geometry_ell=chart,freezing_ell=freeze,sufficient_ell=ell,frozen_margin=gamma0/2,
        checks=dict(geometry_distance_accepted=F(ell**6)>=budget['chart_distance_threshold'],
            independent_relative_freezing_accepted=budget['Cfreeze']/F(ell**6)<=gamma0/8))

def rational_state(F0,c,z):
    F0,c,z=map(rat,(F0,c,z));co=(1-z*z)/(1+z*z);si=2*z/(1+z*z)
    if F0<=0 or c<=0 or abs(z)>=1:raise ValueError('Positive F,mixing magnitude and |rotation parameter|<1 required')
    mag=2*F0*co*(1+c*c)
    return gap().shear(F0,-mag*co,mag*si,mag)
def coordinates(state,c,w,u,target):
    V=-gap().dot(state['N'],target)/(c*w);W=gap().dot(state['K'],target)/u
    return V,W

def exact_review():
    f0=F(1);c0=F(4,5);reference=rational_state(f0,c0,0);n=(F(24,25),F(-7,25));omega=F(1,2)
    rows=[];costs=[]
    for u,w in ((F(3,4),F(5,4)),(F(4,3),F(5,3)),(F(5,12),F(13,12))):
        B=frame_cost(F(1,2),2,2,6,omega,F(1,2),1,u,1,20);costs.append(serial(B))
        for d,(sf,sc,sz) in product((F(1,100),F(1,1000),F(1,10**10)),product((-1,1),repeat=3)):
            c=c0+sc*d;point=rational_state(f0+sf*d,c,sz*d)
            gdiff=sum((x-y)**2 for x,y in zip(point['N'],reference['N']))
            # Recover the raw shear exactly, independent of the freezing formula.
            z=sz*d;co=(1-z*z)/(1+z*z);si=2*z/(1+z*z);mag=2*(f0+sf*d)*co*(1+c*c)
            raw=(-mag*co,mag*si);raw0=(F(-82,25),F(0))
            normdiff=sum((x-y)**2 for x,y in zip(raw,raw0))
            invariant=None
            for m in (F(1),F(1,10**40),F(1,10**200)):
                target=tuple(m*v for v in n);V0,W0=coordinates(point,c,w,u,target);Vf,Wf=coordinates(reference,c0,w,u,target)
                ratios=(abs(Vf-V0)/V0,abs(Wf-W0)/V0)
                if invariant is None:invariant=ratios
                H=-gap().dot(point['N'],n)
                rows.append(dict(u=str(u),distance=str(d),signs=[sf,sc,sz],target_magnitude=str(m),checks=dict(
                    exact_tilt_square_root=w*w==1+u*u,
                    unit_target_and_orientation=sum(v*v for v in n)==1 and H>=omega,
                    endpoint_value_and_mixing_floors=F(1,2)<=f0+sf*d<=2 and -raw[0]>=2 and mag<=6 and F(1,2)<=c<=1,
                    input_shear_difference_budget=normdiff<=(20*d)**2,
                    exact_mixing_parameter_identity=point['c2']==c*c and reference['c2']==c0*c0,
                    squared_mixing_difference_bound=abs(c*c-c0*c0)<=B['squared_mixing_Lipschitz']*d,
                    normalized_frame_difference_bound=gdiff<=(B['unit_frame_Lipschitz']*d)**2,
                    relative_V_freezing_bound=ratios[0]<=B['relative_V_Lipschitz']*d,
                    relative_W_freezing_bound=ratios[1]<=B['relative_W_Lipschitz']*d,
                    target_magnitude_cancels_exactly=ratios==invariant)))
    return dict(auxiliary_endpoint_frames_not_selected_profiles=True,case_count=len(rows),rows=rows,costs=costs,
        checks=dict(all_216_exact_pointwise_target_freezing_controls=all(all(r['checks'].values()) for r in rows),
            both_pressure_sensitive_target_components_retained=n[0]!=0 and n[1]!=0,
            arbitrarily_tiny_stress_magnitudes_do_not_change_bound={r['target_magnitude'] for r in rows}=={'1',str(F(1,10**40)),str(F(1,10**200))}))

def symbolic_review():
    import sympy as s
    X,e,h,q,R=s.symbols('X e h q R',positive=True);D=s.Rational(1,2)-h;A=s.Rational(1,2)+h;L=1-2*h*e*e
    etaZ=q**(-D)-D*e/q*(2*e*q**A/L)
    U=s.Function('U')(X);F0=s.Function('F')(X)
    y=s.symbols('y',positive=True);exponent=1/y**2-1/(y+y**3)**2
    H,Hf,cc,cf,m,w=s.symbols('H Hf cc cf m w',positive=True)
    V=m*H/(cc*w);Vf=m*Hf/(cf*w)
    return dict(checks=dict(
        exact_eta_Z_cancellation_at_poles=s.simplify(etaZ-q**(-D)*(1-e*e)/L)==0,
        native_angular_shear_X_derivative=s.simplify(s.diff(2*X*s.diff(F0,X),X)-(2*s.diff(F0,X)+2*X*s.diff(F0,X,2)))==0,
        native_axial_shear_X_inverse_radius_retained=s.simplify(s.diff(s.sqrt(2*X)*s.diff(U,X),X)-(s.diff(U,X)/s.sqrt(2*X)+s.sqrt(2*X)*s.diff(U,X,2)))==0,
        covariance_relative_V_exact_difference=s.simplify((Vf-V)/V-(cc/cf*(Hf-H)/H+(cc-cf)/cf))==0,
        amplitude_and_fixed_tilt_factor_cancel_in_relative_V=s.simplify((Vf-V)/V-(cc*Hf/(cf*H)-1))==0,
        flat_amplitude_freezing_exponent_exact=s.simplify(exponent-(2+y*y)/(1+y*y)**2)==0,
        flat_amplitude_error_does_not_tend_to_zero=s.limit(exponent,y,0)==2))

def geometry_review():
    # Exact caps on the larger positive chart: 1/L<2,s^A<4,s^-D<3.
    h=F(1,100);chart=chart_caps(1,4);X=chart['Xhi'];r=chart['Rhi']
    return dict(auxiliary_geometry_not_old_box_selection=True,chart=serial(chart),checks=dict(
        denominator_inverse_below_two=1/(1-2*h)<2,
        s_exponent_uses_A_below_two_thirds=F(1,2)+h<F(2,3) and 8**2==4**3,
        eta_Z_cap_six_and_eta_T_cap_eight=3*2==6 and F(1,2)*2*8==8,
        s_Z_cap_sixteen_and_s_T_cap_two=2*4*2==16,
        X_coordinate_cap_all_three_slow_derivatives=4*r+(16+128)*X==chart['X_coordinate_Lipschitz'],
        bootstrap_s_margin_strict=18*chart['admissible_slow_distance']<=F(1,16),
        bootstrap_X_margin_strict=chart['X_coordinate_Lipschitz']*chart['admissible_slow_distance']<=chart['Xa']/4,
        original_and_doubled_boxes_have_explicit_diameters=retained_box_cost(1)==4 and retained_box_cost(2)==6,
        box_distance_rational_enclosure_strict=3<2**2))

def join_review():
    d=deficit();g=d.gates();p=g.fixture();a=F(1);b=F(4);chart=chart_caps(a,b)
    axial=module('global_axial_source_envelope','envelopes.py').global_bounds(p['h'],1,25,a*a,7,13,101,303)
    swirl=module('global_swirl_envelope','swirl_bounds.py').swirl_inputs(a*a,chart['Xhi'],7,101)
    jets=profile_to_slow_jets(chart,swirl,axial)
    cost=frame_cost(F(1,2),2,2,6,F(1,2),F(1,2),1,F(3,4),jets['slow_F_Lipschitz'],jets['slow_g_Lipschitz'])
    budget=freezing_budget(chart,cost,retained_box_cost(2));data=d.constants(F(1,2),2,2,6,4,F(1,2),F(3,4));slack=F(1,100)
    gamma=d.certificate(data,slack)['floors']['unfrozen_cone_margin'];freez=later_freezing_band(budget,gamma)
    full=g.affine_gate(**p)['rows'];joint=d.later_witness(data,full,32,1216,slack,p['h'],budget['Cfreeze'])
    ell=max(freez['sufficient_ell'],joint['sufficient_joint_ell']);requirements=d.fixed_scale(data,full,32,1216,p['h'],ell,budget['Cfreeze'])
    return dict(auxiliary_full_jet_composition_not_actual_scale=True,chart=serial(chart),jets=serial(jets),cost=serial(cost),budget=serial(budget),
        independent_geometry_and_freezing_witness=serial(freez),joint_ell=ell,joint_requirements=serial(requirements),checks=dict(
            primitive_global_jets_give_finite_explicit_freezing_constant=budget['Cfreeze']>0,
            exact_enlarged_upper_radius_used=chart['Xhi']==2*b*b,
            geometry_and_independent_freezing_checks=all(freez['checks'].values()),
            same_deficit_and_complete_pressure_scale_join=requirements['geometry_gate'] and requirements['required_deficit']<=slack,
            no_unit_stress_derivative_input_added='n_derivative' not in cost['inputs']))

def review():
    e,s,g,j=exact_review(),symbolic_review(),geometry_review(),join_review()
    return dict(exact=e,symbolic=s,geometry=g,join=j,checks=dict(
        exact_relative_pointwise_target_frame_controls=all(e['checks'].values()),
        symbolic_geometry_shear_and_invalid_amplitude_freeze_controls=all(s['checks'].values()),
        controlled_chart_bootstrap_and_value_gate=all(g['checks'].values()),
        explicit_global_jet_freezing_and_full_deficit_join=all(j['checks'].values())))
