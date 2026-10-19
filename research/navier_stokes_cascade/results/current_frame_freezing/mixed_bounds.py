"""Current background frame freezing from two endpoint value errors, not their derivatives."""
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
    p=HERE.parent/name/file;s=importlib.util.spec_from_file_location('current_freeze_'+name,p)
    m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
def leading():return module('leading_frame_freezing','freezing_bounds.py')
def deficit():return module('same_tilt_deficit','deficit_bounds.py')

def current_costs(data,JF,JG,Dbox,dchart):
    """Leading profile floors; after profile-normalized value error <=base_radius.

    Both current endpoints are compared with their own leading value. Current
    frames use unchanged u and the inherited positive |c0| branch.
    """
    JF,JG,Dbox,dchart=map(rat,(JF,JG,Dbox,dchart));p=data['inputs']
    if min(JF,JG,Dbox)<0 or dchart<=0:raise ValueError('Nonnegative leading jets/box cost and positive chart distance required')
    low=p['d']/(4*(p['Fbar']+1)*(p['Gbar']+1))
    high=2*((p['Gbar']+1)**2+2*(p['Fbar']+1)*(p['Gbar']+1))/(p['f']*p['t'])
    cmin=min(F(1),low);cmax=max(F(1),high)
    args=(p['f']/2,p['Fbar']+1,p['t']/2,p['Gbar']+1,p['omega']/3,cmin,cmax,p['u'])
    slow=leading().frame_cost(*args,JF,JG);noise=leading().frame_cost(*args,2,2)
    return dict(endpoint_inputs=slow['inputs'],mixing_squared_lower=low,mixing_squared_upper=high,
        slow_coordinate_cost=slow,value_error_cost=noise,
        Cmesh=Dbox*slow['relative_coordinate_Lipschitz'],Cpert=noise['relative_coordinate_Lipschitz'],
        chart_distance_threshold=Dbox/dchart,profile_normalization=F(8))

def mixed_gate(data,cost,E,ell):
    E=rat(E)
    if E<0 or type(ell) is not int or ell<2:raise ValueError('Nonnegative band error and integer ell>=2 required')
    eps=cost['profile_normalization']*E;w=data['inputs']['omega'];k=data['kappa']
    mesh=cost['Cmesh']/F(ell**6);pert=cost['Cpert']*eps;tau=mesh+pert
    background=2*data['Lpsi']*eps/(w*w);freeze=72*k*tau/(w*w);need=max(background,freeze)
    return dict(band_error=E,profile_error=eps,mesh_freezing_error=mesh,perturbation_freezing_error=pert,
        relative_freezing_error_upper=tau,background_deficit=background,freezing_deficit=freeze,
        required_deficit=need,raw_squared_supremum_required_upper=k-need,
        background_base_gate=eps<=data['base_radius'],
        chart_gate=F(ell**6)>=cost['chart_distance_threshold'],deficit_requirement_within_possible_range=need<=k)

def fixed_scale(data,cost,rows,M,N,h,ell):
    d=deficit();M,N=map(rat,(M,N));scale=d.dyadic_scale(h,ell)
    C=d.gates().cstar(rows,M,N);out=mixed_gate(data,cost,C/scale['factor'],ell)
    return dict(Cstar=C,dyadic=scale,shape_capacity=d.shapes().capacity(M,N),
        shape_budgets_possible_for_some_smooth_cutoff=d.shapes().capacity(M,N)>1,**out)

def later_witness(data,cost,rows,M,N,slack,h):
    """Conditional later band; does not change an old field or verify its actual inputs."""
    d=deficit();c=d.certificate(data,slack);gamma=c['floors']['unfrozen_cone_margin']
    if d.shapes().capacity(M,N)<=1:raise ValueError('Shape enclosures cannot describe any smooth cutoff')
    cap=min(c['radius']/8,gamma/(128*cost['Cpert']))
    tail=d.gates().dyadic_witness(rows,M,N,cap,h)
    mesh=d.sixth_root_ceiling(16*cost['Cmesh']/gamma)
    chart=d.sixth_root_ceiling(cost['chart_distance_threshold'])
    ell=max(2,tail['sufficient_ell'],mesh,chart);out=fixed_scale(data,cost,rows,M,N,h,ell)
    return dict(sufficient_ell=ell,background_witness=tail,mesh_ell=mesh,chart_ell=chart,
        band_error_cap=cap,unfrozen_margin=gamma,frozen_margin=gamma/2,fixed_scale=out,
        checks=dict(current_endpoint_value_and_cone_radius=out['profile_error']<=c['radius'],
            positive_chart_geometry=out['chart_gate'],mesh_spends_at_most_one_sixteenth=out['mesh_freezing_error']<=gamma/16,
            value_perturbation_spends_at_most_one_sixteenth=out['perturbation_freezing_error']<=gamma/16,
            combined_freezing_half_margin=out['relative_freezing_error_upper']<=gamma/8,
            both_deficits_fit_supplied_raw_slack=out['required_deficit']<=rat(slack)))

def state_raw(F0,c,z):
    F0,c,z=map(rat,(F0,c,z));co=(1-z*z)/(1+z*z);si=2*z/(1+z*z)
    mag=2*F0*co*(1+c*c)
    return leading().rational_state(F0,c,z),(-mag*co,mag*si),mag
def norm2(v):return sum(x*x for x in v)
def difference(v,w):return tuple(x-y for x,y in zip(v,w))

def endpoint_review():
    l=leading();g=l.gap();n0=(F(24,25),F(-7,25));rows=[]
    for u,w in ((F(3,4),F(5,4)),(F(4,3),F(5,3)),(F(5,12),F(13,12))):
        data=deficit().constants(F(1,2),2,2,6,4,F(1,2),u);c=deficit().certificate(data,F(1,100));gamma=c['floors']['unfrozen_cone_margin']
        cost=current_costs(data,1,20,6,F(1,20000));Cslow=cost['slow_coordinate_cost'];Cnoise=cost['value_error_cost']
        for d,small,(sf,sc,sz) in product((F(1,10**20),F(1,10**50)),(F(1,10**20),F(1,10**50)),product((-1,1),repeat=3)):
            eps=c['radius']*small;f0=1+sf*d/4;c0=F(4,5)+sc*d/4;z0=sz*d/4
            ref0,rawref0,magref0=state_raw(1,F(4,5),0);point0,raw0,mag0=state_raw(f0,c0,z0)
            fv=f0+sf*eps/32;cv=c0+sc*eps/32;zv=z0+sz*eps/32
            fr=1-sf*eps/32;cr=F(4,5)-sc*eps/32;zr=-sz*eps/32
            point,raw,mag=state_raw(fv,cv,zv);ref,rawref,magref=state_raw(fr,cr,zr)
            r=(F(sf)*eps/32,F(sz)*eps/32);target=tuple(x+y for x,y in zip(n0,r));H=-g.dot(point['N'],target);J=g.dot(point['K'],target)
            oldH=-g.dot(point0['N'],n0);oldJ=g.dot(point0['K'],n0)
            B=Cslow['relative_coordinate_Lipschitz']*d+Cnoise['relative_coordinate_Lipschitz']*eps
            invariant=None
            for m in (F(1),F(1,10**120)):
                T=tuple(m*x for x in target);V0,W0=l.coordinates(point,cv,w,u,T);Vf,Wf=l.coordinates(ref,cr,w,u,T)
                ratios=(abs(Vf-V0)/V0,abs(Wf-W0)/V0)
                if invariant is None:invariant=ratios
                rows.append(dict(u=str(u),distance=str(d),profile_error=str(eps),signs=[sf,sc,sz],target_scale=str(m),checks=dict(
                    leading_endpoint_and_global_floor_controls=F(1,2)<=f0<=2 and -raw0[0]>=2 and mag0<=6 and point0['D']>=4 and oldH>=F(1,2) and point0['c2']*oldJ**2<=(data['kappa']-F(1,100))*oldH**2,
                    leading_two_point_jets=abs(f0-1)<=d and norm2(difference(raw0,rawref0))<=(20*d)**2,
                    both_current_endpoint_value_errors=abs(fv-f0)<=eps and abs(fr-1)<=eps and norm2(difference(raw,raw0))<=eps**2 and norm2(difference(rawref,rawref0))<=eps**2,
                    pointwise_relative_target_error=norm2(r)<=eps**2,
                    current_two_point_shear_bound=norm2(difference(raw,rawref))<=(20*d+2*eps)**2,
                    current_endpoint_floors=all(F(1,4)<=f<=3 and -rr[0]>=1 and mm<=7 and ss['D']>=2 and cost['mixing_squared_lower']<=cc*cc<=cost['mixing_squared_upper'] for f,rr,mm,ss,cc in ((fv,raw,mag,point,cv),(fr,rawref,magref,ref,cr))),
                    actual_unit_direction_orientation_spends_omega_over_three=H>0 and H*H>=F(1,6)**2*norm2(target),
                    unfrozen_same_tilt_margin=(1+u*u)*cv*cv*J*J<=(1-gamma)**2*u*u*H*H,
                    mixed_relative_V_bound=ratios[0]<=B,mixed_relative_W_bound=ratios[1]<=B,
                    mixed_freezing_budget_fits_half_margin=B<=gamma/8,
                    frozen_same_tilt_half_margin=Vf>0 and abs(Wf)<=(1-gamma/2)*Vf,
                    target_amplitude_cancels_exactly=ratios==invariant)))
    return dict(auxiliary_endpoint_controls_not_actual_NS_profiles=True,case_count=len(rows),rows=rows,
        checks=dict(all_192_exact_current_endpoint_controls=all(all(row['checks'].values()) for row in rows),
            pointwise_target_has_both_components=all(x!=0 for x in n0)))

def obstruction_review():
    """Smooth small-amplitude frame oscillations: no derivative bound follows from C0."""
    rows=[];l=leading();u=F(3,4);w=F(5,4);c=F(4,5);old,raw0,mag0=state_raw(1,c,0)
    for a in (F(1,10),F(1,100),F(1,1000),F(1,10**8)):
        point,raw,mag=state_raw(1,c,a);ref,rawref,magref=state_raw(1,c,-a);dist=a*a;eps=8*a
        V,W=l.coordinates(point,c,w,u,(F(1),F(0)));Vf,Wf=l.coordinates(ref,c,w,u,(F(1),F(0)))
        relative=abs(Wf-W)/V;ratio=relative/dist
        rows.append(dict(rotation=str(a),distance=str(dist),value_error_upper=str(eps),relative_W_error=str(relative),relative_error_over_distance=str(ratio),checks=dict(
            both_endpoints_C0_close=norm2(difference(raw,raw0))<=eps**2 and norm2(difference(rawref,raw0))<=eps**2,
            exact_nonvanishing_error_per_value_error=relative==4*a*c*w/(u*(1-a*a)),
            ratio_grows_at_least_inverse_amplitude=ratio>1/a,
            both_frames_strictly_positive=point['D']>0 and ref['D']>0)))
    return dict(manufactured_frames_not_a_PDE_counterexample=True,rows=rows,checks=dict(
        value_only_perturbation_cannot_supply_pure_distance_Lipschitz_bound=all(all(row['checks'].values()) for row in rows),
        shrinking_value_error_and_distance_but_growing_slope=all(F(a['relative_error_over_distance'])<F(b['relative_error_over_distance']) and F(a['value_error_upper'])>F(b['value_error_upper']) for a,b in zip(rows,rows[1:]))))

def symbolic_review():
    import sympy as s
    E,L,k,w,Cm,Cp,ell=s.symbols('E L k w Cm Cp ell',positive=True)
    eps=8*E;tau=Cm/ell**6+Cp*eps
    x,z,b=s.symbols('x z b',positive=True);angular=(x+b)**2/((x+b)**2+z*z)-x*x/(x*x+z*z)
    v,d0=s.symbols('v d0',positive=True);gamma=w*w*d0/(9*k)
    return dict(checks=dict(
        normalization_eight_below_safe_three_halves_power=4**3==8**2,
        background_deficit_has_sixteen=s.simplify(2*L*eps/w**2-16*L*E/w**2)==0,
        mixed_freezing_deficit_has_both_terms=s.simplify(8*tau/gamma-72*k*(Cm/ell**6+8*Cp*E)/(w*w*d0))==0,
        split_mesh_and_perturbation_spend_one_eighth=s.Rational(1,16)+s.Rational(1,16)==s.Rational(1,8),
        nonnegative_outer_angular_stress_improves_orientation=s.simplify(angular-z*z*b*(2*x+b)/(((x+b)**2+z*z)*(x*x+z*z)))==0,
        nonnegative_outer_angular_stress_reduces_cone_ratio=s.simplify(z*z/(x+b)**2-z*z/x**2+z*z*b*(2*x+b)/(x*x*(x+b)**2))==0,
        outer_stress_magnitude_ratio_can_be_unbounded=s.limit(s.sqrt((x+b)**2+z*z)/s.sqrt(x*x+z*z),b,s.oo)==s.oo,
        rational_mixing_enclosures_are_safe=all(min(F(1),q)**2<=q<=max(F(1),q)**2 for q in (F(1,100),F(1),F(121,4)))))

def join_review():
    l=leading();old=l.join_review();d=deficit();g=d.gates();p=g.fixture();rows=g.affine_gate(**p)['rows']
    data=d.constants(F(1,2),2,2,6,4,F(1,2),F(3,4));chart=l.chart_caps(1,4)
    costs=current_costs(data,F(old['jets']['slow_F_Lipschitz']),F(old['jets']['slow_g_Lipschitz']),6,chart['admissible_slow_distance'])
    witness=later_witness(data,costs,rows,32,1216,F(1,100),p['h']);ell=witness['sufficient_ell']
    no_pressure=g.affine_gate(**dict(p,leading=(*p['leading'][:5],F(0),F(0),F(0))))['rows']
    erased=fixed_scale(data,costs,no_pressure,32,1216,p['h'],ell)
    wrong=data['inputs']['omega'];actual=costs['endpoint_inputs']['omega']
    return dict(auxiliary_complete_jet_pressure_composition_not_retained_band_acceptance=True,
        leading_profile_jets=old['jets'],costs=serial(costs),witness=serial(witness),checks=dict(
            all_joint_current_frame_gates=all(witness['checks'].values()),
            current_orientation_is_omega_over_three=actual==wrong/3,
            doubled_box_and_positive_chart_retained=costs['chart_distance_threshold']==6/chart['admissible_slow_distance'],
            full_pressure_and_viscosity_rows_increase_actual_error_budget=witness['fixed_scale']['Cstar']>erased['Cstar'] and witness['fixed_scale']['relative_freezing_error_upper']>erased['relative_freezing_error_upper'],
            mixed_value_error_cost_cannot_be_erased=costs['Cpert']>0 and witness['fixed_scale']['perturbation_freezing_error']>0,
            no_current_derivative_or_unit_target_derivative_input='current_derivative' not in costs))

def review():
    e,o,s,j=endpoint_review(),obstruction_review(),symbolic_review(),join_review()
    return dict(endpoints=e,obstruction=o,symbolic=s,join=j,checks=dict(
        exact_current_frame_mixed_freezing_and_same_tilt_margin=all(e['checks'].values()),
        value_only_pure_mesh_obstruction=all(o['checks'].values()),
        band_normalization_mixed_deficit_and_outer_sign_algebra=all(s['checks'].values()),
        complete_jet_pressure_and_current_frame_scale_join=all(j['checks'].values())))
