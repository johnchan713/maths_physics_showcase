"""Same unchanged tilt: normalized deficit radius and coupled scale/freezing requirements."""
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
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x
@lru_cache(None)
def module(name,file):
    p=HERE.parent/name/file;s=importlib.util.spec_from_file_location('deficit_'+name,p)
    m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m
def gap():return module('background_gap_certificate','gap_bounds.py')
def gates():return module('cutoff_cone_scale_gate','gate_bounds.py')
def shapes():return module('cutoff_shape_feasibility','shape_bounds.py')

def tilt_slack(u,rho2_upper):
    u,R=map(rat,(u,rho2_upper))
    if u<=0 or not 0<=R<=1:raise ValueError('Positive unchanged tilt and raw squared-ratio upper bound in [0,1] required')
    k=u*u/(1+u*u);slack=k-R
    return dict(kappa=k,slack=slack,admissible_by_supplied_bound=slack>0)
def constants(f,Fbar,t,Gbar,d,omega,u):
    f,Fbar,t,Gbar,d,omega,u=map(rat,(f,Fbar,t,Gbar,d,omega,u))
    if min(f,t,d,omega,u)<=0 or Fbar<f or Gbar<t or omega>1:
        raise ValueError('Positive ordered leading floors, bounds and 0<omega<=1 required')
    k=u*u/(1+u*u);LD=4*Gbar+2*Fbar+3;LP=Fbar+Gbar+1;Dbar=Gbar**2+2*Fbar*Gbar
    cbar2=Dbar/(2*f*t);Lc=2*LD/(f*t)+2*Dbar*LP/(f*f*t*t);LH=1+2/t
    Lpsi=3*k*LH+4*Lc+3*cbar2*LH
    base=dict(relative_stress=F(1,2),swirl=f/2,angular_shear=t/2,discriminant=d/(2*LD),orientation=omega/(2*LH))
    return dict(inputs=dict(f=f,Fbar=Fbar,t=t,Gbar=Gbar,d=d,omega=omega,u=u),kappa=k,
        LD=LD,LP=LP,Dbar=Dbar,cbar2=cbar2,Lc=Lc,LH=LH,Lpsi=Lpsi,
        base_constraints=base,base_radius=min(base.values()))
def certificate(data,slack):
    """Requires a VALID lower bound 0<slack<=kappa on kappa-rho_actual^2.

    This function does not compute the actual supremum or verify that premise.
    """
    slack=rat(slack);k=data['kappa'];w=data['inputs']['omega'];p=data['inputs']
    if not 0<slack<=k:raise ValueError('Positive raw deficit lower bound no larger than kappa required')
    rcone=w*w*slack/(2*data['Lpsi']);radius=min(data['base_radius'],rcone)
    gamma=w*w*slack/(9*k)
    return dict(slack=slack,radius=radius,cone_radius=rcone,normalized_cone_floor=w*w*slack/2,
        floors=dict(F=p['f']/2,minus_gtheta=p['t']/2,D=p['d']/2,
            normal_theta=p['t']/(2*(p['Gbar']+1)),lambda2=p['f']*p['t']*p['d']/(4*(p['Gbar']+1)**2),
            c2=p['d']/(4*(p['Fbar']+1)*(p['Gbar']+1)),oriented_unit_stress=w/3,
            unfrozen_cone_margin=gamma,frozen_cone_margin=gamma/2))
def dyadic_scale(h,ell):
    h=rat(h)
    if not 0<h<F(1,100) or type(ell) is not int or ell<2:raise ValueError('Admissible h and integer ell>=2 required')
    exponent=2*h*(ell-2);k=exponent.numerator//exponent.denominator
    return dict(exact_exponent=exponent,integer_floor=k,factor=F(2)**k)
def required_deficit(data,error,Cfreeze,ell):
    error,Cfreeze=map(rat,(error,Cfreeze));w=data['inputs']['omega'];k=data['kappa']
    if min(error,Cfreeze)<0 or type(ell) is not int or ell<2:raise ValueError('Nonnegative error/freezing costs and integer ell>=2 required')
    background=2*data['Lpsi']*error/(w*w);freeze=72*k*Cfreeze/(w*w*ell**6)
    need=max(background,freeze)
    return dict(background_deficit=background,freezing_deficit=freeze,required_deficit=need,
        raw_squared_supremum_required_upper=k-need,geometry_gate=error<=data['base_radius'],
        deficit_requirement_within_possible_range=need<=k)
def fixed_scale(data,rows,M,N,h,ell,Cfreeze):
    """Supplied M,N must separately bound the unchanged cutoff; capacity is only necessary."""
    rows=shapes().normalized_rows(rows);M,N=map(rat,(M,N));scale=dyadic_scale(h,ell)
    C=gates().cstar(rows,M,N);e=C/scale['factor'];r=required_deficit(data,e,Cfreeze,ell)
    return dict(Cstar=C,dyadic=scale,error_upper=e,shape_capacity=shapes().capacity(M,N),
        shape_budgets_are_possible_for_some_smooth_cutoff=shapes().capacity(M,N)>1,**r)
def sixth_root_ceiling(q):
    q=rat(q)
    if q<0:raise ValueError('Nonnegative exact sixth-root threshold required')
    if q<=1:return int(q>0)
    lo=0;hi=1<<((q.numerator.bit_length()-q.denominator.bit_length()+6)//6)
    while lo<hi:
        m=(lo+hi)//2
        if m**6*q.denominator>=q.numerator:hi=m
        else:lo=m+1
    return lo
def later_witness(data,rows,M,N,slack,h,Cfreeze):
    Cfreeze=rat(Cfreeze)
    if Cfreeze<0:raise ValueError('Nonnegative freezing constant required')
    if shapes().capacity(M,N)<=1:raise ValueError('Supplied shape budgets cannot bound any smooth cutoff')
    c=certificate(data,slack);old=gates().dyadic_witness(rows,M,N,c['radius'],h)
    ellF=sixth_root_ceiling(8*Cfreeze/c['floors']['unfrozen_cone_margin'])
    ell=max(old['sufficient_ell'],ellF,2);r=fixed_scale(data,rows,M,N,h,ell,Cfreeze)
    return dict(background_witness=old,freezing_ell=ellF,sufficient_joint_ell=ell,
        fixed_scale_requirements=r,checks=dict(
            whole_four_Q_chart_background_gate=r['error_upper']<=c['radius'],
            independent_freezing_gate=Cfreeze/F(ell**6)<=c['floors']['unfrozen_cone_margin']/8,
            both_required_deficits_fit_supplied_lower_bound=r['required_deficit']<=rat(slack)))
def auxiliary_shape_join(data,rows,slack,h,ell,Cfreeze):
    Cfreeze=rat(Cfreeze)
    if Cfreeze<0:raise ValueError('Nonnegative freezing constant required')
    c=certificate(data,slack);scale=dyadic_scale(h,ell);freeze=Cfreeze/F(ell**6)<=c['floors']['unfrozen_cone_margin']/8
    return dict(certificate=c,dyadic=scale,freezing_gate=freeze,
        auxiliary_shape_program=shapes().shape_gate(rows,c['radius']*scale['factor']))

def reference_target(z):
    """Exact unit stress rotated inwards from a fixed-u cone boundary; not a NS profile."""
    z=rat(z);co=(1-z*z)/(1+z*z);si=2*z/(1+z*z)
    return (F(4,5)*co+F(3,5)*si,-F(3,5)*co+F(4,5)*si)
def psi(state,target,kappa):
    H=-gap().dot(state['N'],target);J=gap().dot(state['K'],target)
    return rat(kappa)*H**2-state['c2']*J**2

def perturbation_review():
    old=gap().shear(1,F(-82,25),0,F(82,25));n0=reference_target(F(1,32));controls=[];comparisons=[]
    for u in (F(3,4),F(1),F(2)):
        data=constants(F(1,2),2,2,5,4,F(1,2),u);k=data['kappa'];q0=gap().cone(old,n0,u)
        R=old['c2']*q0['J']**2/q0['H']**2;slack=(k-R)/2;c=certificate(data,slack);e=c['radius']/2
        gamma=slack/(2*k);legacy=gap().certificate(F(1,2),2,2,5,4,F(1,2),gamma,u)
        comparisons.append(dict(u=str(u),checks=dict(
            normalized_Lipschitz_exactly_matches_legacy=data['Lpsi']==legacy['Lphi']/(1+u*u),
            all_five_noncone_constraints_unchanged=data['base_constraints']=={v:legacy['constraints'][v] for v in data['base_constraints']},
            new_radius_never_smaller=c['radius']>=legacy['radius'],
            direct_cone_radius_strictly_larger=c['cone_radius']>legacy['constraints']['same_tilt_cone'],
            new_cone_floor_margin_strictly_larger=c['floors']['unfrozen_cone_margin']>legacy['floors']['cone_margin'])))
        for sf,sm,sg,sx,sz in product((-1,1),repeat=5):
            a=sg*e/32;co=(1-a*a)/(1+a*a);si=2*a/(1+a*a);magnitude=F(82,25)+sm*e/4
            state=gap().shear(1+sf*e/3,-magnitude*co,magnitude*si,magnitude)
            target=(n0[0]+sx*e/4,n0[1]+sz*e/5);q=gap().cone(state,target,u);n2=sum(v*v for v in target)
            before=psi(old,n0,k);after=psi(state,target,k);gammaB=c['floors']['unfrozen_cone_margin']
            controls.append(dict(u=str(u),signs=[sf,sm,sg,sx,sz],checks=dict(
                reference_unit_direction=sum(v*v for v in n0)==1,
                leading_raw_deficit_is_valid=k-R>=slack>0,
                normalized_cone_polynomial_identity=after==q['phi']/(1+u*u),
                normalized_cone_Lipschitz=abs(after-before)<=data['Lpsi']*e,
                normalized_positive_cone_floor=after>=c['normalized_cone_floor'],
                spectral_and_swirl_floors=state['D']>=c['floors']['D'] and state['lambda2']>=c['floors']['lambda2'] and state['c2']>=c['floors']['c2'],
                relative_target_norm_not_reset=F(1,4)<=n2<=F(9,4),
                oriented_unit_target_floor=q['H']>0 and q['H']**2>=c['floors']['oriented_unit_stress']**2*n2,
                same_unchanged_tilt_new_margin=state['c2']*q['J']**2<=k*(1-gammaB)**2*q['H']**2)))
    return dict(auxiliary_rational_shears_not_actual_profile=True,rows=controls,legacy_comparisons=comparisons,
        checks=dict(all_ninety_six_exact_perturbation_controls=all(all(r['checks'].values()) for r in controls),
            normalized_legacy_constants_and_strict_margin_improvement=all(all(r['checks'].values()) for r in comparisons)))

def obstruction_review():
    state=gap().shear(1,F(-82,25),0,F(82,25));u=F(3,4);k=u*u/(1+u*u)
    data=constants(F(1,2),2,2,5,4,F(1,2),u);rows=[]
    for n in (5,10,20,60):
        t=F(1,2**n);before=reference_target(t);after=reference_target(-t)
        a=gap().cone(state,before,u);b=gap().cone(state,after,u);R=state['c2']*a['J']**2/a['H']**2
        slack=k-R;c=certificate(data,slack);error=4*t/(1+t*t)
        rows.append(dict(n=n,slack=str(slack),cone_radius=str(c['cone_radius']),relative_error=str(error),checks=dict(
            original_tilt_admissible_with_positive_deficit=slack>0,
            same_fixed_positive_spectral_gap=state['D']==F(2624,625) and state['c2']==F(16,25),
            exact_unit_directions=sum(v*v for v in before)==sum(v*v for v in after)==1,
            exact_small_relative_direction_error=sum((a-b)**2 for a,b in zip(before,after))==error*error,
            arbitrary_small_direction_error_can_cross_unquantified_margin=a['phi']>0>b['phi'],
            normalized_radius_linear_in_unquantified_deficit=c['cone_radius']/slack==data['inputs']['omega']**2/(2*data['Lpsi']))))
    coarse=module('leading_cone_margin','cone_bounds.py').tilt_gate(F(1,2**40),u*u)
    t=reference_target(F(1,32));c0=gap().cone(state,t,u);R=state['c2']*c0['J']**2/c0['H']**2
    return dict(auxiliary_obstruction_not_counterexample_to_actual_profile=True,rows=rows,checks=dict(
        all_four_fixed_tilt_thin_margin_controls=all(all(r['checks'].values()) for r in rows),
        strictly_decreasing_slack_and_cone_radius=all(F(rows[i+1]['slack'])<F(rows[i]['slack']) and F(rows[i+1]['cone_radius'])<F(rows[i]['cone_radius']) for i in range(3)),
        coarse_gate_failure_does_not_remove_original_admissibility=not coarse['accepted'] and tilt_slack(u,R)['admissible_by_supplied_bound']))

def symbolic_review():
    import sympy as s
    u,H,J,c2,R,delta,w,L,e=s.symbols('u H J c2 R delta w L e',positive=True)
    k=u*u/(1+u*u);gamma=w*w*delta/(9*k);lead=(1+u*u)*H*H*(k-R)
    z=s.symbols('z',real=True);co=(1-z*z)/(1+z*z);si=2*z/(1+z*z)
    hz=(4*co+3*si)/5;jz=(3*co-4*si)/5
    return dict(checks=dict(
        exact_normalized_cone_polynomial=s.simplify((u*u*H*H-(1+u*u)*c2*J*J)/(1+u*u)-(k*H*H-c2*J*J))==0,
        exact_leading_raw_supremum_deficit=s.simplify(lead-(u*u*H*H-(1+u*u)*R*H*H))==0,
        inverse_squared_tilt_coarse_margin_is_raw_slack_over_kappa=s.simplify(1-R/k-(k-R)/k)==0,
        direct_margin_uses_nine_from_unit_relative_target_bound=s.simplify(w*w*delta/(4*k*(1+s.Rational(1,2))**2)-gamma)==0,
        freezing_condition_has_seventy_two=s.simplify(8*s.Symbol('Cfreeze')/gamma-72*k*s.Symbol('Cfreeze')/(w*w*delta))==0,
        required_background_deficit_exact=s.simplify(e-w*w*(2*L*e/(w*w))/(2*L))==0,
        direct_deficit_improves_gamma_half_route=s.simplify(delta-k*(delta/(2*k))*(2-delta/(2*k))-delta**2/(4*k))==0,
        thin_margin_exact_unit_target=s.simplify(hz*hz+jz*jz-1)==0,
        thin_margin_exact_cone_factorization=s.simplify(s.Rational(9,25)*hz*hz-s.Rational(16,25)*jz*jz-2*z*(24*(1-z*z)-14*z)/(25*(1+z*z)**2))==0))

def join_review():
    g=gates();p=g.fixture();rows=g.affine_gate(**p)['rows'];data=constants(F(1,2),2,2,5,4,F(1,2),F(3,4))
    slack=F(1,100);Cfreeze=F(7);w=later_witness(data,rows,32,1216,slack,p['h'],Cfreeze)
    ell=w['sufficient_joint_ell'];j=auxiliary_shape_join(data,rows,slack,p['h'],ell,Cfreeze)
    test=required_deficit(data,F(1,10**20),1000,2)
    erased=dict(p,leading=(*p['leading'][:5],F(0),F(0),F(0)))
    return dict(auxiliary_full_interface_not_actual_scale=True,complete_rows=serial(rows),joint_witness=serial(w),
        auxiliary_shape_join=serial(j),checks=dict(
            same_full_pressure_and_viscosity_rows_retained=rows['stress'][0]>g.affine_gate(**erased)['rows']['stress'][0],
            whole_chart_and_independent_freezing_witness=all(w['checks'].values()),
            sharp_shape_program_receives_new_radius_but_does_not_select_old_cutoff=j['freezing_gate'] and j['auxiliary_shape_program']['possible'],
            freezing_gate_cannot_be_erased_even_when_background_error_tiny=test['geometry_gate'] and test['freezing_deficit']>data['kappa'],
            later_witness_keeps_original_h_and_plus_two=w['background_witness']['sufficient_ell']==2+g.ceiling(F(w['background_witness']['dyadic_exponent'])/(2*p['h']))))

def review():
    p,o,s,j=perturbation_review(),obstruction_review(),symbolic_review(),join_review()
    return dict(perturbations=p,obstruction=o,symbolic=s,join=j,checks=dict(
        exact_normalized_deficit_radius_and_legacy_comparison=all(p['checks'].values()),
        original_admissibility_not_a_quantitative_fixed_scale_slack=all(o['checks'].values()),
        symbolic_normalization_margin_and_freezing_cost=all(s['checks'].values()),
        complete_corrected_width_shape_and_freezing_scale_join=all(j['checks'].values())))
