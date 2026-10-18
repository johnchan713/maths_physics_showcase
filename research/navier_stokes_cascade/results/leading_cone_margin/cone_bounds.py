"""Exact global leading direction bounds and a fixed-tilt sufficient gate."""
from fractions import Fraction as F
from itertools import product
import importlib.util
from pathlib import Path

def rat(x):
    if isinstance(x,(bool,float)):raise TypeError('Exact rational inputs required')
    return F(x)
def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x

def margins(G):
    G=rat(G)
    if G<2:raise ValueError('Retained hierarchy requires G>=2')
    return dict(omega=1/(640*G**13),epsilon=1/(2**18*G**22))

def coordinates(a,b,p1,p2):
    a,b,p1,p2=map(rat,(a,b,p1,p2))
    if a<=0:raise ValueError('Positive a required')
    t=-b/a;v=a+b*b/a;D=p1+t*p2-v;J=p2-t*p1
    if D<=0 or v<2:raise ValueError('Oriented D>0 and v>=2 required')
    return dict(t=t,v=v,D=D,J=J,fourth=2*D*D-(v-2)*J*J,
        orientation2=D*D/(D*D+J*J),raw_ratio2=(v-2)*J*J/(2*D*D))

def tilt_gate(epsilon,u2):
    """No tilt is chosen or changed. A failed sufficient test is not cone failure."""
    epsilon,u2=map(rat,(epsilon,u2))
    if not 0<epsilon<1 or u2<=0:raise ValueError('0<epsilon<1 and the positive EXISTING u^2 required')
    upper=(1+1/u2)*(1-epsilon);delta=1-upper;accepted=delta>0
    return dict(squared_ratio_upper=upper,delta=delta,accepted=accepted,
        gamma=delta/2 if accepted else None,sharp_sufficient_u2_threshold=(1-epsilon)/epsilon,
        padded_u2_threshold=2/epsilon,padded_accepted=u2>=2/epsilon,
        padded_gamma=epsilon/4 if u2>=2/epsilon else None)

def chart_inputs(fmin,fmax,alpha,abar,bbar,nu):
    """All quarter-to-four band factors charged, with 1/5<4^-1-h<1."""
    fmin,fmax,alpha,abar,bbar,nu=map(rat,(fmin,fmax,alpha,abar,bbar,nu))
    if min(fmin,alpha,nu)<=0 or fmax<fmin or abar<alpha or bbar<0:
        raise ValueError('Positive ordered leading floors and bounds required')
    return dict(f=fmin/5,Fbar=5*fmax,t=alpha*fmin/5,
        Gbar=5*fmax*(abar+bbar),d=alpha*nu*fmin*fmin/25)

def compact_review():
    R=F(16);k=F(1,100);omega=k/(10*R**3);epsilon=k/(64*R**6);rows=[]
    for a,b,dt,dz in product(tuple(map(F,(3,4,5))),(F(-1,2),F(0),F(1,2)),
                             tuple(map(F,(1,2,3))),(F(-1,10),F(0),F(1,10))):
        p1=a+dt;p2=-b+dz;r=coordinates(a,b,p1,p2)
        checks=dict(actual_coordinate_budget=a>=1/R and max(map(abs,(a,b,p1,p2)))<=R,
            four_gap_inputs=min(a,r['v']-2,r['D'],r['fourth'])>=k/2,
            exact_raw_deficit_identity=1-r['raw_ratio2']==r['fourth']/(2*r['D']**2),
            D_upper=r['D']<=4*R**3,J_upper=abs(r['J'])<=2*R**3,
            oriented_unit_floor=r['orientation2']>=omega*omega,
            raw_cone_deficit=r['raw_ratio2']<=1-epsilon)
        rows.append(dict(a=str(a),b=str(b),theta_excess=str(dt),axial_excess=str(dz),checks=checks))
    return dict(model='81 exact local coordinate controls, not selected-profile samples',rows=rows,
        compact_omega='kappa/(10R^3)=1/(640G^13)',compact_epsilon='kappa/(64R^6)=1/(2^18 G^22)',
        checks=dict(all_exact_coordinate_and_margin_controls=all(all(r['checks'].values()) for r in rows),
            exact_hierarchy_constants=10*4**3==640 and 64*4**6==2**18))

def exterior_review():
    rows=[];omega=F(1,64000);epsilon=F(1,490000000)
    for pc,p,v,j in product(tuple(map(F,(1,2,100,14000))),(F(10**12),F(2*10**12)),
                            tuple(map(F,(3,10,100))),(F(-1,1000),F(0),F(1,1000))):
        D=pc-v/p;fourth=2*D*D-(v-2)*j*j;H2=D*D/(D*D+j*j);q2=(v-2)*j*j/(2*D*D)
        checks=dict(finite_p_D_lower=D>=F(1,4),D_upper=D<=14000,J_upper=abs(j)<=2000,
            complete_fourth_margin=fourth>=F(4,5),orientation_floor=H2>=omega*omega,
            exact_deficit_identity=1-q2==fourth/(2*D*D),raw_deficit=q2<=1-epsilon)
        rows.append(dict(Pc_over_p=str(pc),p=str(p),v=str(v),J_over_p=str(j),checks=checks))
    return dict(model='72 exact uses of retained compensated-exterior inequalities, not a new heat proof',rows=rows,
        bounds=dict(omega=str(omega),epsilon=str(epsilon)),
        checks=dict(all_finite_radius_direction_controls=all(all(r['checks'].values()) for r in rows),
            exact_finite_p_lower=F(29,100)-F(120,10**12)>F(1,4),
            actual_w_and_shear_ratio_budget=F(1801)+F(15,2)<2000,
            safe_norm_denominator=14000**2+2000**2<16000**2,
            exact_fourth_to_raw_deficit=F(4,5)/(2*14000**2)==epsilon))

def global_review():
    G=F(2);r=margins(G);outer=[]
    for h in (F(1,100),F(1,200),F(1,1000)):
        raw=16*h**9;H2=1/(1+16*h**8)
        outer.append(dict(h=str(h),checks=dict(terminal_raw_below_half=raw<F(1,2),
            terminal_orientation_above_half=H2>F(1,4),global_epsilon_retained=raw<1-r['epsilon'],
            global_orientation_retained=H2>r['omega']**2)))
    return dict(domain='closed [Xa,Xb] with limiting unit direction at both flat edges',
        actual_hierarchy='G>C^4096, C>=2; actual G and C are not evaluated',
        omega='1/(640G^13)',epsilon='1/(2^18 G^22)',inner_orientation_comparison='180 G^13/C^4 >180 C^53244',
        terminal_controls=outer,
        checks=dict(inner_orientation_power_dominance=13*4096-4==53244>0,
            global_epsilon_below_inner_half=r['epsilon']<F(1,2),
            global_omega_below_exterior= r['omega']<F(1,64000),
            global_epsilon_below_exterior=r['epsilon']<F(1,490000000),
            all_terminal_direction_bounds=all(all(x['checks'].values()) for x in outer)))

def tilt_review():
    rows=[]
    for e in (F(1,4),F(1,16),F(1,2**40)):
        threshold=(1-e)/e
        for U in (threshold/2,threshold,2*threshold,2/e,4/e):
            r=tilt_gate(e,U)
            checks=dict(exact_threshold=r['accepted']==(U>threshold),
                strict_boundary_not_accepted=(U!=threshold or not r['accepted']),
                safe_linear_sqrt_margin=(not r['accepted'] or (1-r['gamma'])**2>=r['squared_ratio_upper']),
                padded_margin=(not r['padded_accepted'] or (1-e/4)**2>=r['squared_ratio_upper']))
            rows.append(dict(epsilon=str(e),u_squared=str(U),certificate=serial(r),checks=checks))
    e=margins(2)['epsilon'];coarse=tilt_gate(e,1);actual_raw2=F(1,4);actual_tilted2=2*actual_raw2
    sharp=[dict(u2=str(U),squared_tilted_ratio=str((1+1/U)*F(9,16))) for U in (F(1),F(9,7),F(2))]
    return dict(rows=rows,
        missed_sufficient_control=dict(model='abstract exact cone data, not an actual chosen profile or tilt',
            epsilon=str(e),existing_u2='1',actual_raw2=str(actual_raw2),actual_tilted2=str(actual_tilted2),
            coarse_certificate=serial(coarse)),sharp_boundary_control=sharp,
        checks=dict(all_fixed_tilt_controls=all(all(x['checks'].values()) for x in rows),
            coarse_gate_failure_is_not_actual_cone_failure=not coarse['accepted'] and actual_raw2<=1-e and actual_tilted2<1,
            equality_cannot_supply_strict_margin=F(sharp[1]['squared_tilted_ratio'])==1,
            both_sides_of_sharp_boundary=F(sharp[0]['squared_tilted_ratio'])>1>F(sharp[2]['squared_tilted_ratio'])))

def parent_review():
    p=Path(__file__).resolve().parent.parent/'background_gap_certificate/gap_bounds.py'
    s=importlib.util.spec_from_file_location('leading_cone_parent_gap',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    r=margins(2);u=F(2**21);tilt=tilt_gate(r['epsilon'],u*u)
    chart=chart_inputs(F(1,10),2,F(1,4),4,1,F(1,10))
    cert=m.certificate(**chart,omega=r['omega'],gamma=r['epsilon']/4,u=u)
    return dict(model='rational envelope substitution into the original gap map; NOT the selected profile or tilt',
        direction=serial(r),illustrative_existing_u=str(u),tilt=serial(tilt),chart=serial(chart),
        certificate=serial(cert),
        checks=dict(padded_existing_tilt_certified=tilt['padded_accepted'],
            exact_chart_swirl_floor=chart['f']==F(1,50),
            exact_squared_chart_discriminant_floor=chart['d']==F(1,100000),
            unchanged_parent_radius_is_positive=cert['radius']>0,
            parent_gamma_uses_supplied_global_deficit=cert['inputs']['gamma']==r['epsilon']/4))

def review():
    c=compact_review();e=exterior_review();g=global_review();t=tilt_review();p=parent_review()
    return dict(compact_direction_margins=c,compensated_exterior_direction_margins=e,
        closed_annulus_direction_join=g,existing_tilt_gate=t,scaled_original_gap_interface=p,
        checks=dict(compact_exact_margins=all(c['checks'].values()),exterior_finite_radius_margins=all(e['checks'].values()),
            closed_annulus_join=all(g['checks'].values()),fixed_tilt_and_obstruction_controls=all(t['checks'].values()),
            original_scaled_gap_interface=all(p['checks'].values())))
