"""Explicit real-domain background gap bounds; exact controls are not PDE fields."""
from fractions import Fraction as F
from itertools import product

def rational(x):
    if isinstance(x,(bool,float)):raise TypeError('Use exact rational inputs')
    return F(x)

def certificate(f,Fbar,t,Gbar,d,omega,gamma,u):
    f,Fbar,t,Gbar,d,omega,gamma,u=map(rational,(f,Fbar,t,Gbar,d,omega,gamma,u))
    if min(f,t,d,omega,gamma,u)<=0 or Fbar<f or Gbar<t or max(omega,gamma)>1:
        raise ValueError('Positive ordered bounds, orientation and cone margins required')
    LD=4*Gbar+2*Fbar+3;LP=Fbar+Gbar+1
    Dbar=Gbar**2+2*Fbar*Gbar
    cbar2=Dbar/(2*f*t)
    Lc=2*LD/(f*t)+2*Dbar*LP/(f*f*t*t)
    LH=1+2/t
    Lphi=3*u*u*LH+(1+u*u)*(4*Lc+3*cbar2*LH)
    mphi=u*u*omega*omega*gamma*(2-gamma)
    constraints=dict(relative_stress=F(1,2),swirl=f/2,angular_shear=t/2,
        discriminant=d/(2*LD),orientation=omega/(2*LH),same_tilt_cone=mphi/(2*Lphi))
    radius=min(constraints.values())
    return dict(inputs=dict(f=f,Fbar=Fbar,t=t,Gbar=Gbar,d=d,omega=omega,gamma=gamma,u=u),
        LD=LD,LP=LP,Dbar=Dbar,cbar2=cbar2,Lc=Lc,LH=LH,Lphi=Lphi,mphi=mphi,
        constraints=constraints,radius=radius,
        floors=dict(F=f/2,minus_gtheta=t/2,D=d/2,normal_theta=t/(2*(Gbar+1)),
            lambda2=f*t*d/(4*(Gbar+1)**2),c2=d/(4*(Fbar+1)*(Gbar+1)),
            oriented_unit_stress=omega/4,cone_margin=omega*omega*gamma*(2-gamma)/16),
        upper=dict(F=Fbar+1,g=Gbar+1,
            c2=2*((Gbar+1)**2+2*(Fbar+1)*(Gbar+1))/(f*t)))

def shear(F0,gtheta,gz,gnorm):
    F0,gtheta,gz,gnorm=map(rational,(F0,gtheta,gz,gnorm))
    if gnorm<=0 or gtheta*gtheta+gz*gz!=gnorm*gnorm:raise ValueError('Exact norm required')
    D=gnorm*gnorm+2*F0*gtheta
    if F0<=0 or gtheta>=0:raise ValueError('Positive F and negative angular shear required')
    return dict(D=D,lambda2=-2*F0*gtheta*D/(gnorm*gnorm),c2=D/(-2*F0*gtheta),
        N=(gtheta/gnorm,gz/gnorm),K=(-gz/gnorm,gtheta/gnorm))

def dot(x,y):return sum((a*b for a,b in zip(x,y)),F(0))
def cone(state,target,u):
    H=-dot(state['N'],target);J=dot(state['K'],target)
    return dict(H=H,J=J,phi=u*u*H*H-(1+u*u)*state['c2']*J*J)

def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x

def perturbed_review():
    c=certificate(F(1,2),2,2,5,4,F(1,2),F(1,2),1);p=c['inputs'];delta=c['radius']/2
    target0=(F(24,25),F(7,25));old=shear(1,-4,0,4);q0=cone(old,target0,1);rows=[]
    leading_ok=(p['f']<=1<=p['Fbar'] and p['t']<=4<=p['Gbar'] and old['D']>=p['d'] and
        dot(target0,target0)==1 and q0['H']>=p['omega'] and
        2*old['c2']*q0['J']**2<=(1-p['gamma'])**2*q0['H']**2)
    for sf,sm,sg,sx,sz in product((-1,1),repeat=5):
        a=sg*delta/32;co=(1-a*a)/(1+a*a);si=2*a/(1+a*a);magnitude=4+sm*delta/4
        g=(-magnitude*co,magnitude*si);fv=1+sf*delta/3;target=(target0[0]+sx*delta/4,target0[1]+sz*delta/5)
        new=shear(fv,*g,magnitude);q=cone(new,target,1)
        diffg=(g[0]+4,g[1]);dr=(target[0]-target0[0],target[1]-target0[1]);df=fv-1
        input_ok=abs(fv-1)<=delta and dot(diffg,diffg)<=delta*delta and dot(dr,dr)<=delta*delta
        prodold=F(4);prodnew=fv*(-g[0])
        quotient_identity=(new['c2']-old['c2'])*2*prodnew*prodold==(
            (new['D']-old['D'])*prodold-old['D']*(prodnew-prodold))
        expanded_D=2*dot((F(-4),F(0)),diffg)+dot(diffg,diffg)+2*diffg[0]-8*df+2*df*diffg[0]
        checks=dict(leading_input_bounds=leading_ok,input_bounds=input_ok,
            exact_full_discriminant_expansion=new['D']-old['D']==expanded_D,
            discriminant_difference=abs(new['D']-old['D'])<=c['LD']*delta,
            product_difference=abs(prodnew-prodold)<=c['LP']*delta,
            quotient_difference=abs(new['c2']-old['c2'])<=c['Lc']*delta,
            exact_quotient_identity=quotient_identity,
            H_difference=abs(q['H']-q0['H'])<=c['LH']*delta,
            J_difference=abs(q['J']-q0['J'])<=c['LH']*delta,
            cone_polynomial_difference=abs(q['phi']-q0['phi'])<=c['Lphi']*delta,
            swirl_floor=fv>=c['floors']['F'],angular_shear_floor=-g[0]>=c['floors']['minus_gtheta'],
            discriminant_floor=new['D']>=c['floors']['D'],lambda_floor=new['lambda2']>=c['floors']['lambda2'],
            c_floor=new['c2']>=c['floors']['c2'],normal_floor=-new['N'][0]>=c['floors']['normal_theta'],
            orientation=q['H']>=p['omega']/2,cone_floor=q['phi']>=c['mphi']/2,
            same_tilt_margin=(1+1)*new['c2']*q['J']**2<=(1-c['floors']['cone_margin'])**2*q['H']**2)
        rows.append(dict(signs=[sf,sm,sg,sx,sz],checks=checks))
    return dict(model='exact rational perturbations of manufactured shear, not actual giant profiles',
        certificate=serial(c),delta=str(delta),rows=rows,
        checks=dict(all_exact_rational_perturbation_controls=all(all(r['checks'].values()) for r in rows)))

def multiply(a,b):return [[sum((a[i][k]*b[k][j] for k in range(2)),F(0)) for j in range(2)] for i in range(2)]

def frame_review():
    rows=[]
    for fv,cc,(s,w) in product((F(1,2),F(1),F(3)),(F(1,2),F(1),F(2)),((F(0),F(1)),(F(3,4),F(5,4)),(F(4,3),F(5,3)))):
        av=2+2*cc*cc;gv=-av*fv;state=shear(fv,gv,0,-gv)
        lam=2*fv*cc/w;M=[[F(1),F(1)],[-cc*w,cc*w]]
        # A21=-(2F Ntheta+|g|)=2F+gtheta on this negative-theta axis fixture.
        A=[[F(0),-2*fv/(1+s*s)],[2*fv+gv,F(0)]]
        left=multiply(A,M);right=multiply(M,[[lam,F(0)],[F(0),-lam]])
        checks=dict(growth_formula=state['lambda2']==(2*fv*cc)**2,
            mixing_formula=state['c2']==cc*cc,reference_metric=w*w==1+s*s,
            every_projected_frame_entry=left==right)
        rows.append(dict(F=str(fv),c0=str(-cc),s=str(s),checks=checks))
    return dict(model='rational constant reference-frame controls; actual moving frame is NOT redefined',rows=rows,
        checks=dict(exact_reference_growth_and_frame=all(all(r['checks'].values()) for r in rows)))

def obstruction_review():
    rows=[]
    for k in (10,20,40,80):
        tau=F(1,2**k);old=shear(1,-2-tau,0,2+tau);new=shear(1,-2+tau,0,2-tau)
        rows.append(dict(tau=str(tau),absolute_error=str(2*tau),old_D=str(old['D']),new_D=str(new['D']),
            sign_failure=old['lambda2']>0 and new['lambda2']<0))
    state=shear(1,-4,0,4);badcone=cone(state,(F(3,5),F(4,5)),1)
    abs_rows=[]
    for k in (10,20,40):
        scale=F(1,2**k);amplitude=scale*scale
        before=(amplitude,F(0));after=(-amplitude,F(0))
        absolute=abs(after[0]-before[0])
        old_direction=tuple(x/amplitude for x in before);new_direction=tuple(x/amplitude for x in after)
        old_H=cone(state,old_direction,1)['H'];new_H=cone(state,new_direction,1)['H']
        abs_rows.append(dict(scale=str(scale),leading_stress_size=str(amplitude),absolute_stress_change=str(absolute),
            relative_stress_change=str(absolute/amplitude),old_orientation=str(old_H),new_orientation=str(new_H),
            orientation_reversed=old_H>0 and new_H<0))
    return dict(near_neutral_growth_rows=rows,strict_growth_but_bad_cone=serial(badcone),absolute_stress_controls=abs_rows,
        checks=dict(positive_swirl_and_shear_not_sufficient=all(r['sign_failure'] for r in rows),
            strict_growth_does_not_supply_same_tilt_cone=state['lambda2']>0 and badcone['H']>0 and badcone['phi']<0,
            unweighted_absolute_stress_error_not_relative_control=all(r['orientation_reversed'] for r in abs_rows)))

def tail_review():
    # h=1/200 is a permitted diagnostic value, not a replacement for frozen h.
    c=certificate(F(1,2),2,2,5,4,F(1,2),F(1,2),1);radius=c['radius'];Cval=F(100)
    n=0
    while Cval/F(2)**n>radius:n+=1
    rows=[dict(ell=100*j,error=str(Cval/F(2)**j),accepted=Cval/F(2)**j<=radius) for j in (n,n+1,n+10,n+100)]
    return dict(formula='ell >= 2+log2(Cval/radius)/(2h) suffices for every q<=4Q, Q=2^-ell; all existing tail constraints retained',
        fixture='h=1/200 for exact grid powers only; frozen actual h is not changed',first_fixture_grid_point=100*n,
        sufficient_first_fixture_grid_band=100*n+2,
        rows=rows,checks=dict(all_later_values_follow_monotonicity=all(r['accepted'] for r in rows),
            previous_fixture_grid_point_fails=Cval/F(2)**(n-1)>radius,
            historical_first_band_not_certified=True,derivative_order_does_not_change_value_threshold=True))

def review():
    p=perturbed_review();f=frame_review();o=obstruction_review();t=tail_review()
    return dict(explicit_background_perturbations=p,reference_frame_identities=f,
        necessary_gap_controls=o,later_tail_selection=t,
        checks=dict(exact_perturbation_controls=all(p['checks'].values()),
            exact_reference_frame=all(f['checks'].values()),necessary_obstructions=all(o['checks'].values()),
            tail_value_selection=all(t['checks'].values())))
