"""High-precision original-phase controls with exact integer rounding.

These are manufactured accepted-input examples, not the selected giant
background, its actual epsilon schedule or continuum supremum certificates.
"""
from fractions import Fraction as F
from functools import lru_cache
from mpmath import mp
from frame_budget import inputs,budgets


def nearest_nonzero(x,k):
    x=F(x)
    if type(k) is not int or k<1:raise ValueError('Positive integer carrier required')
    z=k*x;n=z.numerator//z.denominator
    candidates={n-1,n,n+1,n+2,-1,1}-{0}
    selected=min(candidates,key=lambda q:(abs(F(q)-z),abs(q),q))
    return F(selected,k)


def phase_data(S,sign,fraction):
    if type(S) is not int or S<1 or sign not in (-1,1) or not 0<=F(fraction)<=1:raise ValueError('Finite positive S, sign and pulse fraction required')
    p=inputs();b=budgets(p);r0=F(1,6057072);L=2*r0*S;k=2*S*S;eps=F(1,k*k);v=F(fraction)*L
    R0=F(2);R=R0+p['Dbox']/S**3;Bs=F(16,25);u=p['u'];K=(F(3,5),F(-4,5));g0=(F(-8,5),F(-6,5))
    pt=R0*Bs*(K[0]-sign*u*g0[0]/(4*L));pz=Bs*(K[1]-sign*u*g0[1]/(4*L));angular=nearest_nonzero(pt,k)
    F0=F(1);Fv=F0+p['Db']*eps*eps;FR0=F(-4,5);GR0=F(-6,5)
    FR=FR0+p['Db']*eps*eps;GR=GR0-p['Db']*eps*eps;FZ=p['Mz']/2;GZ=-p['Mz']/3
    WR=angular*FR+pz*GR;WZ=angular*FZ+pz*GZ
    n=(sign*Bs*u/2-v*WR,angular/R,pz-eps*v*WZ);np=(-WR,F(0),-eps*WZ)
    return dict(parameters=p,budget=b,S=S,L=L,k=k,epsilon=eps,v=v,R=R,R0=R0,Bs=Bs,u=u,
        K=K,g0=g0,ptilde=pt,p=angular,pz=pz,F0=F0,F=Fv,FR=FR,GR=GR,FZ=FZ,GZ=GZ,
        n=n,nprime=np,sref=sign*u*(F(1,2)+F(fraction)),lambda0=F(4,5),c0=F(-1,2))


@lru_cache(None)
def phase_review():
    rows=[];good=True;frame=True;reference_frame=True;damping=True;column=True;midpoint=True;identity=True
    max_identity=mp.mpf(0)
    with mp.workdps(80):
        def m(x):
            x=F(x);return mp.mpf(x.numerator)/x.denominator
        for S in (10**14,10**16):
            for sign in (-1,1):
                for fraction in (F(0),F(1,2),F(1)):
                    d=phase_data(S,sign,fraction);b=d['budget'];L=m(d['L']);u=m(d['u']);v=m(d['v'])
                    n=mp.matrix(list(map(m,d['n'])));npv=mp.matrix(list(map(m,d['nprime'])));t=mp.sqrt(n[1]**2+n[2]**2)
                    a=mp.matrix([n[1]/t,n[2]/t]);na=mp.matrix([a[1],-a[0]]);s=n[0]/t;w=mp.sqrt(1+s*s)
                    tp=a[0]*npv[1]+a[1]*npv[2];sp=npv[0]/t-s*tp/t
                    ap=mp.matrix([(npv[i+1]-a[i]*tp)/t for i in range(2)]);om=(ap.T*na)[0]
                    c0=m(d['c0']);Fv=m(d['F']);g=mp.matrix([m(d['R']*d['FR']),m(d['GR'])])
                    K=mp.matrix([[0,-2*Fv,0],[2*Fv+g[0],0,0],[g[1],0,0]])
                    An=-K+n*(n.T*K-npv.T)/(n.T*n)[0]
                    U=mp.matrix([[1,0],[-s*a[0],na[0]],[-s*a[1],na[1]]])
                    Up=mp.matrix([[0,0],[-sp*a[0]-s*om*na[0],-om*a[0]],[-sp*a[1]-s*om*na[1],-om*a[1]]])
                    M=mp.matrix([[1,1],[c0*w,-c0*w]]);Mi=M**-1
                    Mp=mp.matrix([[0,0],[c0*s*sp/w,-c0*s*sp/w]])
                    B=U*M;Bp=Up*M+U*Mp;left=Mi*mp.matrix([[1,0,0],[0,na[0],na[1]]])
                    actual=left*(An*B-Bp)
                    H=mp.matrix([[s*((a.T*g)[0]-sp)/(1+s*s),(2*Fv*na[0]-s*om)/(1+s*s)],
                        [-(2*Fv*na[0]+(na.T*g)[0])+s*om,0]])
                    error=mp.norm(actual-(Mi*H*M-Mi*Mp));max_identity=max(max_identity,error);identity &= error<mp.mpf('1e-60')
                    sr=m(d['sref']);lam0=m(d['lambda0']);lam=lam0/mp.sqrt(1+sr*sr)
                    target=mp.diag([lam,-lam]);E=mp.norm(actual-target)
                    wr=mp.sqrt(1+sr*sr);srp=sign*u/L
                    Mr=mp.matrix([[1,1],[c0*wr,-c0*wr]])
                    Mrp=mp.matrix([[0,0],[c0*sr*srp/wr,-c0*sr*srp/wr]])
                    Br=U*Mr;Brp=Up*Mr+U*Mrp
                    leftr=Mr**-1*mp.matrix([[1,0,0],[0,na[0],na[1]]])
                    Er=mp.norm(leftr*(An*Br-Brp)-target)
                    reference_frame &= L*Er<=m(b['CEref'])*(1+mp.mpf('1e-50'))
                    # Frobenius is a stronger control than the claimed operator norm.
                    frame &= L*E<=m(b['CE'])*(1+mp.mpf('1e-50'))
                    damp=m(d['epsilon'])*d['k']**2*(n.T*n)[0]
                    dr=lam0*(1+sr*sr)/(1+u*u)**mp.mpf('1.5')
                    damping &= L*abs(damp-dr)<=m(b['Cd'])*(1+mp.mpf('1e-50'))
                    kr=mp.matrix(list(map(m,d['K'])));nr=mp.matrix([n[1],n[2]])
                    normalerr=mp.norm(nr-m(d['Bs'])*kr)
                    good &= normalerr*L<=m(b['Dt']) and abs(n[0]-m(d['Bs'])*sr)*L<=m(b['Dr'])
                    good &= t>=m(b['beta']) and d['L']>=b['normal_acceptance'] and d['L']>=b['Riccati_acceptance']
                    Nref=mp.matrix([kr[1],-kr[0]])
                    worstcol=mp.mpf(0)
                    for qsign in (-1,1):
                        q=qsign*m(b['CE'])/(2*m(b['mu'])*L)
                        direction=-s*a+c0*w*(1-q)/(1+q)*na
                        ideal=-sr*kr+c0*mp.sqrt(1+sr*sr)*Nref
                        diff=mp.norm(direction-ideal);column &= diff*L<=m(b['D0'])
                        worstcol=max(worstcol,diff*L/m(b['D0']))
                    instant=-sr*kr+c0*mp.sqrt(1+sr*sr)*Nref
                    center=-sign*u*kr+c0*mp.sqrt(1+u*u)*Nref
                    midpoint &= mp.norm(instant-center)<=m(b['D1'])*abs(v-L/2)/L+mp.mpf('1e-60')
                    rows.append(dict(S=S,sign=sign,pulse_fraction=str(fraction),k=d['k'],
                        normal_error_budget_ratio=float(normalerr*L/m(b['Dt'])),
                        frame_error_budget_ratio=float(E*L/m(b['CE'])),reference_frame_error_budget_ratio=float(Er*L/m(b['CEref'])),
                        damping_error_budget_ratio=float(abs(damp-dr)*L/m(b['Cd'])),
                        column_error_budget_ratio=float(worstcol)))
    return dict(rows=rows,maximum_direct_reduced_identity_error=float(max_identity),
        checks=dict(original_full_phase_primitive_and_floor_budgets=bool(good),
            direct_and_reduced_frame_generator_match=bool(identity),explicit_operator_error_budget_in_phase_controls=bool(frame),
            manuscript_reference_mixing_operator_budget_in_full_phase_controls=bool(reference_frame),
            explicit_damping_budget_in_phase_controls=bool(damping),both_frame_coordinates_obey_column_budget=bool(column),
            midpoint_direction_Lipschitz_budget=bool(midpoint),
            actual_selected_schedule_and_giant_fields_not_evaluated=True))
