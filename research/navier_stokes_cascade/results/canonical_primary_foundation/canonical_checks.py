"""Independent physical ODE/covariance diagnostics, not continuum certificates."""
from functools import lru_cache
import importlib.util
from math import asinh,exp,sqrt,log
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import quad,solve_ivp


@lru_cache(None)
def phase_module():
    p=Path(__file__).resolve().parent.parent/'frame_column_foundation'
    old=sys.path[:];sys.path.insert(0,str(p))
    try:
        spec=importlib.util.spec_from_file_location('canonical_original_phase',p/'phase_checks.py')
        m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
    finally:sys.path[:]=old


def cutoff(t):
    d=abs(t-.5)
    if d<=.2:return 1.
    if d>=1/3:return 0.
    z=(d-.2)/(1/3-.2);a,b=exp(-1/z),exp(-1/(1-z));return b/(a+b)


def logP(v,L,u,lam0):
    a=.5+v/L
    return lam0*L*((asinh(u*a)-asinh(u))/u-((a-1)+u*u*(a**3-1)/3)/(1+u*u)**1.5)


def coefficients(data,v):
    L=float(data['L']);n0=np.array(data['n'],dtype=float);npv=np.array(data['nprime'],dtype=float)
    n=n0+v*npv;g=np.array([data['R']*data['FR'],data['GR']],dtype=float);fv=float(data['F'])
    K=np.array([[0,-2*fv,0],[2*fv+g[0],0,0],[g[1],0,0]])
    An=-K+np.outer(n,n@K-npv)/(n@n)
    mag=np.linalg.norm(n[1:]);a=n[1:]/mag;b=np.array([a[1],-a[0]]);s=n[0]/mag
    magp=a@npv[1:];sp=npv[0]/mag-s*magp/mag;ap=(npv[1:]-a*magp)/mag
    om=ap@b;U=np.array([[1.,0.],[-s*a[0],b[0]],[-s*a[1],b[1]]])
    Up=np.array([[0.,0.],[-sp*a[0]-s*om*b[0],-om*a[0]],[-sp*a[1]-s*om*b[1],-om*a[1]]])
    u=float(data['u']);sign=1 if data['sref']>0 else -1;sr=sign*u*(.5+v/L)
    wr=sqrt(1+sr*sr);wa=sqrt(1+s*s);c=float(data['c0']);lam0=float(data['lambda0'])
    Mr=np.array([[1.,1.],[c*wr,-c*wr]]);Ma=np.array([[1.,1.],[c*wa,-c*wa]])
    Mrp=np.array([[0.,0.],[c*sr*sign*u/(L*wr),-c*sr*sign*u/(L*wr)]])
    Br=U@Mr;Ba=U@Ma;Brp=Up@Mr+U@Mrp
    left=np.linalg.solve(Mr,np.array([[1.,0.,0.],[0.,b[0],b[1]]]))
    lam=lam0/wr;dr=lam0*(1+sr*sr)/(1+u*u)**1.5
    d=float(data['epsilon'])*data['k']**2*(n@n);rate=lam-dr
    normalized=An-(d+rate)*np.eye(3)
    reduced=left@(An@Br-Brp)-(d+rate)*np.eye(2)
    return dict(n=n,np=npv,K=K,An=An,Ba=Ba,Br=Br,Brp=Brp,left=left,normalized=normalized,reduced=reduced,d=d,rate=rate)


@lru_cache(None)
def physical_review():
    rows=[];ode=True;pressure=True;missing_pressure=True;mass=True;rebuilt=True;different=True
    transverse=quad(lambda x:exp(-2*x*x/(.25-x*x)) if abs(x)<.5 else 0,-.5,.5,epsabs=1e-12)[0]
    for L in (40,160):
        normalized=[]
        for sign in (-1,1):
            data=phase_module().phase_data(L*3028536,sign,0);c0=coefficients(data,0)
            initial=np.r_[c0['Ba'][:,0],c0['Br'][:,0]]
            def rhs(v,y):
                A=coefficients(data,v)['normalized'];return np.r_[A@y[:3],A@y[3:]]
            physical=solve_ivp(rhs,(0,L),initial,method='Radau',rtol=3e-11,atol=3e-13,dense_output=True)
            frame=solve_ivp(lambda v,z:coefficients(data,v)['reduced']@z,(0,L),(1.,0.),method='Radau',rtol=3e-11,atol=3e-13,dense_output=True)
            if not physical.success or not frame.success:raise RuntimeError('Independent ODE failed')
            error=0.;constraint=0.;balance=0.;omission=0.;difference=0.
            for v in np.linspace(0,L,81):
                c=coefficients(data,v);y=physical.sol(v);t=y[3:]
                error=max(error,np.linalg.norm(t-c['Br']@frame.sol(v)))
                constraint=max(constraint,abs(c['n']@t))
                pinorm=1j*(c['n']@c['K']@t-c['np']@t)/(data['k']*(c['n']@c['n']))
                residual=c['normalized']@t+c['rate']*t+c['K']@t+c['d']*t+1j*data['k']*c['n']*pinorm
                balance=max(balance,np.linalg.norm(residual));difference=max(difference,np.linalg.norm(y[:3]-t))
                omission=max(omission,np.linalg.norm(c['n']*(c['np']@t)/(c['n']@c['n'])))
            ode &= error<2e-8 and constraint<2e-9
            pressure &= balance<2e-9
            missing_pressure &= omission>1e-8
            u=float(data['u']);lam=float(data['lambda0']);J=1.25
            def integral(which,i):
                def f(w):
                    v=L/2+sqrt(L)*w;t=physical.sol(v)[which:which+3]
                    return cutoff(v/L)**2*exp(2*logP(v,L,u,lam))*t[0]*t[i+1]
                return J*transverse*quad(f,-sqrt(L)/3,sqrt(L)/3,points=[0],epsabs=2e-10,epsrel=2e-10,limit=150)[0]
            ha=np.array([integral(0,i) for i in range(2)]);hr=np.array([integral(3,i) for i in range(2)])
            normalized.append((ha,hr))
            mass &= bool(np.all(np.isfinite(ha)) and np.all(np.isfinite(hr)) and np.linalg.norm(hr)>0)
            different &= difference>1e-10
            rows.append(dict(length=L,sign=sign,independent_frame_physical_error=float(error),normal_constraint=float(constraint),
                complete_pressure_balance_error=float(balance),maximum_normalized_pulse_difference=float(difference),
                missing_nprime_pressure_defect=float(omission),
                actual_normalized_Haar_column=ha.tolist(),canonical_normalized_Haar_column=hr.tolist(),
                normalized_column_difference=float(np.linalg.norm(hr-ha))))
        Ha=np.column_stack([x[0] for x in normalized]);Hr=np.column_stack([x[1] for x in normalized])
        target=np.array([.8,.6]);ya=np.linalg.solve(Ha,target);yr=np.linalg.solve(Hr,target)
        defect=np.linalg.norm(Hr@ya-target);rebuilt &= bool(np.all(ya>0) and np.all(yr>0) and np.linalg.norm(Hr@yr-target)<1e-12)
        different &= defect>1e-12
        for row in rows[-2:]:row.update(stale_root_principal_covariance_defect=float(defect),
            rebuilt_principal_covariance_error=float(np.linalg.norm(Hr@yr-target)),
            relative_positive_root_difference=float(np.max(np.abs(np.sqrt(yr/ya)-1))))
    return dict(rows=rows,normalized_transverse_integral=float(transverse),checks=dict(
        independent_reference_frame_and_physical_ODE_agree=bool(ode),
        full_moving_normal_pressure_satisfies_principal_balance=bool(pressure),
        omitting_nprime_pressure_is_detected=bool(missing_pressure),
        true_cutoff_and_angular_half_Haar_columns_retained=bool(mass),
        independently_rebuilt_positive_roots_match_same_target=bool(rebuilt),
        resetting_frames_and_reusing_roots_detectably_changes_pulse_and_covariance=bool(different),
        small_fixtures_do_not_certify_giant_value_acceptance_or_continuum_suprema=True))
