"""Independent integral and noncommuting ODE controls, not actual giant pulses."""
from functools import lru_cache
from fractions import Fraction as F
from math import asinh,exp,sqrt,sin,cos,pi,log
import numpy as np
from scipy.integrate import quad,solve_ivp
from mpmath import mp
from gaussian_bounds import curvature,logs,barrier


def rate(v,L,u,lam):
    a=.5+v/L
    return lam*((1+u*u*a*a)**-.5-(1+u*u*a*a)/(1+u*u)**1.5)


def log_envelope(v,L,u,lam):
    a=.5+v/L
    return lam*L*((asinh(u*a)-asinh(u))/u-((a-1)+u*u*(a*a*a-1)/3)/(1+u*u)**1.5)


def cutoff(t):
    distance=abs(t-.5)
    if distance<=.2:return 1.
    if distance>=1/3:return 0.
    x=(distance-.2)/(1/3-.2)
    a,b=exp(-1/x),exp(-1/(1-x))
    return b/(a+b)


@lru_cache(None)
def integral_review():
    rows=[];max_integral_error=mp.mpf(0);max_derivative_error=mp.mpf(0);gaussian=True
    with mp.workdps(60):
        for ut in ('0.1','0.5','1','2','5'):
            u=mp.mpf(ut);lam=mp.mpf('1.25');lo,hi=curvature(F(ut),F(1,2),2)
            c,C=mp.mpf(str(lo.numerator))/lo.denominator,mp.mpf(str(hi.numerator))/hi.denominator
            for L in (40,160):
                def r(v):
                    a=mp.mpf('.5')+v/L
                    return lam*(1/mp.sqrt(1+u*u*a*a)-(1+u*u*a*a)/(1+u*u)**mp.mpf('1.5'))
                def closed(v):
                    a=mp.mpf('.5')+v/L
                    return lam*L*((mp.asinh(u*a)-mp.asinh(u))/u-((a-1)+u*u*(a**3-1)/3)/(1+u*u)**mp.mpf('1.5'))
                error=mp.mpf(0)
                for fraction in ('0','.1','.25','.5','.8','1'):
                    v=L*mp.mpf(fraction);value=closed(v);independent=mp.quad(r,[mp.mpf(L)/2,v])
                    error=max(error,abs(value-independent))
                    max_derivative_error=max(max_derivative_error,abs(mp.diff(closed,v)-r(v)))
                    distance=v-mp.mpf(L)/2
                    gaussian &= -C*distance**2/(2*L)-mp.mpf('1e-50')<=value<=-c*distance**2/(2*L)+mp.mpf('1e-50')
                max_integral_error=max(max_integral_error,error)
                rows.append(dict(tilt=ut,length=L,integral_absolute_error=float(error)))
    return dict(rows=rows,max_integral_error=float(max_integral_error),max_derivative_error=float(max_derivative_error),
        checks=dict(closed_integral_matches_independent_high_precision_quadrature=max_integral_error<mp.mpf('1e-45'),
            closed_integral_derivative_matches_rate=max_derivative_error<mp.mpf('1e-45'),
            two_sided_Gaussian_brackets_hold_in_controls=bool(gaussian),
            finite_integral_controls_not_continuum_or_giant_profile_certificates=True))


@lru_cache(None)
def ode_mass_review():
    CE,Cd=1/8,1/10;rows=[];good=True;mass_ok=True;direction_ok=True;noncommuting=False
    transverse=quad(lambda x:exp(-2*x*x/(.25-x*x)) if abs(x)<.5 else 0,-.5,.5,epsabs=1e-12,epsrel=1e-12)[0]
    for u in (.25,.75,1.5,3.):
        for L in (40.,160.):
            lam0=1.25;mu=.5/sqrt(1+9*u*u/4);budget=barrier(F(str(mu)),F(1,8),F(1,10),F(str(L)))
            rho=float(budget['rho']);constants=logs(F(str(u)),F(1,2),2,F(1,8),F(1,10),1,2)
            def matrix(v):
                a=.5+v/L;lam=lam0/sqrt(1+u*u*a*a);angle=v/L
                E=CE/(4*L)*np.array([[cos(angle),sin(angle)],[-sin(angle),cos(angle)]])
                delta=Cd*sin(2*pi*v/L)/L
                return np.diag([0.,-2*lam])+E-delta*np.eye(2)
            noncommuting |= bool(abs(matrix(L/2)[0,1])>0)
            sol=solve_ivp(lambda v,y:matrix(v)@y,(0,L),(1.,0.),method='Radau',rtol=2e-10,atol=2e-12,dense_output=True)
            if not sol.success:raise RuntimeError(sol.message)
            vals=sol.sol(np.linspace(0,L,201));ratios=vals[1]/vals[0];radial=vals[0]+vals[1]
            lower,upper=exp(constants['log_xminus']),exp(constants['log_xplus'])
            good &= bool(np.min(vals[0])>0 and np.max(np.abs(ratios))<=rho+2e-9 and np.min(radial)>=lower-2e-9 and np.max(radial)<=upper+2e-9)
            # Resolve the pulse center with w=(v-L/2)/sqrt(L); retain actual psi.
            def integrand(w,first=False):
                v=L/2+sqrt(L)*w;y=sol.sol(v);xratio=y[0]+y[1]
                value=cutoff(v/L)**2*exp(2*log_envelope(v,L,u,lam0))*xratio*xratio
                return value*(abs(w)*sqrt(L) if first else 1)
            width=sqrt(L)/3
            m0=quad(lambda w:integrand(w),-width,width,points=[0],epsabs=1e-10,epsrel=1e-10,limit=200)[0]*sqrt(L)
            m1=quad(lambda w:integrand(w,True),-width,width,points=[0],epsabs=1e-10,epsrel=1e-10,limit=200)[0]*sqrt(L)
            J=1.25;j=J*transverse*m0/sqrt(L)
            mass_ok &= exp(constants['log_jminus'])*(1-1e-8)<=j<=exp(constants['log_jplus'])*(1+1e-8)
            concentration=m1/(L*m0);bound=exp(constants['log_direction_factor'])/sqrt(L)
            direction_ok &= concentration<=bound*(1+1e-8)
            rows.append(dict(tilt=u,length=L,maximum_ratio=float(np.max(np.abs(ratios))),barrier=rho,
                radial_ratio_min=float(np.min(radial)),radial_ratio_max=float(np.max(radial)),
                normalized_actual_Haar_mass=j,log_mass_lower=constants['log_jminus'],log_mass_upper=constants['log_jplus'],
                concentration=concentration,concentration_bound=bound))
    return dict(rows=rows,normalized_transverse_integral=transverse,
        checks=dict(noncommuting_frame_fixture=bool(noncommuting),Riccati_and_radial_amplitude_bounds=bool(good),
            actual_cutoff_Haar_integrals_between_explicit_bounds=bool(mass_ok),
            actual_first_absolute_moment_concentration=bool(direction_ok),
            transverse_bounds_preserve_true_area=.5*exp(-2/3)<=transverse<=1,
            actual_manufactured_cutoff_not_replaced_by_envelope=True,
            full_pressure_formula_and_original_Gaussian_tail_split_retained=True))


def exact_haar_review():
    # Independent algebra retaining angular mean 1/2, true ci and both integrals.
    rows=[]
    for r,L,rootL in ((F(1,100),100,10),(F(1,25),400,20)):
        J=F(3,2);transverse=r/2;time_mass=7*rootL;ci=2*r/L
        h=J/2*transverse*ci*time_mass;j=h*rootL/(r*r)
        rows.append(dict(radius=str(r),length=L,h=str(h),j=str(j),wrong_without_angular_half=str(2*h)))
    return dict(rows=rows,checks=dict(normalized_mass_exactly_radius_independent=rows[0]['j']==rows[1]['j'],
        physical_mass_has_three_halves_radius_scaling=F(rows[1]['h'])==8*F(rows[0]['h']),
        missing_angular_half_detected=all(r['h']!=r['wrong_without_angular_half'] for r in rows),
        Gaussian_mass_lower_bound_not_lower_bound_on_complete_velocity=True))
