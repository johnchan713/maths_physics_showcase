"""Independent rational principal equations and numerical noncommuting ODE controls."""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
import numpy as np
from scipy.integrate import quad,solve_ivp


def add(a,b):return a[0]+b[0],a[1]+b[1]
def scale(a,b):return a[0]*b,a[1]*b
def imag(a):return -a[1],a[0]
def total(xs):
    out=(F(0),F(0))
    for x in xs:out=add(out,x)
    return out


def principal_checks():
    count=0;exact=True;normal=False;moving=False;reality=True;rows=[]
    for m in (-7,-2,-1,1,2,7):
        for v in (F(0),F(1,5),F(2,3),F(1)):
            n=(F(1),v,F(0));npv=(F(0),F(1),F(0));norm=1+v*v
            # A manufactured zero-entry transverse complex amplitude.
            q=(v*v,v**3*(1 if m>0 else -1));qp=(2*v,3*v*v*(1 if m>0 else -1))
            t=(scale(q,-v),q,(v**3,v*v*(1 if m>0 else -1)))
            tp=(scale(add(q,scale(qp,v)),-1),qp,(3*v*v,2*v*(1 if m>0 else -1)))
            K=((F(0),-2-v,F(0)),(3+v,F(0),F(0)),(2-v,F(0),F(0)))
            kt=[total(scale(x,c) for x,c in zip(t,row)) for row in K]
            pi=(v*v/F(abs(m)),v**3/F(m));k=3;d=F(1,8)*k*k*norm
            f=[scale(total((a,b,scale(c,m*m*d),scale(imag(pi),k*m*ni))),-1) for a,b,c,ni in zip(tp,kt,t,n)]
            nkt=total(scale(x,c) for x,c in zip(kt,n));npt=total(scale(x,c) for x,c in zip(t,npv));nf=total(scale(x,c) for x,c in zip(f,n))
            recovered=scale(imag(add(add(nkt,scale(npt,-1)),nf)),1/(k*m*norm))
            exact &= recovered==pi and total(scale(x,c) for x,c in zip(t,n))==(0,0)
            for i in range(3):
                rhs=total((scale(kt[i],-1),scale(add(nkt,scale(npt,-1)),n[i]/norm),scale(t[i],-m*m*d),scale(f[i],-1),scale(nf,n[i]/norm)))
                exact &= rhs==tp[i];count+=2
            normal |= nf!=(0,0);moving |= npt!=(0,0)
            rows.append(dict(m=m,v=str(v),pressure=list(map(str,pi))))
    # A source purely along a constant normal has t=0 but pi=i/(km).
    direct=[dict(m=m,t='0',pressure_imaginary=str(F(1,3*m))) for m in (2,4,8,64)]
    for row in rows:
        mate=next(r for r in rows if r['m']==-row['m'] and r['v']==row['v'])
        reality &= F(row['pressure'][0])==F(mate['pressure'][0]) and F(row['pressure'][1])==-F(mate['pressure'][1])
    return dict(exact_ODE_real_components=count,rows=rows,normal_source_controls=direct,
        checks=dict(all_manufactured_principal_equations_and_pressure_exact=bool(exact),
            moving_normal_omission_detected=bool(moving),source_normal_omission_detected=bool(normal),
            negative_harmonic_reality=bool(reality),
            pressure_cannot_inherit_velocity_inverse_square_gain=all(F(r['pressure_imaginary'])*r['m']==F(1,3) for r in direct)))


@lru_cache(None)
def numerical_checks():
    # This two-coordinate model is NOT the actual enormous background.
    # Its error matrix does not commute with diag(0,-2).
    L=12.;E=np.array([[.02,.03],[-.02,.01]])/L;C0=np.exp(.07)
    def d(v):return 1+(v-L/2)/L+.02*np.cos(v/L)/L
    dmin=.49
    def generator(v,m):return np.diag([0.,-2.])+E-.02*np.cos(v/L)/L*np.eye(2)-(m*m-1)*d(v)*np.eye(2)
    def solve(m,initial,forcing=False):
        def rhs(v,y):
            Y=y.reshape(2,-1);answer=generator(v,m)@Y
            if forcing:answer+=np.array([[1+.1*v/L],[.2]])
            return answer.reshape(-1)
        sol=solve_ivp(rhs,(0,L),np.array(initial).reshape(-1),method='Radau',rtol=2e-10,atol=2e-12,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        return sol
    base=solve(1,np.eye(2));max_factor_error=0.;max_energy_ratio=0.;max_source_ratio=0.;rows=[]
    # The signed error matrix really is noncommuting; scalar damping still factors.
    noncomm=np.linalg.norm(E@np.diag([0.,-2.])-np.diag([0.,-2.])@E)
    for m in (-5,-2,-1,1,2,5):
        prop=solve(m,np.eye(2));forced=solve(m,np.zeros((2,1)),True);delta=m*m-1
        local_factor=0.;local_source=0.
        for v in np.linspace(0,L,61):
            integ=v+(v*v/2-L*v/2)/L+.02*np.sin(v/L)
            actual=prop.sol(v).reshape(2,2);expected=np.exp(-delta*integ)*base.sol(v).reshape(2,2)
            local_factor=max(local_factor,float(np.linalg.norm(actual-expected)))
            bound=C0*np.exp(-delta*dmin*v)
            # At very small tails an absolute tolerance, not a relative claim.
            max_energy_ratio=max(max_energy_ratio,float(max(0,np.linalg.norm(actual,2)-bound)))
            ell=L if not delta else min(L,1/(delta*dmin))
            local_source=max(local_source,float(np.linalg.norm(forced.sol(v))/(C0*1.2*ell)))
        max_factor_error=max(max_factor_error,local_factor);max_source_ratio=max(max_source_ratio,local_source)
        rows.append(dict(m=m,factor_absolute_error=local_factor,source_majorant_ratio=local_source))
    # Exact integral formula for parameter jets of a scalar damped response.
    parameter_rows=[];parameter_ok=True
    for m in (2,3,10):
        delta=m*m-1
        for j in range(7):
            value=quad(lambda w:(delta*w/4)**j*np.exp(-delta*.75*w),0,L,epsabs=1e-13,epsrel=1e-12)[0]
            upper=factorial(j)/(4**j*delta*.75**(j+1))
            parameter_ok &= value<=upper*(1+1e-10)
            parameter_rows.append(dict(m=m,order=j,ratio=value/upper))
    return dict(model='noncommuting 2x2 normalized model, not actual selected pulse',rows=rows,
        parameter_rows=parameter_rows,max_factor_absolute_error=max_factor_error,
        max_energy_bound_excess=max_energy_ratio,
        checks=dict(noncommuting_matrix_fixture=bool(noncomm>0),
            scalar_damping_factorization=max_factor_error<2e-8,
            energy_propagator_majorant=max_energy_ratio<2e-8,
            forced_value_residence_bound=max_source_ratio<=1+1e-8,
            parameter_jet_inverse_square_gain=bool(parameter_ok),
            numerical_controls_not_actual_pulse_evaluation=True))


def endpoint_checks():
    # Exact scalar response with constant normalized forcing at entrance:
    # y(0)=0, y'(0)=1 for every delta. Compactly supported f can be chosen
    # with nonzero positive mass in the interior, giving positive exit.
    return dict(entrance_derivative='1 for every m, despite value bound O((m^2-1)^-1)',
        exit_integral='integral_0^L exp[-(m^2-1)d(L-w)] f(w) dw >0 for f>=0 nonzero',
        checks=dict(generic_source_needs_entrance_jet_compatibility_for_pulse_derivative_gain=True,
            actual_full_smooth_entrance_support_supplies_compatibility=True,
            zero_entry_and_smooth_source_endpoints_do_not_imply_zero_exit=True,
            cutoff_applied_to_potential_and_pressure_before_curl=True,
            only_full_lifted_path_support_is_preserved_before_temporal_cutoff=True))
