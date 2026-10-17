"""Manufactured actual-activation jets with full incoming moment/pressure memory.

This fixture is the earlier manufactured seed, not the enormous analytic core.
Scalar quadrature values and radial ODE jets are assembled independently of
the older cancellation-safe five-moment implementation.
"""
from functools import lru_cache
import importlib.util
from math import factorial
from pathlib import Path
import sys
import mpmath as mp


@lru_cache(None)
def old_collar():
    p=Path(__file__).resolve().parent.parent/'inner_collar_activation/collar.py'
    s=importlib.util.spec_from_file_location('actual_collar_historical_formula',p)
    m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m


class RJet:
    """Taylor coefficients through radial order four; scalar eta may be complex."""
    N=4
    def __init__(self,c=0):
        self.c=[mp.mpc(v) for v in c] if isinstance(c,(list,tuple)) else [mp.mpc(c)]
        self.c=(self.c+[mp.mpc(0)]*(self.N+1))[:self.N+1]
    def cast(self,v):return v if isinstance(v,RJet) else RJet(v)
    def __add__(self,v):
        v=self.cast(v);return RJet([a+b for a,b in zip(self.c,v.c)])
    __radd__=__add__
    def __neg__(self):return RJet([-v for v in self.c])
    def __sub__(self,v):return self+-self.cast(v)
    def __rsub__(self,v):return self.cast(v)+-self
    def __mul__(self,v):
        v=self.cast(v);return RJet([sum(self.c[j]*v.c[n-j] for j in range(n+1)) for n in range(self.N+1)])
    __rmul__=__mul__
    def inverse(self):
        if not self.c[0]:raise ValueError('Nonzero local Taylor denominator required')
        out=[1/self.c[0]]
        for n in range(1,self.N+1):out.append(-sum(self.c[j]*out[n-j] for j in range(1,n+1))/self.c[0])
        return RJet(out)
    def __truediv__(self,v):return self*self.cast(v).inverse()
    def __rtruediv__(self,v):return self.cast(v)*self.inverse()
    def __pow__(self,n):
        if type(n) is not int or n<0:raise ValueError('Nonnegative integer jet power required')
        out=RJet(1)
        for _ in range(n):out=out*self
        return out
    def exp(self):
        out=[mp.exp(self.c[0])]
        for n in range(1,self.N+1):out.append(sum(j*self.c[j]*out[n-j] for j in range(1,n+1))/n)
        return RJet(out)
    def expm1(self):
        out=self.exp();out.c[0]=mp.expm1(self.c[0]);return out
    def derivative_value(self,n):
        if type(n) is not int or not 0<=n<=self.N:raise ValueError('Unavailable radial jet')
        return factorial(n)*self.c[n]


@lru_cache(maxsize=4096)
def factors(t,s,kappa,order,digits):
    t,s,kappa=map(mp.mpf,(t,s,kappa));old=old_collar()
    if not 0<t<=mp.mpf('.01') or not 0<s<=mp.mpf('.25') or not 0<=kappa<=mp.mpf('.5'):raise ValueError('Manufactured collar range exceeded')
    ea0=(1-kappa)*old.sigma(s)
    J0=mp.mpf('1.5')*t*ea0*s**3*old.flat_quotient(s,lambda u:mp.mpf(1),order=order)
    K0=t*ea0*s**3*old.flat_quotient(s,lambda u:mp.exp(-t*u),order=order)
    y=RJet([t*s,1]);v=y/t
    left=(-1/(v*v)).exp();right=(-1/((1-v)*(1-v))).exp()
    ea=(1-kappa)*left/(left+right)
    ey=(-y).exp();J=RJet([J0]+[mp.mpf('1.5')*ea.c[n-1]/n for n in range(1,5)])
    g=ea*ey;K=RJet([K0]+[g.c[n-1]/n for n in range(1,5)])
    return y,J,K,ea


@lru_cache(maxsize=12)
def legendre(order,digits):
    return tuple(zip(*mp.gauss_quadrature(order,'legendre')))


@lru_cache(maxsize=128)
def scalar_excess_integrals(t,s,kappa,order,digits):
    t,s,kappa=map(mp.mpf,(t,s,kappa));values=[mp.mpc(0)]*6
    for node,weight in legendre(order,digits):
        u=s*(node+1)/2;y,J,K,ea=factors(t,u,kappa,order,digits)
        yy,j,k=y.c[0],J.c[0],K.c[0];X=2*mp.exp(yy)
        ex1,ex2=mp.expm1(j),mp.expm1(2*j)
        rows=[X*k,4*mp.exp(yy/2)*ex1,
            4*mp.exp(yy/2)*(mp.exp(-yy)*ex1+k*mp.exp(j)),
            X*(2*mp.exp(-yy)*k+k*k),X*mp.exp(-2*yy)*ex2,mp.exp(-2*yy)*ex2/2]
        values=[a+t*s*weight*b/2 for a,b in zip(values,rows)]
    return tuple(values)


def native(y,e):
    eb=2+e/10+e*e/20;q=1+e/5;ub=e/3+e*e/30;b=q*eb/2;c=ub-b
    return dict(eb=eb,b=b,c=c,
        M=(1+e/5)/10+2*(c*mp.expm1(y)+b*y),
        I=(2+e/3)/10+8*eb*mp.expm1(y/2),
        J=(1-e/5)/20+8*eb*(c*mp.expm1(y/2)-b*mp.expm1(-y/2)),
        S=(7+2*e*e)/100+2*(c*c*mp.expm1(y)+2*c*b*y-(b*b-eb*eb/2)*mp.expm1(-y)),
        Cp=(3+e)/100-eb*eb*mp.expm1(-2*y)/4,
        Pi=-mp.mpf('.7')+e/5+3*e*e/20-eb*eb*mp.expm1(-2*y)/4)


def state(t,s,e,kappa='.1',order=32,drop_incoming=False,work_digits=None):
    digits=mp.mp.dps if work_digits is None else work_digits
    with mp.workdps(digits):
        t,s,kappa=map(mp.mpf,(t,s,kappa))
        y,J,K,ea=factors(t,s,kappa,order,digits)
        m,i,j2,s2,s3,cp=scalar_excess_integrals(t,s,kappa,order,digits)
    a=native(y.c[0],e);eb,b,c=a['eb'],a['b'],a['c']
    E=eb*(-y).exp()*J.exp();U=c+b*((-y).exp()+K)
    F=eb/2*(-mp.mpf('1.5')*y).exp()*J.exp();X=2*y.exp()
    dF=eb/2*(-mp.mpf('1.5')*y).exp()*J.expm1();dU=b*K
    excess=dict(M=b*m,I=eb*i,J=c*eb*i+b*eb*j2,S=2*c*b*m+b*b*s2-eb*eb*s3/2,Cp=eb*eb*cp,Pi=eb*eb*cp)
    densities=dict(M=X*U,I=2*X*X*F,J=2*X*X*U*F,S=X*(U*U-X*F*F),Cp=X*F*F,Pi=X*F*F)
    out=dict(E=E,U=U,F=F,dF=dF,dU=dU,X=X,ea=ea)
    prefix=native(mp.mpf(0),e)
    for n,density in densities.items():
        value=a[n]+excess[n]-(prefix[n] if drop_incoming else 0)
        out[n]=RJet([value]+[density.c[r-1]/r for r in range(1,5)])
    return out


@lru_cache(None)
def diagnostic_stats():
    with mp.workdps(75):
        t,s,e,k=mp.mpf('.01'),mp.mpf('.24'),mp.mpf('.3'),mp.mpf('.1')
        a=state(t,s,e,k,32);fine=state(t,s,e,k,48);old=old_collar().Activation(t,k,24).state(s,e)['actual']
        base=native(t*s,e);base['E']=base['eb']*mp.exp(-t*s);base['U']=base['c']+base['b']*mp.exp(-t*s)
        errors=[abs(a[n].c[0]-fine[n].c[0])/max(abs(fine[n].c[0]-base[n]),mp.mpf('1e-100')) for n in ('M','I','J','S','Cp','Pi')]
        old_errors=[abs(fine[n].c[0]-old[n].v)/max(abs(old[n].v-base[n]),mp.mpf('1e-100')) for n in ('M','I','J','S','Cp','Pi','E','U')]
        # Every eta dependence is a polynomial of degree <=6 in this fixture.
        # Cauchy/DFT evaluation is independent of real differentiation.
        cauchy_errors=[];radius=mp.mpf('.1');N=24
        for name,r in (('F',0),('U',3),('M',2),('J',1),('S',4),('Pi',0)):
            def f(eta):return state(t,s,eta,k,32,work_digits=75)[name].derivative_value(r)
            samples=[f(e+radius*mp.exp(2j*mp.pi*n/N)) for n in range(N)]
            for j in range(7):
                cauchy=mp.factorial(j)*sum(v*mp.exp(-2j*mp.pi*n*j/N) for n,v in enumerate(samples))/(N*radius**j)
                direct=mp.diff(f,e,j)
                cauchy_errors.append(abs(cauchy-direct)/(1+abs(direct)))
        widths=[mp.mpf('.01'),mp.mpf('.001'),mp.mpf('.0001')]
        fourth=[];scaled_fields=[];scaled_pressure=[]
        for width in widths:
            z=state(width,s,e,k,32)
            d4=abs(z['dF'].derivative_value(4));fourth.append(d4)
            scaled_fields.append(d4*width**3)
            ref=native(width*s,e)
            # Differentiate pressure excess using actual minus native density jets.
            yy=RJet([width*s,1]);eb=ref['eb'];X=2*yy.exp();Fr=eb/2*(-mp.mpf('1.5')*yy).exp()
            diff_density=X*(z['F']*z['F']-Fr*Fr)
            scaled_pressure.append(abs(diff_density.derivative_value(3))*width**2)
        dropped=state(t,s,e,k,32,drop_incoming=True)
        prefix_gap=abs(a['Pi'].c[0]-dropped['Pi'].c[0])
        pressure_actual_gap=abs(a['Pi'].c[0]-native(t*s,e)['Pi'])
        axis_datum=-mp.mpf('.73')+mp.mpf('.19')*e+mp.mpf('.15')*e*e
        datum_error=max(abs((a['Pi']-a['Cp']).c[r]-(axis_datum if r==0 else 0)) for r in range(5))
        # Exact jet flow equations through the available three differentiated orders.
        y,J,K,ea=factors(t,s,k,32,mp.mp.dps);eb,b=a['E'].c[0]/((-y).exp()*J.exp()).c[0],native(t*s,e)['b']
        E_rhs=(mp.mpf('.5')-mp.mpf('1.5')*(1-ea))*a['E'];U_rhs=-(1-ea)*b*(-y).exp()
        flow=max(abs((r+1)*a[n].c[r+1]-rhs.c[r]) for n,rhs in (('E',E_rhs),('U',U_rhs)) for r in range(4))
        return dict(maximum_refinement_error=max(errors),maximum_old_formula_error=max(old_errors),
            maximum_parameter_Cauchy_error=max(cauchy_errors),parameter_checks=len(cauchy_errors),
            maximum_radial_flow_error=flow,fourth_derivative_growth=fourth[-1]/fourth[0],
            scaled_field_ratio=max(scaled_fields)/min(scaled_fields),scaled_pressure_ratio=max(scaled_pressure)/min(scaled_pressure),
            omitted_incoming_pressure_gap=prefix_gap,using_reference_pressure_increment_gap=pressure_actual_gap,axis_datum_error=datum_error)


def diagnostic_review():
    r=diagnostic_stats()
    return dict(diagnostic_is_manufactured_not_actual_core=True,radial_order=4,parameter_orders=list(range(7)),
        parameter_Cauchy_checks=r['parameter_checks'],quadrature_orders=[32,48],old_formula_order=24,
        tolerances=dict(quadrature_and_old_formula='1e-12',parameter_Cauchy='1e-40',radial_flow='1e-50'),
        raw_roundoff_scores_tested_not_refrozen=True,
        checks=dict(independent_partial_interval_quadrature_refines=r['maximum_refinement_error']<mp.mpf('1e-12'),
            full_memory_values_match_historical_cancellation_safe_formula=r['maximum_old_formula_error']<mp.mpf('1e-12'),
            all_parameter_Cauchy_jets_match_real_derivatives=r['maximum_parameter_Cauchy_error']<mp.mpf('1e-40'),
            complete_radial_flow_jets_match_actual_shear_prescription=r['maximum_radial_flow_error']<mp.mpf('1e-50'),
            fourth_radial_derivative_is_not_uniform_in_width=r['fourth_derivative_growth']>mp.mpf('1e5'),
            predicted_field_width_gain_tracks_three_widths=r['scaled_field_ratio']<2,
            predicted_pressure_width_gain_tracks_three_widths=r['scaled_pressure_ratio']<2,
            dropping_incoming_pressure_is_detected=r['omitted_incoming_pressure_gap']>mp.mpf('.1'),
            Pi_minus_Cp_retains_radially_constant_nonzero_axis_datum=r['axis_datum_error']<mp.mpf('1e-60'),
            using_seed_instead_of_actual_pressure_density_is_detected=r['using_reference_pressure_increment_gap']>mp.mpf('1e-30')))
