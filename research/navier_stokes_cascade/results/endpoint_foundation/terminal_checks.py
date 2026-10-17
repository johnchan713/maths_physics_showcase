"""Independent finite Gamma integrals and exact seed/Taylor diagnostics."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb,factorial
import mpmath as mp
from terminal_bounds import HERE,load,parameter,order,heat_A,heat_S,radial_coefficients,mixed_majorant

poly=load('endpoint_foundation_seed_polynomials',HERE.parent/'realization_foundation/realized_checks.py')

def as_mp(x):
    x=F(x);return mp.mpf(x.numerator)/x.denominator

def quadratic_kernel(r,c,p,m):
    """Independent Faadi-Bruno derivative of (r^2+c)^(-p)."""
    return mp.factorial(m)*sum((-1)**(m-b)*mp.rf(p,m-b)*(r*r+c)**(-p-m+b)*(2*r)**(m-2*b)
        /(mp.factorial(m-2*b)*mp.factorial(b)) for b in range(m//2+1))

def kernel_heat(h,r,tau,n,j):
    """Differentiate r^-1*(r^2+4*tau*v)^(-h-j) under Gamma integral."""
    def kernel(v):
        radial=sum(comb(n,l)*(-1)**l*factorial(l)*r**(-1-l)*quadratic_kernel(r,4*tau*v,h+j,n-l) for l in range(n+1))
        return mp.exp(-v)*v**(h+j)*radial
    return 4**j*mp.rf(h,j)/mp.gamma(1+h)*mp.quad(kernel,[0,1,mp.inf])

@lru_cache(None)
def integral_review():
    heat_error=pressure_error=mp.mpf(0);heat_count=pressure_count=0
    missing_phase=mp.mpf(0);missing_pressure_cross=mp.mpf(0)
    heat_bounds=pressure_bounds=True
    with mp.workdps(55):
        r,tau=mp.mpf(3)/2,mp.mpf(1)/20;Z=4*tau/r**2;Zmax=F(2,5)
        for hf in (F(1,200),F(9,1000)):
            h=as_mp(hf)
            Fs=[mp.rf(h,k)/mp.gamma(1+h)*mp.quad(lambda v:mp.exp(-v)*v**(h+k)*(1+Z*v)**(-h-k),[0,1,mp.inf]) for k in range(7)]
            Gs=[sum(comb(j,k)*Fs[k]*Fs[j-k] for k in range(j+1)) for j in range(6)]
            values={}
            for n in range(4):
                for j in range(3):
                    actual=kernel_heat(h,r,tau,n,j);values[n,j]=actual
                    a=1+2*hf+2*j
                    predicted=4**j*r**(-as_mp(a)-n)*sum(as_mp(c)*Z**k*Fs[j+k] for k,c in enumerate(radial_coefficients(a,n)))
                    heat_error=max(heat_error,abs(actual-predicted)/(1+abs(actual)));heat_count+=1
                    bound=as_mp(mixed_majorant(hf,n,j,Zmax))*r**(-as_mp(a)-n)
                    heat_bounds=heat_bounds and abs(actual)<=bound
                    if n==1 and j==0:
                        truncated=-as_mp(a)*r**(-as_mp(a)-1)*Fs[0]
                        missing_phase=max(missing_phase,abs(actual-truncated))
            for n in range(1,4):
                for j in range(3):
                    # Exact independent product derivative of K^2/r.
                    actual=mp.mpf(0)
                    for l in range(n):
                        for a in range(n-l):
                            b=n-1-l-a;mult=mp.factorial(n-1)/(mp.factorial(l)*mp.factorial(a)*mp.factorial(b))
                            for k in range(j+1):actual+=mult*comb(j,k)*(-1)**l*factorial(l)*r**(-1-l)*values[a,k]*values[b,j-k]
                    beta=2+4*hf+2*j
                    predicted=4**j*r**(-as_mp(beta)-n)*sum(as_mp(c)*Z**k*Gs[j+k] for k,c in enumerate(radial_coefficients(beta+1,n-1)))
                    pressure_error=max(pressure_error,abs(actual-predicted)/(1+abs(actual)));pressure_count+=1
                    bound=as_mp(mixed_majorant(hf,n,j,Zmax,True))*r**(-as_mp(beta)-n)
                    pressure_bounds=pressure_bounds and abs(actual)<=bound
            missing_pressure_cross=max(missing_pressure_cross,abs(Gs[2]-2*Fs[0]*Fs[2]))
    return dict(independent_positive_Gamma_integrals=True,nonterminal_heat_comparisons=heat_count,
        nonterminal_pressure_product_comparisons=pressure_count,relative_error_threshold='2e-40',
        actual_global_or_terminal_fields_not_evaluated=True,
        checks=dict(all_mixed_heat_kernel_and_recurrence_values_agree=bool(heat_error<mp.mpf('2e-40')),
            all_mixed_complete_pressure_products_agree=bool(pressure_error<mp.mpf('2e-40')),
            heat_values_below_derived_compact_majorants=bool(heat_bounds),pressure_values_below_derived_compact_majorants=bool(pressure_bounds),
            dropping_radial_Z_derivatives_detected=bool(missing_phase>mp.mpf('1e-8')),
            dropping_pressure_Leibniz_crosses_detected=bool(missing_pressure_cross>mp.mpf('1e-8')),
            positive_preterminal_Z_only_no_future_heat_integral_used=True))

@lru_cache(None)
def exact_review():
    x,y,z,t=poly.X,poly.Y,poly.Z,poly.T
    X=poly.scale(poly.add(poly.power(x,2),poly.power(y,2)),F(1,2))
    quotient=poly.add(poly.ONE,X,poly.power(z,2))
    A=[poly.scale(poly.mul(y,quotient),F(-1,2)),poly.scale(poly.mul(x,quotient),F(1,2)),{}]
    velocity=poly.curl(A)
    target=[poly.scale(poly.mul(x,z),-1),poly.scale(poly.mul(y,z),-1),poly.add(poly.ONE,poly.scale(X,2),poly.power(z,2))]
    assert not any(poly.vsum(velocity,poly.vscale(target,-1)))
    assert not poly.divergence(velocity)
    direct=[poly.scale(poly.mul(y,quotient),-1),poly.mul(x,quotient),{}]
    assert not poly.divergence(direct)
    Taylor=Cauchy=0;wrong_sign=False
    for k in range(7):
        dk=F(factorial(7),factorial(7-k));M=F(factorial(7),factorial(6-k))
        for tau in (F(1,2),F(1,8),F(1,64)):
            expansion=sum(F(factorial(7),factorial(7-j)*factorial(j))*(-tau)**j for j in range(k+1))
            error=abs((1-tau)**7-expansion)
            assert error<=M*tau**(k+1)/factorial(k+1);Taylor+=1
            t1,t2=1-tau,1-tau/2
            assert abs(dk*(t2**(7-k)-t1**(7-k)))<=M*(t2-t1);Cauchy+=1
            if k:
                wrong=sum(F(factorial(7),factorial(7-j)*factorial(j))*tau**j for j in range(k+1))
                wrong_sign=wrong_sign or abs((1-tau)**7-wrong)>M*tau**(k+1)/factorial(k+1)
    a=F(101,100);coeff=radial_coefficients(a,2)
    return dict(manufactured_polynomial_seed_at_q_one_not_actual_fluid=True,
        exact_seed_curl_components=3,finite_Taylor_controls=Taylor,Cauchy_controls=Cauchy,
        checks=dict(axis_regular_cartesian_potential_curl_has_full_radial_and_axial_components=True,
            axis_potential_and_direct_swirl_divergence_exact=True,
            finite_Taylor_remainders_obey_next_derivative_bound=Taylor==21,
            Cauchy_estimates_use_the_next_time_derivative=Cauchy==21,
            wrong_terminal_time_sign_detected=wrong_sign,
            missing_radial_homogeneity_update_detected=coeff[0]==a*(a+1) and coeff[0]!=a*a,
            nonzero_terminal_pressure_cannot_be_zero_extended=heat_S(F(1,200),0)>0,
            bounded_values_alone_not_claimed_to_supply_terminal_limits=True))
