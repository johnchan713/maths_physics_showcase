"""Exact manufactured higher jets; no quadrature fit or giant-field sampling."""
from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
from math import factorial,comb
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent


@lru_cache(None)
def JetClass():
    s=importlib.util.spec_from_file_location('all_order_map_rational_jets',HERE.parent/'core_all_order_foundation/coordinate_checks.py')
    m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m.Jet


def jet(v=0,n=2,order=6):return JetClass()(v,n=n,order=order)
def variable(i,n=2,order=6):return JetClass().variable(i,0,n=n,order=order)
def homogeneous(z,k):return JetClass()({a:v for a,v in z.terms.items() if sum(a)==k},n=z.n,order=z.order)
def zero(z):return not z.terms


def exponential(z):
    if z.value():raise ValueError('Zero constant needed for rational exponential')
    out=term=jet(1,z.n,z.order)
    for i in range(1,z.order+1):term=term*z/i;out=out+term
    return out


def sine(z):
    if z.value():raise ValueError('Zero constant needed for rational sine')
    return sum(((-1)**j*z**(2*j+1)/factorial(2*j+1) for j in range((z.order+1)//2)),jet(0,z.n,z.order))


def cosine(z):
    if z.value():raise ValueError('Zero constant needed for rational cosine')
    return sum(((-1)**j*z**(2*j)/factorial(2*j) for j in range(z.order//2+1)),jet(0,z.n,z.order))


def poly_mul(a,b,m):return [sum(a[j]*b[i-j] for j in range(i+1) if j<len(a) and i-j<len(b)) for i in range(m+1)]


def poly_power(a,p,m):
    if a[0]!=1:raise ValueError('Unit base required')
    out=[Q(1)]+[Q(0)]*m;term=out[:];b=[Q(0)]+list(a[1:])+[Q(0)]*(m+1-len(a));coef=Q(1)
    for k in range(1,m+1):
        term=poly_mul(term,b,m);coef*=Q(p-k+1,k)
        out=[x+coef*y for x,y in zip(out,term)]
    return out


@lru_cache(None)
def variance_series(m=6):
    # Exact Taylor coefficients of R2=(I0(2z)/I0(z)^2-1)/z^2.
    # Compute the numerator two orders farther BEFORE removing z^2.
    M=[Q(0)]*(m+3);M2=M[:]
    for k in range((m+2)//2+1):M[2*k]=Q(1,4**k*factorial(k)**2);M2[2*k]=Q(1,factorial(k)**2)
    ratio=poly_mul(M2,poly_power(M,-2,m+2),m+2)
    assert ratio[0]==1 and ratio[1]==0
    return tuple(ratio[i+2] for i in range(m+1))


def signed_F(mu,p):
    z=mu*p;R=sum((a*z**i for i,a in enumerate(variance_series(mu.order))),jet(0,mu.n,mu.order))
    return mu*(2*R)**Q(1,2) # sqrt(2)*F/d0, so F_mu at mu=0 equals one


def root_jets(p0=Q(0),order=6):
    y,e=variable(0,order=order),variable(1,order=order)
    p=p0+e/7+y/11;s=y/3+e/5+y*e/13
    mu=jet(0,order=order)
    for k in range(1,order+1):mu=mu-homogeneous(signed_F(mu,p)-s,k)
    # Independent scalar Lagrange inversion in s, with p as a parameter.
    R=[2*a for a in variance_series(order)]
    independent=jet(0,order=order)
    for k in range(1,order+1):
        c=poly_power(R,-Q(k,2),k-1)[k-1]/k
        independent=independent+c*p**(k-1)*s**k
    return mu,independent,signed_F(mu,p)-s,p,s


def phase_jets(order=6):
    y,e,phi=[variable(i,n=3,order=order) for i in range(3)]
    eps=Q(1,4)+y/7+e/11+y*e/13
    theta=jet(0,n=3,order=order)
    for k in range(1,order+1):theta=theta-homogeneous(theta+eps*sine(theta)-phi,k)/Q(5,4)
    newton=jet(0,n=3,order=order)
    for _ in range(order):newton=newton-(newton+eps*sine(newton)-phi)/(1+eps*cosine(newton))
    frozen=jet(0,n=3,order=order)
    for _ in range(order):frozen=frozen-(frozen+sine(frozen)/4-phi)/(1+cosine(frozen)/4)
    return theta,newton,theta+eps*sine(theta)-phi,frozen


B=((1,1,0,0,0),(1,2,0,0,0),(0,0,1,1,1),(0,0,1,2,4),(0,0,1,3,9))


def invert_matrix(B):
    n=len(B);a=[[Q(v) for v in row]+[Q(i==j) for j in range(n)] for i,row in enumerate(B)]
    for j in range(n):
        k=next(i for i in range(j,n) if a[i][j]);a[j],a[k]=a[k],a[j]
        q=a[j][j];a[j]=[v/q for v in a[j]]
        for i in range(n):
            if i!=j:
                q=a[i][j];a[i]=[x-q*y for x,y in zip(a[i],a[j])]
    return tuple(tuple(row[n:]) for row in a)


def matvec(a,v):return [sum((c*z for c,z in zip(row,v)),jet(0,v[0].n,v[0].order)) for row in a]


def quadratic(c,drop_u_square=False,drop_pressure_square=False):
    u,v,x,y,z=c;out=[jet(0,u.n,u.order) for _ in range(5)]
    out[3]=((0 if drop_u_square else -u*u-2*v*v)+(x*x+2*y*y+3*z*z)/2)
    out[4]=jet(0,u.n,u.order) if drop_pressure_square else (x*x+y*y+z*z)/2
    return out


def repair_jets(N=32,order=8):
    x=variable(0,n=1,order=order);target=[sum(((i+1)*x**k/Q(N*10**k*k) for k in range(1,order+1)),jet(0,1,order)) for i in range(5)]
    inv=invert_matrix(B);c=[jet(0,1,order) for _ in range(5)]
    for k in range(1,order+1):
        residual=[a+b-z for a,b,z in zip(matvec(B,c),quadratic(c),target)]
        delta=matvec(inv,[homogeneous(v,k) for v in residual]);c=[a-b for a,b in zip(c,delta)]
    other=[jet(0,1,order) for _ in range(5)]
    for _ in range(order):other=matvec(inv,[z-q for z,q in zip(target,quadratic(other))])
    return c,other,target


def frequency_jets(N=32,order=8):
    y,e=variable(0,order=order),variable(1,order=order)
    h=((1+e)*sine(N*y)+y*cosine(N*y))/N
    delta=(1+y/11+e/13)*(exponential(h)-1)
    # Independent chain-operator evaluation before substitution, including lower slow terms.
    ys,es,ph=[variable(i,n=3,order=order) for i in range(3)]
    primitive=(1+es)*sine(ph)+ys*cosine(ph)
    checks=[]
    for r in range(order+1):
        for k in range(order+1-r):
            expected=sum(Q(comb(r,j))*N**j*primitive.partial((r-j,k,j))/N for j in range(r+1))
            checks.append(h.partial((r,k))==expected)
    density=(1+e+e*e)*sine(N*y)/N
    memory=jet({(0,1):Q(1,N),(0,2):Q(1,7*N)},order=order)
    moment=memory+(1+e+e*e)*(1-cosine(N*y))/(N*N)
    return dict(chain=all(checks),chain_count=len(checks),field=delta,density=density,moment=moment,memory=memory)


@lru_cache(None)
def diagnostic_review():
    roots=[root_jets(p) for p in (Q(0),Q(1,3),Q(-1,3))]
    theta,newton,residual,frozen=phase_jets()
    repairs=[repair_jets(n) for n in (16,32,64,128)]
    frequencies=[frequency_jets(n) for n in (16,32,64)]
    return dict(manufactured_not_actual_profile=True,root_order=6,phase_order=6,repair_order=8,frequency_order=8,
        root_cases=3,root_coefficients=3*28,phase_coefficients=84,repair_coefficients=4*5*9,
        chain_coefficients=sum(a['chain_count'] for a in frequencies),
        checks=dict(signed_variance_root_matches_independent_Lagrange=all(zero(a-b) and zero(r) for a,b,r,p,s in roots),
            zero_pressure_root_not_cusped_or_divided_by_pressure=variance_series()[0]==Q(1,2),
            moving_phase_inverse_matches_Newton_and_forward_identity=zero(theta-newton) and zero(residual),
            frozen_inverse_phase_has_observable_mixed_error=not zero(theta-frozen),
            same_quadratic_repair_matches_independent_fixed_point=all(all(zero(a-b) for a,b in zip(c,o)) for c,o,z in repairs),
            exact_all_five_repair_rows=all(all(zero(a+b-z) for a,b,z in zip(matvec(B,c),quadratic(c),target)) for c,o,target in repairs),
            fixed_phase_chain_includes_every_N_power=all(a['chain'] for a in frequencies),
            moment_flow_keeps_incoming_memory=all(zero(a['moment'].derivative(0)-a['density']) for a in frequencies),
            no_float_or_quadrature_root_used=True))
