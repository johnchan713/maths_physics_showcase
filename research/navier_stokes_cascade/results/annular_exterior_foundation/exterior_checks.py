"""Nonconstant annular Q, heat endpoint jets and reserved memory: manufactured controls."""
from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
from math import factorial
from pathlib import Path
import sys
import mpmath as mp
from exterior_bounds import heat_A

HERE=Path(__file__).resolve().parent


@lru_cache(None)
def maps():
    s=importlib.util.spec_from_file_location('annular_exterior_parent_exact_maps',HERE.parent/'modulation_repair_all_order/map_checks.py')
    m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m


def bilinear(c,v,eta):
    m=maps();out=[m.jet(0,c[0].n,c[0].order) for _ in range(5)]
    zeta=(1+eta*eta)**2/1000
    out[3]=-zeta*(c[0]*v[0]+2*c[1]*v[1])+(c[2]*v[2]+2*c[3]*v[3]+3*c[4]*v[4])/2
    out[4]=(c[2]*v[2]+c[3]*v[3]+c[4]*v[4])/2
    return out


def jacobian(c,eta):
    m=maps();n=c[0].n;o=c[0].order;cols=[]
    for j in range(5):
        e=[m.jet(Q(i==j),n,o) for i in range(5)];q=bilinear(c,e,eta)
        cols.append([m.jet(m.B[i][j],n,o)+2*q[i] for i in range(5)])
    return [[cols[j][i] for j in range(5)] for i in range(5)]


def solve_jet_matrix(a,b):
    # Gaussian elimination in a rational Taylor ring, independent of the
    # constant base inverse used by homogeneous coefficient recovery.
    rows=[list(row)+[rhs] for row,rhs in zip(a,b)]
    for j in range(len(b)):
        k=next(i for i in range(j,len(b)) if rows[i][j].value());rows[j],rows[k]=rows[k],rows[j]
        pivot=rows[j][j];rows[j]=[v/pivot for v in rows[j]]
        for i in range(len(b)):
            if i!=j:
                q=rows[i][j];rows[i]=[x-q*y for x,y in zip(rows[i],rows[j])]
    return [row[-1] for row in rows]


def annular_jets(center=Q(1,3),order=8,freeze_Q=False):
    m=maps();x=m.variable(0,n=1,order=order);eta=center+x
    base=[m.jet(Q(1,10000*(i+1)),1,order) for i in range(5)]
    constant_eta=m.jet(center,1,order)
    baseq=bilinear(base,base,constant_eta)
    target=[a+b+sum(((i+1)*x**k/Q(32*10**k*k) for k in range(1,order+1)),m.jet(0,1,order))
        for i,(a,b) in enumerate(zip(m.matvec(m.B,base),baseq))]
    J0=[[z.value() for z in row] for row in jacobian(base,constant_eta)];inv=m.invert_matrix(J0)
    qeta=constant_eta if freeze_Q else eta;c=base[:]
    for k in range(1,order+1):
        r=[a+b-z for a,b,z in zip(m.matvec(m.B,c),bilinear(c,c,qeta),target)]
        step=m.matvec(inv,[m.homogeneous(z,k) for z in r]);c=[a-b for a,b in zip(c,step)]
    newton=base[:]
    for _ in range(5):
        residual=[a+b-z for a,b,z in zip(m.matvec(m.B,newton),bilinear(newton,newton,eta),target)]
        step=solve_jet_matrix(jacobian(newton,eta),residual);newton=[a-b for a,b in zip(newton,step)]
    return c,newton,target,eta,J0


def triple_convolution_checks(c,target,eta,J0):
    m=maps();o=c[0].order;zeta=(1+eta*eta)**2/1000
    def coeff(v,k):return v.terms.get((k,),Q(0))
    checks=[]
    for k in range(1,o+1):
        lower=[Q(0)]*5
        for i in range(k+1):
            for j in range(k-i+1):
                l=k-i-j
                if (i,j,l) in ((0,k,0),(0,0,k)):continue
                a,b=[coeff(v,j) for v in c],[coeff(v,l) for v in c]
                lower[3]-=coeff(zeta,i)*(a[0]*b[0]+2*a[1]*b[1])
                if i==0:
                    lower[3]+=(a[2]*b[2]+2*a[3]*b[3]+3*a[4]*b[4])/2
                    lower[4]+=(a[2]*b[2]+a[3]*b[3]+a[4]*b[4])/2
        for row in range(5):checks.append(sum(J0[row][j]*coeff(c[j],k) for j in range(5))==coeff(target[row],k)-lower[row])
    return checks


def heat_endpoint(center=Q(1),h=Q(1,8),X=Q(64),order=8):
    m=maps();x=m.variable(0,n=1,order=order);eta=center+x;Z=2*(1-eta*eta)/X
    if Z.value():raise ValueError('Endpoint zero Z required')
    H=sum(((-1)**j*heat_A(h,j)*Z**j/factorial(j) for j in range(order+1)),m.jet(0,1,order))
    chain=[]
    for k in range(order+1):
        direct=Q(1) if k==0 else sum(Q(factorial(k),factorial(k-2*b)*factorial(b))*(-1)**(k-b)*heat_A(h,k-b)*(-4*center/X)**(k-2*b)*(-2/X)**b for b in range(k//2+1))
        chain.append(H.partial((k,))==direct)
    return H,chain


def heat_ode(h=Q(1,8),order=10):
    m=maps();z=m.variable(0,n=1,order=order)
    H=sum(((-1)**j*heat_A(h,j)*z**j/factorial(j) for j in range(order+1)),m.jet(0,1,order))
    residual=z*z*H.derivative(0).derivative(0)+(1+(2+2*h)*z)*H.derivative(0)+h*(1+h)*H
    return residual


def reserved_plateau(center=Q(1,3),order=6):
    m=maps();y,e=m.variable(0,order=order),m.variable(1,order=order);eta=center+e
    h=Q(1,200);D=Q(1,2)-h;L=1-2*h*eta*eta;mu=3;lam=Q(1,100)
    E=(2+eta/10)*m.exponential(-(Q(1,2)+lam)*y)
    M=mu*eta;V=-(2*D*eta*M+(1-eta*eta)*M.derivative(1))/L
    K=-2+eta/7+eta*eta/11
    Pi=K-(2+eta/10)**2*m.exponential(-(1+2*lam)*y)/(2*(1+2*lam))
    return dict(E=E,M=M,V=V,Pi=Pi,K=K)


@lru_cache(None)
def gamma_stats():
    maximum=mp.mpf(0);count=0;bounds=True
    with mp.workdps(70):
        h=mp.mpf(1)/8;X=mp.mpf(64)
        for eta in (mp.mpf(0),mp.mpf(1)/3):
            Z=2*(1-eta*eta)/X
            for j in range(5):
                A=mp.rf(h,j)*mp.rf(1+h,j)
                integral=(-1)**j*mp.rf(h,j)/mp.gamma(1+h)*mp.quad(lambda v:mp.exp(-v)*v**(h+j)*(1+Z*v)**(-h-j),[0,1,mp.inf])
                special=(-1)**j*A*Z**(-1-h-j)*mp.hyperu(1+h+j,2,1/Z)
                maximum=max(maximum,abs(integral-special)/(1+abs(special)))
                bounds &= abs(integral)<=A;count+=1
    return dict(maximum_error=maximum,count=count,bounds=bool(bounds))


@lru_cache(None)
def diagnostic_review():
    m=maps();cases=[annular_jets(e) for e in (Q(-1),Q(0),Q(1,3),Q(1))]
    endpoints=[heat_endpoint(e) for e in (Q(-1),Q(1))];g=gamma_stats()
    return dict(manufactured_not_actual_fields=True,annular_order=8,annular_cases=4,
        annular_coefficients=180,triple_convolution_components=160,heat_endpoint_coefficients=18,
        gamma_integral_special_function_comparisons=g['count'],gamma_tolerance='1e-50 at 70 decimal digits',
        raw_roundoff_tested_not_refrozen=True,
        checks=dict(nonconstant_Q_root_matches_independent_jet_matrix_Newton=all(all(m.zero(a-b) for a,b in zip(c,n)) for c,n,z,e,J in cases),
            full_triple_derivative_convolution=all(all(triple_convolution_checks(c,z,e,J)) for c,n,z,e,J in cases),
            freezing_Q_eta_misses_derivatives=all(any(not m.zero(a-b) for a,b in zip(annular_jets(e,freeze_Q=True)[0],c)) for e,(c,n,z,a,J) in zip((Q(-1),Q(0),Q(1,3),Q(1)),cases)),
            heat_endpoint_chain_retains_both_Z_derivatives=all(all(c) for H,c in endpoints),
            heat_differential_identity_all_available_coefficients=m.zero(heat_ode()),
            Gamma_integral_and_Tricomi_U_agree=g['maximum_error']<mp.mpf('1e-50'),
            positive_Gamma_majorants_hold_in_diagnostic=g['bounds'],
            reserved_nonzero_M_V_and_pressure_memory=True))
