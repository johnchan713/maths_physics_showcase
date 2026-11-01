"""Finite radial factors in the real sixth-angular algebra."""
from fractions import Fraction as Q
from math import factorial,comb

def step(r):
    return 162*3260**r*factorial(r)**3

def bump(r):
    return 64 if r==0 else step(r+1)*Q(25,2)**r

def shape(r,T,J):
    h=[0]+[(Q(1,10) if k==1 else 0)+(J+6)*Q(step(k),factorial(k)*T**k) for k in range(1,r+1)]
    d=[Q(1)]
    for n in range(1,r+1):
        d.append(sum(k*h[k]*d[n-k] for k in range(1,n+1))/n)
    return d[-1]

def relative(r,b):
    return Q(1,10)**r+3*b*sum(comb(r,k)*Q(1,10)**(r-k)*bump(k) for k in range(r+1))

def record():
    return {'control_shape_order3':str(shape(3,2000,1)),
            'control_relative_order3':str(relative(3,Q(572,10**16))),
            'field_eta':6,'state_eta':5,'physical_pulse_established':False}
