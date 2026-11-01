"""Exact finite controls for the two separate log changes of variable."""
from fractions import Fraction as Q
from math import factorial

def stirling(n):
    row=[1]
    for k in range(n):
        row=[k*v+w for v,w in zip(row+[0],[0]+row)]
    return row

def step(j):
    return 162*3260**j*factorial(j)**3

def taper(n):
    if n==0:return Q(4)
    return 4*sum(Q(v*step(j),64**j) for j,v in enumerate(stirling(n)) if j)

def exponential(H):
    out=[Q(1)]
    for n in range(1,len(H)):
        out.append(sum(j*H[j]*out[n-j] for j in range(1,n+1))/n)
    return out

def physical(i,R,raw):
    return sum(v*a for v,a in zip(stirling(i),raw))*Q(1,R**i*factorial(i))

def record():
    H=[0,1]+[Q(step(j-1),factorial(j)) for j in range(2,5)]
    return {'control_order':4,'exp_coefficients':list(map(str,exponential(H))),
            'physical_taper_control':str(physical(4,2,[taper(j) for j in range(5)])),
            'physical_pulse_established':False}
