"""Exact finite controls for the post-taper source and switch bounds."""
from fractions import Fraction as Q
from math import comb,factorial

def switch(i,delta):
    return 162*3260**i*factorial(i)**3/delta**i if i else 1

def slopes(r,delta,Xi,Vp,Vn):
    t=[switch(i,delta) for i in range(r+1)]
    return (delta*10*Vp*sum(comb(r,i)*t[i] for i in range(r+1))+Q(4,5)*t[r],
            delta*19*Xi*Vn*sum(comb(r,i)*t[i]*2**(r-i) for i in range(r+1))/2)

def record():
    return {'manufactured_r3_bounds':list(map(str,slopes(3,Q(1,16),110,1,1))),
            'source_powers':[-1,0,1],'actual_moments_are_reference':False}
