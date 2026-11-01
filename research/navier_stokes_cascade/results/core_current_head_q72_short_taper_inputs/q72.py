"""Order-dependent factors in the short-taper bound."""
from fractions import Fraction as Q
from math import factorial

def slope(m,radius):
    return sum((i+1)*Q(1,radius)**(i+1) for i in range(m+1))

def switch(m,t):
    return 1+sum(162*3260**i*factorial(i)**2*Q(1,t)**i for i in range(1,m+1))

def densities(x,f,u):
    return x*u,2*x*x*f,2*x*x*f*u,x*(u*u+x*f*f),x*f*f

def record():
    return {'manufactured_slope_factor':str(slope(4,Q(1,16))),
            'manufactured_switch_factor':str(switch(4,Q(1,32))),
            'trace_reused_at_Xa_only':True,'physical_pulse_established':False}
