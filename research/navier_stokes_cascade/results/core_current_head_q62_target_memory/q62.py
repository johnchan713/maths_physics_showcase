"""Normalized finite-jet pieces for the unstripped forward stress."""
from fractions import Fraction


def heat(l,h,B):
    if l < 0 or not 0 < h < 1 or B <= 0:
        raise ValueError("invalid heat trace")
    a = 1
    c = 2*h*(1+h)
    for j in range(1,l+1):
        a *= (h+j-1)*(h+j)*Fraction(1,j)
        c += a*6**j/Fraction(B)**(j-1)
    return c, 1+c/B


def primitive(i,j,trace,length,density):
    if min(i,j)<0:
        raise ValueError("negative jet index")
    return (trace+length*density(0,j) if not i else
            density(i-1,j)*Fraction(1,i))


def stress(X,r,h,t,E,U,Ex,Ux,M,I,J,S,P,Me,Ie,Je,Se,Pe):
    D,A=Fraction(1,2)-h,Fraction(1,2)+h
    L,d=1-2*h*t*t,1-t*t
    w=X-2*D*t*M-d*Me
    v=(1-h)*I-D*t*Ie-d*Je+2*(h-D)*t*J
    return (((-E*w+v/r)/L-E+2*X*Ex)/r,
            ((-w*U+D*(M-t*Me)+4*h*t*S-d*Se+X*(4*A*t*P-d*Pe))/L+2*X*Ux)/r)


def record():
    return {"heat_control":[str(v) for v in heat(3,Fraction(1,100),2)],
            "target_trace_order_offset":2,"local_velocity_order_offset":99998,
            "actual_velocity_enclosed":False}
