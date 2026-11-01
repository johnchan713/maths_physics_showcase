"""Rational envelopes; E is an outward upper bound for exp(G)."""
from fractions import Fraction as Q


def radii(G,E,d):
    if G<2 or E<1 or not 0<d<=1:
        raise ValueError("invalid cap; caller must also establish E>=exp(G)")
    G,E,d=map(Q,(G,E,d))
    A=4*E;C=3*E*A;b=1/(2*E**4)
    R=b/(16*A*C*C)
    tau=min(Q(1),R/(4*G+2));M=2*d*(G+1)*C;m=d/(32*E**4)
    a=min(tau/4,m*tau*tau/(64*M));P=4*M/tau
    sig=min(a,m*a/(4*(1+P)))
    return dict(A=A,C=C,b=b,R=R,tau=tau,M=M,m=m,a=a,P=P,sigma=sig)


def record():
    r=radii(2,9,Q(1,32)) # exp(2)<9, using e<3.
    return {"control_radius":str(r['sigma']),"physical_inputs_enclosed":False,
            "divides_by_pressure_or_variance":False}
