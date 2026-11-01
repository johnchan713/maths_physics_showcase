"""Positive disk bounds; C,F,Y are outward bounds supplied as in README."""
from fractions import Fraction as Q

def moments(R,B,P,C,Y,root2B,root2R):
    if not (0<R<=B and P>0 and C>=P and Y>=0): raise ValueError("invalid interval")
    q=B-R
    H=Q(5,4)*P*R*root2R
    return (6*B,H+2*C*root2B*q,6*H+12*C*root2B*q,
            R*(36+Q(5,3)*P**2)+(36+2*C**2)*q,10*P**2+2*C**2*Y)

def state(B,F,h,m):
    M,_,_,S,Pi=m
    D=Q(1,2)-h; A=Q(1,2)+h
    W=B+(3*D+13)*M
    N=6*W+7*D*M+(6*h+13)*S+B*(6*A+13)*Pi
    return 2009,48*F,8*F*N

def record():
    m=moments(2,2,16,16,0,2,2)
    return {'matched_trace_control':[str(x) for x in m],
            'state_control':[str(x) for x in state(2,Q(1,16),Q(1,100),m)],
            'mixed_radial_inputs_closed':False}
