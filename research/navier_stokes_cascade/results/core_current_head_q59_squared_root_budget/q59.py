"""Finite positive majorants for a squared-root solve and a homogeneous path."""
from math import comb


def squared(k,h,t):
    if len(t)!=len(h) or not t or k<=0 or min(*t,*h)<0:
        raise ValueError("nonnegative matching jet archives and positive inverse bound required")
    y=[]
    for j,value in enumerate(t):
        y.append(k*(value+sum(h[i]*y[j-i] for i in range(1,j+1))))
    return y


def path(d,J,K,L,entrance,generator):
    if not 0<=d<=J or K<1 or L<0:
        raise ValueError("invalid triangular domain")
    u={}
    for a in range(d+1):
        u[a,0]=K*(entrance(a)+L*sum(comb(a,i)*generator(i,0)*u[a-i,0]
                                     for i in range(1,a+1)))
    for b in range(1,J+1):
        for a in range(min(d,J-b)+1):
            u[a,b]=sum(comb(a,i)*comb(b-1,j)*generator(i,j)*u[a-i,b-1-j]
                       for i in range(a+1) for j in range(b))
    return u


def record():
    from budget import plan
    out={}
    for r in (0,1):
        rows=plan(0,r)['shells']
        out[str(r)]=[dict(slow=1+x['radial'],total=1+x['radial']+x['aux']) for x in rows]
    return dict(shell_indices_at_zero=out,physical_jet_enclosures_verified=False)
