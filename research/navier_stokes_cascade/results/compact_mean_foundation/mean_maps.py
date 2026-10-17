"""Exact manufactured controls; continuum mapping estimates are in README.md."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product


def trim(p):
    p=list(map(F,p))
    while len(p)>1 and not p[-1]:p.pop()
    return p or [F(0)]


def add(a,b):return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
def scale(a,c):return trim([F(c)*x for x in a])
def mul(a,b):
    p=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):p[i+j]+=x*y
    return trim(p)
def power(a,n):
    out=[F(1)]
    for _ in range(n):out=mul(out,a)
    return out
def diff(a):return trim([i*x for i,x in enumerate(a)][1:])
def value(a,x):return sum(c*F(x)**i for i,c in enumerate(a))
def primitive(a):
    p=[F(0)]+[c/F(i+1) for i,c in enumerate(a)]
    p[0]=-value(p,1)
    return trim(p)
def integral(a,weight=0):return sum(c*(F(2)**(i+weight+1)-1)/F(i+weight+1) for i,c in enumerate(a))
def divide_R(a):
    if a[0]:raise ValueError('Manufactured polynomial must vanish at R=0')
    return trim(a[1:])


def solve(a,b):
    n=len(b)
    if len(a)!=n or any(len(row)!=n for row in a):raise ValueError('Square matrix required')
    q=[[F(x) for x in row]+[F(t)] for row,t in zip(a,b)]
    for i in range(n):
        pivot=next((j for j in range(i,n) if q[j][i]),None)
        if pivot is None:raise ValueError('Singular moment block')
        q[i],q[pivot]=q[pivot],q[i];v=q[i][i];q[i]=[x/v for x in q[i]]
        for j in range(n):
            if j!=i:
                v=q[j][i];q[j]=[x-v*y for x,y in zip(q[j],q[i])]
    return [row[-1] for row in q]


def vandermonde(nodes):return [[x**j for j in range(len(nodes))] for x in nodes]
def determinant(a):
    if len(a)==2:return a[0][0]*a[1][1]-a[0][1]*a[1][0]
    if len(a)!=3:raise ValueError('Two or three rows required')
    return sum((-1)**j*a[0][j]*determinant([[r[k] for k in range(3) if k!=j] for r in a[1:]]) for j in range(3))


def blocks(g,q,a=1,moments=(1,1,1,1,1)):
    """q=exp(-2*d*lambda), g=exp(d); positive row moments supplied exactly."""
    g,q,a=map(F,(g,q,a))
    if g<=1 or not 0<q<=1 or not a or len(moments)!=5 or min(moments)<=0:raise ValueError('Invalid fixed patch data')
    angular=vandermonde([g*g,q/(g*g),q]);axial=vandermonde([g,g*q])
    return [[v*c for v in row] for row,c in zip(angular,(moments[0],2*a*moments[1],-a*moments[2]))],[[v*c for v in row] for row,c in zip(axial,(moments[3],a*moments[4]))]


def flat_order(alpha,target,loss,kappa=F(1,100000)):
    alpha,target,loss,kappa=map(F,(alpha,target,loss,kappa))
    if kappa<=0 or loss<0:raise ValueError('Positive kappa and nonnegative fixed loss required')
    # Strict margin leaves logarithmic powers absorbable on the SAME domain.
    z=(target+loss-alpha)/kappa
    return max(1,z.numerator//z.denominator+1)


@lru_cache(None)
def matrix_checks():
    count=0;rows=True;limits=True;conditioning=True
    for g,q,a in product((F(3,2),F(2)),(F(1,2),F(3,4),F(15,16),F(255,256)),(F(1,3),F(2),F(-2))):
        moments=tuple(F(i+2,i+1) for i in range(5));At,Az=blocks(g,q,a,moments)
        for P,Jt,Jz in ((F(1),F(-2),F(3)),(F(-7),F(5),F(-1))):
            targets=([0,-P,-Jz],[0,-Jt]);sol=[solve(At,targets[0]),solve(Az,targets[1])]
            for matrix,b,x in zip((At,Az),targets,sol):
                rows &= [sum(v*c for v,c in zip(row,x)) for row in matrix]==list(map(F,b));count+=len(b)
            # Zero flux fixes s0=-g*s1, independently of the target's sign.
            s0,s1=sol[1]
            conditioning &= s0==-g*s1 and (1-q)*s1==Jt/(a*moments[4]*g)
        nodes=[g*g,q/(g*g),q];V=vandermonde(nodes)
        rows &= determinant(V)==(nodes[1]-nodes[0])*(nodes[2]-nodes[0])*(nodes[2]-nodes[1])
        limits &= determinant(blocks(g,1,a,moments)[0])!=0 and determinant(blocks(g,1,a,moments)[1])==0
    return dict(moment_rows=count,rows=bool(rows),angular_limit_nonsingular=bool(limits),sharp_axial_scaled_identity=bool(conditioning))


@lru_cache(None)
def pressure_checks():
    bump=mul(power([-1,1],4),power([2,-1],4));rho0=scale(bump,1/integral(bump));rhoZ=mul(rho0,[-F(3,2),1])
    Z=F(1,3);rho=add(rho0,scale(rhoZ,Z));g0=mul(bump,[2,-1,1]);gZ=mul(bump,[-1,3]);g=add(g0,scale(gZ,Z))
    P,PZ=integral(g),integral(gZ)
    p=primitive(add(g,scale(rho,-P)))
    pZ=primitive(add(gZ,add(scale(rhoZ,-P),scale(rho,-PZ))))
    c=integral(rho,2)/2;cZ=integral(rhoZ,2)/2
    moving=integral(pZ,1)==-integral(gZ,2)/2+c*PZ+cZ*P
    # All six entries of a real two-quadrature tensor, including axial ones.
    b=mul([0,1],bump);aa=[mul(b,v) for v in ([Z,1],[-Z,2],[2*Z,0,1])];az=[b,scale(b,-1),scale(b,2)]
    bb=[mul(b,v) for v in ([1,1],[2,0,1],[-1,2])];k=F(3)
    W={};WZ={}
    for i in range(3):
        for j in range(i,3):
            W[i,j]=scale(add(mul(aa[i],aa[j]),scale(mul(bb[i],bb[j]),1/(k*k))),F(1,2))
            WZ[i,j]=scale(add(mul(az[i],aa[j]),mul(aa[i],az[j])),F(1,2))
    eps=F(1,7);gr=add(scale(diff(W[0,0]),-1),add(divide_R(add(W[1,1],scale(W[0,0],-1))),scale(WZ[0,2],-eps)))
    PP=integral(gr);JJz=integral(W[2,2],1)-integral(gr,2)/2;JJt=integral(W[1,2],2)
    expectedP=integral(divide_R(add(W[1,1],scale(W[0,0],-1))))-eps*integral(WZ[0,2])
    expectedJ=integral(add(W[2,2],scale(add(W[0,0],W[1,1]),-F(1,2))),1)+eps*integral(WZ[0,2],2)/2
    return dict(unit_bump=integral(rho)==1 and integral(rhoZ)==0,compact_pressure=value(p,1)==value(p,2)==0,
        pressure_moment=integral(p,1)==-integral(g,2)/2+c*P,moving_pressure_moment=moving,
        moving_product_omission_nonzero=bool(cZ*P),complete_initial_defects=PP==expectedP and JJz==expectedJ,
        axial_covariances_nonzero=bool(JJt) and bool(integral(W[2,2],1)),
        radial_axial_flux_omission_nonzero=bool(eps*integral(WZ[0,2])),tensor_entries=len(W),
        pressure_moment_omission_nonzero=bool(integral(gr,2)/2))


@lru_cache(None)
def inverse_checks():
    norms=0;normok=True
    for k1,k2 in product(range(-8,9),repeat=2):
        if not (k1 or k2):continue
        for A,B in ((k1+k2,-k2),(-k1+k2,k1)):
            n=A*A-2*B*B;normok &= type(n) is int and abs(n)>=1;norms+=1
    count=0;path=True;missing=False
    source=mul(power([-1,1],4),power([2,-1],4))
    for tau in (F(-7),F(2),F(13,3)):
        H=[F(0)];d=source[:]
        for j in range(len(source)):
            H=add(H,scale(d,F((-1)**j)/tau**(j+1)));d=diff(d)
        for order in range(7):
            path &= add(diff(H),scale(H,tau))==source;missing |= diff(H)!=source;count+=1
            H,source=diff(H),diff(source)
        source=mul(power([-1,1],4),power([2,-1],4))
    choices=[(a,t,l) for a,t,l in product((F(-2),F(1),F(3,2)),(F(0),F(4),F(11)),(F(0),F(5,2)))]
    margins=[a+flat_order(a,t,l)*F(1,100000)-l-t for a,t,l in choices]
    return dict(quadratic_integer_norm_controls=norms,divisors_nonzero=bool(normok),
        shifted_primitive_jet_controls=count,shifted_path_identity=bool(path),omitting_shift_detected=bool(missing),
        strict_flatness_choices=len(margins),strict_flatness_margins=all(x>0 for x in margins))


def normalization_checks():
    count=0;rows=True;swirl=True
    for h,lam in product((F(1,200),F(1,1000)),(F(1,5),F(1,8),F(1,10000))):
        A=F(1,2)+h
        # The dR contribution is 1/2, each velocity -A, each R power 1/2.
        physical=[F(3,2)-A,1-A,-2*A,F(3,2)-2*A,1-2*A]
        row_factors=[A-F(3,2),A-1,2*A,2*A-F(3,2),2*A-1]
        rows &= all(x+y==0 for x,y in zip(physical,row_factors));count+=5
        swirl &= -A+(1+2*lam)/2==lam-h
    return dict(normalization_rows=count,all_five_rows=bool(rows),nonconstant_s_swirl_power=bool(swirl))


def reviews():
    return dict(exact_signed_five_map=matrix_checks(),complete_pressure_and_covariance_moments=pressure_checks(),
        exact_divisor_path_and_common_domain_controls=inverse_checks(),physical_row_normalizations=normalization_checks())
