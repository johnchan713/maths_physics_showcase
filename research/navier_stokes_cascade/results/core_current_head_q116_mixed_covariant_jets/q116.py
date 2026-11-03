"""All fixed mixed jets: conditional tangent recurrence, not native certification."""
from fractions import Fraction as Q
from math import comb
from itertools import product as cartesian
import importlib.util
from pathlib import Path
p=Path(__file__).resolve().parent.parent/'core_current_head_q115_covariant_first_jet/q115.py'
s=importlib.util.spec_from_file_location('q116_first',p);first=importlib.util.module_from_spec(s);s.loader.exec_module(first)

def jets(a):
    a=list(a)
    if not a or any(isinstance(x,(bool,float)) or Q(x)<0 for x in a):
        raise ValueError('Nonempty exact nonnegative raw jet bounds required')
    return list(map(Q,a))

def convolution(a,b):
    if len(a)!=len(b) or not a:raise ValueError('Same nonempty jet ceiling required')
    return [sum(comb(j,i)*a[i]*b[j-i] for i in range(j+1)) for j in range(len(a))]

def coefficients(G,n,np,K,p,d,*,q0,T,L,dmin,entrance_flat_through):
    arrays=list(map(jets,(G,n,np,K,p,d)));G,n,np,K,p,d=arrays
    if any(len(a)!=len(G) for a in arrays):raise ValueError('One common finite mixed-jet ceiling required')
    if any(isinstance(x,(bool,float)) for x in (q0,T,L,dmin)):raise ValueError('Exact fixed geometry required')
    q0,T,L,dmin=map(Q,(q0,T,L,dmin));N=len(G)-1
    if min(q0,L,dmin)<=0 or T<1:raise ValueError('Positive original floors/length and tangent bound >=1 required')
    if type(entrance_flat_through) is not int or entrance_flat_through<max(-1,N-1):
        raise ValueError('Declare whole mixed source entrance flatness through N-1; norms do not prove it')
    square=convolution(n,n);inv=[q0]
    for j in range(1,N+1):
        inv.append(q0*sum(comb(j,i)*square[i]*inv[j-i] for i in range(1,j+1)))
    h=convolution(n,inv);Pi=[Q(1)]+convolution(h,n)[1:]
    A=[x+y for x,y in zip(convolution(Pi,K),convolution(h,np))]
    Qlift=(K[0]*n[0]+np[0])*q0
    source=convolution(Pi,G);C=L+Q(4,3)/dmin
    V=[];Z=[]
    for j in range(N+1):
        z=sum(comb(j,i)*n[i]*V[j-i] for i in range(1,j+1))
        lower=sum(comb(j,i)*(A[i]+p[i])*V[j-i] for i in range(1,j+1))
        damping=sum(comb(j,i)*d[i]*V[j-i] for i in range(1,j+1))
        V.append(T*(source[j]+L*Qlift*z+L*lower+C*damping)+h[0]*z);Z.append(z)
    return dict(velocity_over_residence=V,constraint_over_residence=Z,
                inverse_square_jets=inv,h_jets=h,projection_jets=Pi,A0_jets=A,
                source=G,normal=n,normal_prime=np,matrix=K,
                length=L,damping_floor=dmin,m2_residence_upper=C,
                entrance_flatness_declared_through=entrance_flat_through,
                pressure_zero_extension_flatness_sufficient=entrance_flat_through>=N,
                actual_native_or_entrance_contracts_certified=False)

def pressure(G,n,np,K,inv,U,m):
    m=first.harmonic(m)
    arrays=list(map(jets,(G,n,np,K,inv,U)))
    if any(len(a)!=len(arrays[0]) for a in arrays):raise ValueError('Same pressure jet ceiling required')
    G,n,np,K,inv,U=arrays
    row=[a+b for a,b in zip(convolution(n,K),np)]
    numerator=[a+b for a,b in zip(convolution(row,U),convolution(n,G))]
    return [x/m for x in convolution(inv,numerator)]

def specialize(out,m):
    m=first.harmonic(m);ell=first.residence(m,out['length'],out['damping_floor'])
    U=[ell*v for v in out['velocity_over_residence']]
    return dict(velocity=U,residence=ell,
                pressure_times_carrier=pressure(out['source'],out['normal'],out['normal_prime'],out['matrix'],out['inverse_square_jets'],U,m),
                potential_times_carrier=[x/m for x in convolution(out['h_jets'],U)],
                actual_source_or_high_jet_contracts_certified=False)

def partitions(alpha):
    if not alpha or any(type(j) is not int or j<0 for j in alpha):raise ValueError('Nonnegative mixed multiindex required')
    for beta in cartesian(*(range(j+1) for j in alpha)):
        gamma=tuple(a-b for a,b in zip(alpha,beta));c=1
        for a,b in zip(alpha,beta):c*=comb(a,b)
        yield beta,gamma,c

def projected_forcing(alpha,n,u,A,d,p,Pi,g,Qlift,m):
    """Exact matrix identity on supplied raw multi-jets; diagnostic backend only."""
    first.harmonic(m);zero=(0,)*len(alpha);I=A[zero]*0
    for j in range(3):I[j,j]=1
    Z=0;comm=n[zero]*0;source=n[zero]*0
    for beta,gamma,c in partitions(alpha):
        source+=c*Pi[beta]*g[gamma]
        if any(beta):
            Z+=c*(n[beta].T*u[gamma])[0]
            comm+=c*(A[beta]-(m*m*d[beta]+p[beta])*I)*u[gamma]
    return dict(forcing=Qlift*Z+Pi[zero]*comm-Pi[zero]*source,constraint=Z)

def record():
    return dict(status='all-fixed-mixed-covariant-jet-residence-recurrence',
                velocity='U_j <= ell_m V_j; V_j independent of m for uniform native inputs',
                source_jets_required='through N; entrance flat through N-1 for velocity',
                pressure_zero_extension_requires_source_flat_through_N=True,
                full_normal_pressure_and_mixed_moving_normal_retained=True,
                actual_high_jet_inputs_or_radial_gain_certified=False,
                full_cone_pulse_residual_or_blowup=False)
