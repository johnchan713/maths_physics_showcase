"""Finite jet calculus; input bounds are hypotheses, not physical certificates."""
from fractions import Fraction as Q
from math import comb

def jets(x):
    x=list(x)
    if not x or any(isinstance(v,(float,bool)) or Q(v)<0 for v in x):
        raise ValueError('nonempty exact nonnegative bounds required')
    return list(map(Q,x))

def product(a,b):
    a,b=jets(a),jets(b)
    if len(a)!=len(b):raise ValueError('same jet ceiling required')
    return [sum(comb(j,i)*a[i]*b[j-i] for i in range(j+1)) for j in range(len(a))]

def pressure(a,b,u,g,omega):
    # a=(K^T n-nprime)/|n|^2, b=n/|n|^2, u=t/P, g=f/P.
    omega=Q(omega)
    if omega==0:raise ValueError('nonzero harmonic frequency required')
    x,y=product(a,u),product(b,g)
    if len(x)!=len(y):raise ValueError('same jet ceiling required')
    return [(v+w)/abs(omega) for v,w in zip(x,y)]

def entrance_jets(H,g):
    # Exact scalar model u'=H*u-g, u(0)=0, raw derivatives, not Taylor coefficients.
    if len(H)!=len(g) or not g:raise ValueError('matching finite jets required')
    u=[0]
    for j in range(len(g)):
        u.append(sum(comb(j,i)*H[i]*u[j-i] for i in range(j+1))-g[j])
    return u

def record():
    return dict(status='first-real-full-source-entrance-and-normal-reduction',
                finite_entrance_equivalence='f[0:r]=0 iff t[0:r+1]=0 and pi[0:r]=0',
                normal_source_pressure_retained=True,extra_source_jets_for_pressure=0,
                normal_jet_ceiling='r+1',raw_exit_vanishes_from_flat_source=False,
                actual_extended_source_inequality_certified=False,
                complete_residual_improvement_or_blowup=False)
