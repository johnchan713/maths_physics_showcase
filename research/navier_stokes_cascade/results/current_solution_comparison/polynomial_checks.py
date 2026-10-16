"""Exact instantaneous forced difference identities; not whole-space solutions."""
from fractions import Fraction as Q


class Poly:
    """A small rational three-variable polynomial, with no sampled derivatives."""
    def __init__(self, terms=0):
        if not isinstance(terms,dict):terms={(0,0,0):Q(terms)}
        self.terms={k:Q(v) for k,v in terms.items() if v}

    def __add__(self, other):
        other=other if isinstance(other,Poly) else Poly(other)
        d=dict(self.terms)
        for k,v in other.terms.items():d[k]=d.get(k,0)+v
        return Poly(d)

    __radd__=__add__
    def __neg__(self):return Poly({k:-v for k,v in self.terms.items()})
    def __sub__(self,other):return self+-as_poly(other)
    def __rsub__(self,other):return as_poly(other)+-self
    def __mul__(self,other):
        other=as_poly(other);d={}
        for k,v in self.terms.items():
            for l,w in other.terms.items():
                m=tuple(a+b for a,b in zip(k,l));d[m]=d.get(m,0)+v*w
        return Poly(d)
    __rmul__=__mul__
    def derivative(self,i):
        d={}
        for k,v in self.terms.items():
            if k[i]:
                l=list(k);l[i]-=1;d[tuple(l)]=v*k[i]
        return Poly(d)
    def __eq__(self,other):return self.terms==as_poly(other).terms


def as_poly(p):return p if isinstance(p,Poly) else Poly(p)
def gradient(p):return [p.derivative(i) for i in range(3)]
def divergence(v):return sum(v[i].derivative(i) for i in range(3))
def laplacian(p):return sum(p.derivative(i).derivative(i) for i in range(3))
def advect(v,w):return [sum(v[j]*w[i].derivative(j) for j in range(3)) for i in range(3)]
def double_divergence(v):return sum((v[i]*v[j]).derivative(i).derivative(j) for i in range(3) for j in range(3))
def subtract(v,w):return [a-b for a,b in zip(v,w)]


def forced_difference_review():
    x=Poly({(1,0,0):1});y=Poly({(0,1,0):1});z=Poly({(0,0,1):1})
    r2=x*x+y*y
    u=[-y,x,Poly()];v=[x*x-y*y,-2*x*y,Poly()]
    p=3*x*x+y*z;pi=-Q(1,2)*r2*r2-Q(1,2)*r2;P=p+pi
    f=[4*x,Poly(),Poly()]
    # Solve each instantaneous momentum equation for its time jet.
    def time_jet(a,b):return [laplacian(a[i])-advect(a,a)[i]-gradient(b)[i]+f[i] for i in range(3)]
    ut=time_jet(u,p);vt=time_jet(v,P);w=subtract(v,u);wt=subtract(vt,ut)
    ddg=double_divergence(v)-double_divergence(u)
    diff=[wt[i]+advect(v,w)[i]+advect(w,u)[i]-laplacian(w[i])+gradient(pi)[i] for i in range(3)]
    return dict(field='u=(-y,x,0); v=(x^2-y^2,-2xy,0); instantaneous time jets from same f=(4x,0,0)',
        pressure='p=3x^2+yz; P=p-(x^2+y^2)^2/2-(x^2+y^2)/2',
        exact_force_divergence='4',exact_difference_stress_double_divergence='8*(x^2+y^2)+2',
        no_global_finite_energy_trajectory_claim=True,
        checks=dict(both_velocities_solenoidal=divergence(u)==0 and divergence(v)==0,
            both_instantaneous_time_jets_solenoidal=divergence(ut)==0 and divergence(vt)==0,
            common_force_has_nonzero_divergence=divergence(f)==4,
            individual_forced_pressure_u=-laplacian(p)==double_divergence(u)-divergence(f),
            individual_forced_pressure_v=-laplacian(P)==double_divergence(v)-divergence(f),
            individual_unforced_pressure_shortcut_fails=-laplacian(p)!=double_divergence(u),
            force_cancels_in_difference_pressure=-laplacian(pi)==ddg,
            retaining_force_divergence_again_in_difference_is_wrong=-laplacian(pi)!=ddg-divergence(f),
            full_advective_difference_equation=all(a==0 for a in diff),
            difference_pressure_gradient_nonzero=any(a!=0 for a in gradient(pi)),
            difference_stress_double_divergence_nonzero=ddg==8*r2+2))
