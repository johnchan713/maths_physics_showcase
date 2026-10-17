"""Independent exact physical Taylor jets versus normalized-operator words.

Manufactured fields use rational h and a small Lambda, NOT the actual core.
No floating-point differentiation or sampling of the enormous solution.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import factorial


class Jet:
    def __init__(self,terms=0,n=4,order=4):
        self.n,self.order=n,order
        if not isinstance(terms,dict):terms={(0,)*n:Q(terms)}
        self.terms={k:Q(v) for k,v in terms.items() if v and sum(k)<=order}

    @classmethod
    def variable(cls,i,value,n=4,order=4):
        k=[0]*n;k[i]=1
        return cls({(0,)*n:Q(value),tuple(k):Q(1)},n,order)

    def cast(self,v):return v if isinstance(v,Jet) else Jet(v,self.n,self.order)
    def __add__(self,v):
        v=self.cast(v)
        if v.n!=self.n:raise ValueError('Different jet dimensions')
        d=dict(self.terms)
        for k,x in v.terms.items():d[k]=d.get(k,0)+x
        return Jet(d,self.n,min(self.order,v.order))
    __radd__=__add__
    def __neg__(self):return Jet({k:-v for k,v in self.terms.items()},self.n,self.order)
    def __sub__(self,v):return self+-self.cast(v)
    def __rsub__(self,v):return self.cast(v)+-self
    def __mul__(self,v):
        v=self.cast(v);order=min(self.order,v.order);d={}
        if v.n!=self.n:raise ValueError('Different jet dimensions')
        for k,a in self.terms.items():
            for l,b in v.terms.items():
                m=tuple(x+y for x,y in zip(k,l))
                if sum(m)<=order:d[m]=d.get(m,0)+a*b
        return Jet(d,self.n,order)
    __rmul__=__mul__
    def __pow__(self,p):
        p=Q(p);c=self.value()
        if p.denominator==1 and p>=0:
            out=Jet(1,self.n,self.order)
            for _ in range(int(p)):out=out*self
            return out
        if not c:raise ValueError('Nonzero base needed for fractional or inverse power')
        if p.denominator!=1 and c!=1:raise ValueError('Rational fractional jets here require unit base')
        coefficient=c**int(p) if p.denominator==1 else Q(1)
        z=self*(1/c)-1;term=Jet(1,self.n,self.order);out=term;binomial=Q(1)
        for k in range(1,self.order+1):
            term=term*z;binomial*=Q(p-k+1,k);out=out+binomial*term
        return coefficient*out
    def __truediv__(self,v):return self*self.cast(v)**-1
    def derivative(self,i):
        if not self.order:raise ValueError('Derivative exceeds available Taylor order')
        d={}
        for k,v in self.terms.items():
            if k[i]:
                l=list(k);l[i]-=1;d[tuple(l)]=v*k[i]
        return Jet(d,self.n,self.order-1)
    def value(self):return self.terms.get((0,)*self.n,Q(0))
    def partial(self,alpha):
        if len(alpha)!=self.n or min(alpha)<0 or sum(alpha)>self.order:raise ValueError('Unavailable jet order')
        return self.terms.get(tuple(alpha),Q(0))*product_factorials(alpha)


def product_factorials(alpha):
    out=1
    for a in alpha:out*=factorial(a)
    return out


def profile(Y,eta,h):
    D=Q(1,2)-h;L=1-2*h*eta*eta
    F=1+Y/10+eta/7+Y*eta*eta/13
    U=eta+Y/11+Y*Y*eta/17
    average=eta+Y/22+Y*Y*eta/51
    average_eta=1+Y*Y/51
    v0=(2*eta*U-2*D*eta*average-(1-eta*eta)*average_eta)/L
    Pi=-2+Y/19+eta*eta/23
    return F,U,v0,Pi


def implicit_root(z,t,h):
    q=Jet(1,z.n,z.order)
    for _ in range(4):
        residual=q-z*z*q**(2*h)-(1-t)
        q=q-residual/(1-2*h*z*z*q**(2*h-1))
    residual=q-z*z*q**(2*h)-(1-t)
    if residual.terms:raise ValueError('Newton jet did not solve the full implicit equation')
    return q


def physical_field(point=(Q(1,4),Q(1,5),Q(1,3),Q(1,9)),h=Q(1,200)):
    x,y,z,t=[Jet.variable(i,v) for i,v in enumerate(point)]
    q=implicit_root(z,t,h);A=Q(1,2)+h;D=Q(1,2)-h
    Y=3*(x*x+y*y)/(2*q);eta=z*q**-D
    F,U,v0,Pi=profile(Y,eta,h)
    velocity=[x*v0/(2*q)-y*q**(-A-Q(1,2))*F,
              y*v0/(2*q)+x*q**(-A-Q(1,2))*F,q**-A*U]
    return velocity,q**(-2*A)*Pi,q


def normalized_field(h=Q(1,200)):
    xi,upsilon,eta=[Jet.variable(i,v,n=3) for i,v in enumerate((Q(1,4),Q(1,5),Q(1,3)))]
    Y=Q(3,2)*(xi*xi+upsilon*upsilon);F,U,v0,Pi=profile(Y,eta,h);A=Q(1,2)+h
    return [[(Q(1,2),xi*v0/2),(A,-upsilon*F)],
            [(Q(1,2),upsilon*v0/2),(A,xi*F)],[(A,U)],[(2*A,Pi)]]


def operator_word(terms,word,h=Q(1,200),advance_beta=True,omit_time_eta=False):
    D=Q(1,2)-h
    coords=[Jet.variable(i,v,n=3) for i,v in enumerate((Q(1,4),Q(1,5),Q(1,3)))]
    xi,upsilon,eta=coords;L=1-2*h*eta*eta
    for axis in word:
        out=[]
        for beta,G in terms:
            if axis<2:H=G.derivative(axis);shift=Q(1,2)
            else:
                euler=xi*G.derivative(0)+upsilon*G.derivative(1)
                if axis==2:
                    H=((1-eta*eta)*G.derivative(2)-eta*euler-2*beta*eta*G)/L;shift=D
                elif axis==3:
                    eta_term=0 if omit_time_eta else D*eta*G.derivative(2)
                    H=(beta*G+euler/2+eta_term)/L;shift=Q(1)
                else:raise ValueError('Unknown Cartesian/time derivative')
            out.append((beta+shift if advance_beta else beta,H))
        terms=out
    return sum(G.value() for beta,G in terms)


@lru_cache(None)
def coordinate_stats():
    u,p,q=physical_field();fields=u+[p];normalized=normalized_field()
    indices=[a for a in product(range(5),repeat=4) if sum(a)<=4]
    errors=[]
    for alpha in indices:
        word=tuple(i for i,a in enumerate(alpha) for _ in range(a))
        errors.extend(fields[i].partial(alpha)-operator_word(normalized[i],word) for i in range(4))
    freeze_gap=operator_word(normalized[3],(3,3),advance_beta=False)-p.partial((0,0,0,2))
    eta_gap=operator_word(normalized[2],(3,),omit_time_eta=True)-u[2].partial((0,0,0,1))
    commute=operator_word(normalized[3],(2,3,0))-operator_word(normalized[3],(0,3,2))
    # Exact full residual of the manufactured field, retaining pressure and axial viscosity.
    basis=[tuple(int(i==j) for i in range(4)) for j in range(4)]
    values=[v.value() for v in u]
    terms=dict(time=[v.partial(basis[3]) for v in u],
        advection=[sum(values[j]*v.partial(basis[j]) for j in range(3)) for v in u],
        pressure=[p.partial(basis[j]) for j in range(3)],
        viscosity=[-sum(v.partial(tuple(2*x for x in basis[j])) for j in range(3)) for v in u])
    residual=[sum(v[i] for v in terms.values()) for i in range(3)]
    axis_u,axis_p,axis_q=physical_field((0,0,Q(1,3),Q(1,9)))
    divergence=sum(axis_u[i].partial(basis[i]) for i in range(3))
    off_axis_divergence=sum(u[i].derivative(i) for i in range(3))
    return dict(multiindices=len(indices),identities=len(errors),maximum_exact_error=max(abs(e) for e in errors),
        implicit_jet_nonzero_orders=sum(bool(k[2]+k[3]) for k in q.terms),
        frozen_exponent_gap=abs(freeze_gap),missing_time_eta_gap=abs(eta_gap),commutation_error=abs(commute),
        axis_divergence=divergence,axis_radial_derivative=axis_u[0].partial(basis[0]),
        off_axis_divergence_error=max([abs(v) for v in off_axis_divergence.terms.values()]+[Q(0)]),
        pressure_omission_gap=max(abs(v) for v in terms['pressure']),
        axial_viscosity_omission_gap=max(abs(v.partial((0,0,2,0))) for v in u),
        manufactured_full_residual=max(abs(v) for v in residual))


def diagnostic_review():
    r=coordinate_stats()
    return dict(diagnostic_h='1/200',diagnostic_Lambda=3,maximum_total_derivative_order=4,
        multiindices=r['multiindices'],exact_coordinate_identities=r['identities'],
        fixture_is_not_actual_core=True,
        checks=dict(implicit_root_has_higher_axial_and_time_jets=r['implicit_jet_nonzero_orders']>4,
            all_four_fields_match_independent_physical_Taylor_jets=r['maximum_exact_error']==0,
            normalized_operators_commute_only_with_exponent_updates=r['commutation_error']==0 and r['frozen_exponent_gap']>0,
            omitted_parameter_time_derivative_is_detected=r['missing_time_eta_gap']>0,
            Cartesian_axis_and_off_axis_divergence_are_exactly_zero=r['axis_divergence']==0 and r['off_axis_divergence_error']==0,
            axis_radial_quotient_not_arbitrarily_set_to_zero=r['axis_radial_derivative']!=0,
            pressure_gradient_and_axial_viscosity_both_nonzero=r['pressure_omission_gap']>0 and r['axial_viscosity_omission_gap']>0,
            finite_manufactured_field_not_misidentified_as_solution=r['manufactured_full_residual']>0))
