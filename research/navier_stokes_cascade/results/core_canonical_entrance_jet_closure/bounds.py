"""Nonzero canonical-entrance jets from the retained phase/propagator budgets."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb,factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,HERE.parent/path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
reference=load('canonical_entrance_reference','core_reference_scaled_normal_floor/bounds.py')
inverse=reference.transfer.inverse
P,exact,order,product=inverse.P,inverse.exact,inverse.order,inverse.product

class CanonicalEntranceBudget(inverse.ForcedBudget):
    """Bounds only: no ODE solution, source or Gaussian floor is materialized.

    tangential_inverse_square bounds 1/|n_tan|^2, not just 1/|n|^2.
    The production factory derives both from the same original certificate.
    """
    def __init__(self,base,*,tangential_inverse_square,cmax,u):
        q,c,u=map(exact,(tangential_inverse_square,cmax,u))
        if q<=0 or c<=0 or u<0 or base.harmonic!=1:
            raise ValueError('Positive tangential floor/cmax, fixed nonnegative tilt and primary |m|=1 required')
        self.entry_q,self.entry_c,self.entry_u=q,c,u
        super().__init__(normal=base.normal,moving_normal=base.moving_normal,matrix=base.matrix,
            source=lambda n:P(),profile_rate=base.profile_rate,tangent=base.tangent,
            length=base.length,inverse_normal_square=base.inverse_normal_square,harmonic=1)

    @lru_cache(None)
    def entry_reciprocal(self,n):
        order(n)
        if not n:return P(self.entry_q)
        # |D^I |n_tan|^2| <= (normal * normal)_|I|; no equality
        # between full and tangential squares is asserted.
        return self.entry_q*sum((comb(n,i)*self.square(i)*self.entry_reciprocal(n-i)
                                for i in range(1,n+1)),P())

    @lru_cache(None)
    def entry_inverse_sqrt(self,n):
        order(n)
        if not n:return P(max(Q(1),self.entry_q))
        # D(J^-1/2)=-(1/2)J^-1/2 J^-1 DJ; raw mixed Leibniz.
        return Q(1,2)*product(self.entry_inverse_sqrt,
            lambda j:product(self.entry_reciprocal,lambda k:self.square(k+1),j),n-1)

    @lru_cache(None)
    def entry_vector(self,n):
        order(n)
        quadratic=product(lambda j:self.data('normal',j),
            lambda k:product(lambda i:self.data('normal',i),self.entry_reciprocal,k),n)
        rotated=product(lambda j:self.data('normal',j),self.entry_inverse_sqrt,n)
        C=self.entry_c*(1+self.entry_u/2)
        return P(int(n==0))+quadratic+C*rotated

    def entrance_log(self,n):
        order(n)
        return P() if not n else Q(1,2)*self.length*self.data('profile_rate',n)

    @lru_cache(None)
    def entrance_gaussian(self,n):
        order(n)
        if not n:return P(1)
        return sum((comb(n-1,i-1)*self.entrance_log(i)*self.entrance_gaussian(n-i)
                    for i in range(1,n+1)),P())

    @lru_cache(None)
    def initial(self,n):
        order(n)
        return product(self.entrance_gaussian,self.entry_vector,n)

    @lru_cache(None)
    def velocity(self,n):
        order(n)
        parameter=self.ambient()*(self.initial(n)+self.length*sum(
            (comb(n,i)*self.generator(i)*self.velocity(n-i) for i in range(1,n+1)),P()))
        if not n:return parameter
        pulse=sum((comb(n-1,i)*self.generator(i)*self.velocity(n-1-i)
                   for i in range(n)),P())
        return parameter+pulse

    def wave(self,n,cutoff,dimension=5):
        """RootBudget wave slot: J_n(cutoff*t)/P, not J_n(t/P)."""
        order(n)
        if type(dimension) is not int or dimension<1:raise ValueError('Positive fixed coefficient-chart dimension required')
        return Q(dimension**n,factorial(n))*product(cutoff,self.velocity,n)

def from_original(p,selector_inputs,ell,*,F,G,btilde,path):
    """Same original accepted tail and canonical entrance, no new W_j input."""
    phase,base,old,c,tc=reference.forced_from_reference(p,selector_inputs,ell,
        F=F,G=G,btilde=btilde,path=path,source=lambda n:P(),harmonic=1)
    h=CanonicalEntranceBudget(base,tangential_inverse_square=c['inverse_normal_square'],
        cmax=p['cmax'],u=p['u'])
    return phase,h,c,tc

def root_wave_join(template,homogeneous,cutoff,dimension=5):
    """Replace a NORM slot only; no physical canonical root is rebuilt."""
    root=inverse.roots
    if not isinstance(template,root.RootBudget) or not isinstance(homogeneous,CanonicalEntranceBudget):
        raise TypeError('Original root-budget and derived canonical-entrance recipe required')
    return root.RootBudget(H=template.H,target=template.target,
        wave=lambda n:homogeneous.wave(n,cutoff,dimension),curl=template.curl,
        loghalf=template.loghalf,K0=template.K0,beta=template.beta,R=template.R,
        mass=template.mass,width=template.width)

def certificate():
    return dict(parent_commit='9e8f33555745e8a5a07eeeeadb5c10d8798c5a74',
        result='explicit canonical homogeneous mixed-jet/pressure/potential/root-wave recipe from retained full phase bounds and original nonzero entrance',
        entrance='Br(0)eplus=er-nr*ntan/|ntan|^2+c0*wr(0)*rot(ntan)/|ntan|; wr(0)=sqrt(1+u^2/4), frozen c0,u',
        gaussian='D^j[P(0)b0]/P(0) uses Bell/Leibniz products; log P(0)=-integral_0^(L/2) p, L frozen',
        recurrence='U0=Kamb I0; Uj=Kamb[Ij+L sum_(i=1)^j binom(j,i) Ahat_i U_(j-i)]+sum_(i=0)^(j-1) binom(j-1,i) Ahat_i U_(j-1-i)',
        join='Wj=(dimension^j/j!)*(cutoff*U)_j; original RootBudget H,target,curl,loghalf,value gaps and roots unchanged',
        pressure='homogeneous source is exactly zero; moving-normal pressure retained. Forced blocks still retain unprojected n.f',
        canonical_and_inherited_primaries_identified=False,
        actual_background_jet_suprema_full_covariance_residual_or_blowup_certified=False)
