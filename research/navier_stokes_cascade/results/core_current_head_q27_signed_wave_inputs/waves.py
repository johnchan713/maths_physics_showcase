"""First real signed wave: explicit four-input adapter, bounds not fields."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,HERE.parent/path)
    obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj);return obj

roots=load('q27_roots','core_primary_root_budget/budgets.py')
prefix=load('q27_prefix','core_signed_prefix_source_budget/budgets.py')
P=roots.P

def transfer(value):
    """Historical helpers own distinct classes; preserve exact coefficients."""
    return P(value.c)

def rebind_wave(template,wave):
    return roots.RootBudget(H=template.H,target=template.target,wave=wave,
        curl=template.curl,loghalf=template.loghalf,K0=template.K0,
        beta=template.beta,R=template.R,mass=template.mass,width=template.width,
        trace=template.trace)

class SignedWaveInputs:
    """Same canonical roots, two wave NORM slots for regular and cut grammar.

    xi supplies normalized Taylor-sum jets of the actual Haar signed target.
    phase supplies raw component/entry-sum bounds, including moving normal.
    No updated positive root or new wave-jet callback is accepted.
    """
    def __init__(self,template,*,cut_wave,xi,xi_pole,phase,inverse_normal_square):
        if not isinstance(template,roots.RootBudget):raise TypeError('Original RootBudget required')
        if type(xi_pole) is not int or xi_pole<0:raise ValueError('Nonnegative integer finite target pole required')
        self.q=roots.exact(inverse_normal_square)
        if self.q<=0:raise ValueError('Certified positive full-normal reciprocal bound required')
        self.regular=template;self.cut=rebind_wave(template,cut_wave)
        self.xi=xi;self.xi_pole=xi_pole;self.phase=phase

    def phase_jet(self,name,j):return transfer(getattr(self.phase,name)(j))

    def product(self,a,b,j):
        return sum((comb(j,i)*a(i)*b(j-i) for i in range(j+1)),P())

    @lru_cache(None)
    def normal_square(self,j):
        return self.product(lambda i:self.phase_jet('normal',i),lambda i:self.phase_jet('normal',i),j)

    @lru_cache(None)
    def reciprocal_square(self,j):
        if not j:return P(self.q)
        return self.q*sum((comb(j,i)*self.normal_square(i)*self.reciprocal_square(j-i)
                           for i in range(1,j+1)),P())

    def pressure_numerator(self,j):
        return self.product(lambda i:self.phase_jet('normal',i),
                            lambda i:self.phase_jet('matrix',i),j)+self.phase_jet('moving_normal',j)

    def pressure_multiplier(self,j):
        return self.product(self.reciprocal_square,self.pressure_numerator,j)

    @lru_cache(None)
    def principal(self,j,cut=False):
        root=self.cut if cut else self.regular
        return root.amplitude(j)+root.signed(j,lambda i:transfer(self.xi(i)))

    @lru_cache(None)
    def curl(self,j,cut=False):
        root=self.cut if cut else self.regular
        return sum((root.data('curl',i)*(self.principal(j-i,cut)+(j-i+1)*self.principal(j-i+1,cut))
                    for i in range(j+1)),P())

    @lru_cache(None)
    def through(self,n):
        roots.order(n)
        # Two canonical signs. A real cos/sin contributes l1 mass one
        # over its +/- Fourier pair, so no additional conjugate factor.
        T=2*sum((self.principal(j) for j in range(n+1)),P())
        R=2*sum((self.curl(j) for j in range(n+1)),P())
        U=2*sum((self.principal(j,True)+self.curl(j,True) for j in range(n+1)),P())
        multiplier=sum((Q(5**j,factorial(j))*self.pressure_multiplier(j) for j in range(n+1)),P())
        return {'T':T,'R':R,'Pi':multiplier*T,'U':U}

    def source_wave(self,n):
        return {k:prefix.P(v.c) for k,v in self.through(n).items()}

    def pole(self,n):
        roots.order(n)
        return max(6*n+9,6*n+self.xi_pole+6)

def record():
    return {'result':'conditional-first-signed-complete-wave-input-closure',
        'wave_orders':{'T':'1/2','R':'1-kappa','Pi':'1','U':'1/2'},
        'wave_weight':'sqrt(zeta)*P',
        'common_pole':'max(6*n+9,6*n+d+6)',
        'root_and_target_order_for_wave_slot':'n+1',
        'prefix_wave_order_for_source_j':'j+15',
        'root_and_target_order_for_source_j':'j+16 (before expanding target dependencies)',
        'pressure_row':'(n^T K-nprime^T)/|n|^2',
        'signed_square_and_cutoff_curl_retained':True,
        'independent_wave_T_R_Pi_U_callbacks_removed':True,
        'actual_lower_premises_certified':False,
        'complete_first_source_or_residual_proved':False}
