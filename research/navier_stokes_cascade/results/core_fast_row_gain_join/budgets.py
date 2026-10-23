"""Fast-time cancellation, complete pressure differences, and a finite A1 recipe."""
from fractions import Fraction as Q
from functools import lru_cache
from math import ceil,comb,factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


scalar=load('fast_row_scalar',HERE.parent/'core_actual_scalar_gain_join/budgets.py')
fresh,prefix=scalar.fresh,scalar.prefix
P,order,exact=prefix.P,prefix.order,prefix.exact


def difference(n,old,inc,fixed):
    """Full row change EXCEPT c Nt of the tangential increment.

    Fields through n+3; same wave/base/rho, including the full T0 pressure.
    For fast fields that term cancels the old Q1 rows up to exact A1 fz.
    For a slow increment it is identically zero. Output order B+1-2k.
    The radial source/pressure themselves are bounded only at order B.
    """
    order(n)
    if set(old)!=set(prefix.STATE) or set(inc)!=set(prefix.STATE):raise ValueError('Three actual mean components required')
    x,u,z=(old[k] for k in prefix.STATE);X,U,Z=(inc[k] for k in prefix.STATE)
    b,V,G,g,c,r,rho,kp,km=(fixed[k] for k in ('b','V','G','graph','c','r','rho','Kp','Kmom'))
    d=Q(n+4);dr=g*d;time=d*(1+c);vis=d*d*(g*g+1)+r*dr+r*r
    S=2*b*X+2*x*X+X*X
    H=b*Z+G*X+x*Z+X*z+X*Z
    T=2*V*U+2*u*U+U*U
    gr=time*X+(dr+r)*S+d*H+r*T+vis*X
    defect=km*gr;pm=kp*(gr+rho*defect)
    rt=b*U+V*X+x*U+X*u+X*U
    at=G*U+V*Z+z*U+Z*u+Z*U
    az=2*G*Z+2*z*Z+Z*Z+pm
    et=d*U+(dr+2*r)*rt+d*at+vis*U
    ez=d*Z+(dr+r)*H+d*az+vis*Z
    return dict(gr=gr,P=defect,pm=pm,Etheta=et,Ez=ez,viscosity=vis*(U+Z),
                pressure_gradient=d*pm,input_order=n+3)


def covariance_rows(n,D,fixed):
    """All complete covariance entries and their reconstructed pressure."""
    order(n);d=Q(n+4);dr=fixed['graph']*d;r=fixed['r']
    gr=(dr+d+2*r)*D
    pp=fixed['Kmom']*gr;pm=fixed['Kp']*(gr+fixed['rho']*pp)
    return dict(gr=gr,P=pp,pm=pm,Etheta=(dr+d+2*r)*D,Ez=(dr+d+r)*D+d*pm,input_order=n+3)


def rates(kappa=Q(1,100000)):
    k=exact(kappa)
    if not 0<k<Q(1,4):raise ValueError('Retained 0<kappa<1/4 required')
    B=1-k;H=B+1-2*k;C=B+Q(1,2);D=C-k
    q=max(1,(H-B)//k+1)
    groups=dict(fast_slow_time=B+1,fast_radial_transport=B+1-k,
                fast_viscosity=H,fast_mean_products=2*B+1-k,
                fast_pressure=B+1,old_slow_change=H,
                sourced_covariance=C-k,retained_A1=B+q*k)
    return dict(B=B,H=H,C=C,D=D,q=q,groups=groups,postwave=min(groups.values()),
                old_row_margin=H-D,A1_margin=B+q*k-H,
                fresh_fast=D,fresh_slow=C,feedback=min(D,C),next_regular_source=D,
                conditional_on_original_actual_input_and_map_contracts=True,
                complete_residual_or_unconditional_actual_inputs_certified=False)


def flat_recipe(n,kappa=Q(1,100000)):
    """An explicit finite polynomial expression; enormous jets are NOT evaluated.

    B_n=N_n(k_e chi_U), A_J=N_J(f_e) are the ORIGINAL raw kernel
    profile seminorm recipes in core_fast_inverse_budget (8). They must
    not be inferred from an arbitrary operator norm. Q1 costs two; 6^n
    safely converts raw jets to the five-variable factorial through norm.
    """
    order(n);r=rates(kappa);q=r['q'];J=n+2*q+3
    return dict(output_order=n,integration_order=q,raw_source_order=J,
                scalar=2*6**n*8,two_power=2*n+J+q+3,three_power=q,
                source_factorial=J,
                original_factors=['L','N_n(k_e chi_U)','N_J(f_e)','before.rows(J).Ez',
                                  'Eedge(2p,width)','inverse edge floor on the original chi_U hull'],
                floor_factor='(1+before.time_inverse(0))^ceil(original_rho_g*q)',
                epsilon_order=str(r['B']+q*exact(kappa)),
                graph_floor='M^-q < c_i^(-original_rho_g*q) epsilon^(q*kappa)',
                extra_temporal_inverse_loss=0,
                treatment='bound the exact A1 remainder; never forward an additive tail to an inverse')


def materialize_flat(recipe,*,L,Bn,AJ,source,inverse_floor,rho_g,weight_factors,max_order=100):
    """Moderate exact diagnostic only. max_order is not an analytical gate."""
    order(max_order);J=recipe['raw_source_order'];q=recipe['integration_order']
    if J>max_order:raise ValueError('Keep the finite high-jet polynomial recipe symbolic')
    L,rho_g=map(exact,(L,rho_g))
    if L<=0 or rho_g<=0:raise ValueError('Original positive hull and graph exponent required')
    floor_power=inverse_floor*0+1
    for _ in range(ceil(rho_g*q)):floor_power=floor_power*(1+inverse_floor)
    return (recipe['scalar']*2**recipe['two_power']*3**recipe['three_power']*factorial(J)*L
            *Bn*AJ*source*weight_factors*floor_power)


class FastRowGain:
    """One-way join to the SAME old and sourced prefixes, no fine-row input."""
    def __init__(self,means,kappa=Q(1,100000)):
        self.means=means;self.before=means.before;self.generated=means.generated
        self.scalars=scalar.ScalarGainBudget(means);self.orders=rates(kappa);self.kappa=exact(kappa)
        p=Q(means.reweighting['original_pole']);L=Q(means.reweighting['L'])
        self.embedding=max(Q(1),L)**ceil(12*p)

    @lru_cache(None)
    def regular_rows(self,n):
        order(n);N=n+3;b=self.before
        f=difference(n,b.data('initial',N),b.fast(N),b.data('fixed',N))
        s=difference(n,b.postfast(N),b.slow_increment(N),b.data('fixed',N))
        D=self.embedding*P(self.generated.covariance(N).c)
        w=covariance_rows(n,D,self.means.data('fixed',N))
        # Q1=I-Haar has norm at most two. Preserve separate summands.
        out={key:fresh.monomial_upper(2*self.embedding*(f[key]+s[key]+w[key])) for key in ('Etheta','Ez')}
        return dict(rows=out,old_fast_remainder=f,old_slow_change=s,complete_wave_change=w,
                    input_order=N,regular_order=self.orders['D'])

    def postwave_recipe(self,n):
        a=self.regular_rows(n)
        flat=flat_recipe(n,self.kappa)
        p=Q(self.means.reweighting['original_pole']);w=Q(self.means.reweighting['width']);gap=Q(self.means.reweighting['gap'])
        flat['original_geometry']=dict(self.means.reweighting,source_pole=str(2*p),
            source_unweighting='Eedge(2p,width)',output_inverse_floor_power=ceil((w*w+4)/(gap*gap)))
        return dict(regular={k:v.record() for k,v in a['rows'].items()},
                    exact_A1=flat,order=str(self.orders['D']),
                    tail_ledger='original exact additive ledger unchanged and excluded from sourced inputs')

    def fast_recipe(self,n):
        """Same weighted fast map; all actual A1 velocities included."""
        order(n);m=self.means
        coeff=comb(n+5,5)*24*2**(n+4)*m.time_inverse(n)*factorial(n+4)
        rad=comb(n+5,5)*m.KT1(n+1)*24*2**(n+5)*m.time_inverse(n)*factorial(n+5)
        return dict(tangential_factor=coeff.record(),A1_velocity_factor=(1+m.KA1(n)).record(),
                    radial_factor=rad.record(),tangential_rows=self.postwave_recipe(n+4),
                    radial_rows=self.postwave_recipe(n+5),tangential_order=str(self.orders['D']),
                    radial_order=str(self.orders['D']+1))

    def inverse_source(self,*args,**kwargs):
        raise NotImplementedError('Finite A1 kernel/source polynomial recipes must be instantiated before using a materialized improved inverse callback')
