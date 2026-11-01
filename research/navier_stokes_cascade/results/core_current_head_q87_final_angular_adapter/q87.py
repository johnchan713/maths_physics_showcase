"""Final angular slot using same-state Haar and actual mean differences."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial
import importlib.util
_path=Path(__file__).resolve().parent.parent/'core_current_head_q86_final_mean_adapter/q86.py'
_spec=importlib.util.spec_from_file_location('q87_means',_path)
means=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(means)
full=means.full;convert=full.convert

def angular_difference(j,old,inc,fixed):
    full.order(j)
    x,u,z=(old[k] for k in ('x','u','z'));X,U,Z=(inc[k] for k in ('x','u','z'))
    b,V,G,g,r=(fixed[k] for k in ('b','V','G','graph','r'))
    g=1+2*g # slow derivative plus |v_r|_1<2 graph components
    d=j+4;dr=g*d;vis=d*d*(g*g+1)+r*dr+r*r
    rt=b*U+V*X+x*U+X*u+X*U;at=G*U+V*Z+z*U+Z*u+Z*U
    return d*U+(dr+2*r)*rt+d*at+vis*U

def native_factor(j,Rmax,Umin):
    """Normalized R input -> factorial-through U output, for 0<1/d<3."""
    full.order(j);H,lo=Q(Rmax),Q(Umin)
    if min(H,lo)<=0:raise ValueError('Positive same-hull coordinate enclosures required')
    jets=[Q(0)]+[H*comb(k+2,2)/lo**k+4*int(k==1) for k in range(1,j+1)]
    values=[Q(1)]
    for k in range(1,j+1):values.append(sum(jets[a]*values[k-a] for a in range(1,k+1)))
    return 243*factorial(j)*sum(values)

class FinalAngular:
    def __init__(self,final_means):
        if not isinstance(final_means,means.FinalMeans):raise TypeError('Q86 corrected first-real adapter required')
        self.means=final_means
    def native_coefficient(self,j,*,Umin):
        # Umin must enclose the native image of THIS same physical hull.
        H=self.means.current.current.profile.Rmax
        return native_factor(j,H,Umin)*self.coefficient(j)['Etheta']
    def coefficient(self,j):
        full.order(j);f=self.means
        # Recheck the invariant; do not silently accept a replaced state.
        check=means.FinalMeans(f.postslow)
        if check.state is not f.state:raise ValueError('Covariance/pressure state changed')
        if j+3>f.postslow.max_output:raise ValueError('Angular differences require factory reserve j+3')
        current=f.current;N=j+3;b=current.initial.budget
        fixed={k:convert(v) for k,v in b._data(N,b.pole).items()}
        old=current.old_full(N);fast=current.fast_unweighted(N);slow=current.slow(N,False)
        first=angular_difference(j,old,fast,fixed)
        after={k:old[k]+fast[k] for k in old}
        second=angular_difference(j,after,slow,fixed)
        haar=f.state.target.source.haar
        edge=full.refresh.maps.update.signed.maps.edge_absorption
        width=f.state.target.width
        initial=edge(2*haar.pole,width)*convert(haar.rows(j)['Etheta'])
        signed=means.covariance.defect(f.state,j+1,target_rate=Q(2))
        signed_unweighted=edge(signed['pole'],width)*signed['signed']
        H=initial+((j+1)*(fixed['graph']+1)+2*fixed['r'])*signed_unweighted
        return dict(Etheta=243*(H+first+second),Haar=H,fast=first,slow=second,
                    field_request=N,covariance_request=j+1,
                    angular_rate=2-Q(2,100000))

def record():
    k=Q(1,100000);C0=2-2*k
    return dict(status='same-state-final-angular-cancellation-adapter',
                rates=[str(v) for v in (C0,3-6*k,3-5*k)],
                Q84_angular_request=200000,factory_reserve=200003,profile_reserve=200030,
                physical_inequalities_certified=False,full_angular_rate_improved=False,
                unconditional_pulse_residual_or_blowup=False)
