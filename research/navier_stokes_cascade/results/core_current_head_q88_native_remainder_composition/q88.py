"""Native first-real nonlinear remainder composition over positive EP bounds."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial
import importlib.util
_path=Path(__file__).resolve().parent.parent/'core_current_head_q87_final_angular_adapter/q87.py'
_spec=importlib.util.spec_from_file_location('q88_angular',_path)
angular=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(angular)
means=angular.means;full=angular.full;EP=full.EP;convert=full.convert

def geometry(N,r,H):
    full.order(N);r,H=Q(r),Q(H)
    if not 0<r<=H:raise ValueError('Ordered same-profile positive radial hull required')
    lo=min(Q(1),r**3);hi=max(Q(1),H**3);q=max(Q(1),1/lo);B=max(Q(1),H**2,1/r)
    return dict(Umin=lo,Umax=hi,r=(N+1)*comb(N+2,2)*q**N/r,g=3*B*(N+1)**2*q**N)

def polynomial(m,a):
    full.order(m)
    keys=('b0','V','G','x','u','z','W','Z','e','r','g','c','T')
    if set(a)!=set(keys):raise ValueError('All Q84 norm slots required')
    b,V,G,x,u,z,W,Z,e,r,g,c,T=(convert(a[k]) for k in keys)
    d=m+7;D=3*d*(1+g);av=2*V*r;aq=6*V*r*r+2*d*g*V*r
    J=b*u+V*x+x*u+Z;L=G*u+V*z+z*u+W
    F=e+d*u+d*L+(D*D+r*D+r*r+d*d)*u
    H=b*z+G*x+x*z+W
    rest=d*(1+4*c)*x+r*(2*b*x+x*x)+d*H+2*r*Z+d*d*x
    return 2*(aq*T*J+av*T*F+r*u*u+rest)

class NativeRemainder:
    def __init__(self,postslow):
        self.means=means.FinalMeans(postslow);self.angular=angular.FinalAngular(self.means)
    def inputs(self,m):
        full.order(m);f=self.means
        if means.FinalMeans(f.postslow).state is not f.state:raise ValueError('First-real pressure state changed')
        if m+7>f.postslow.max_output:raise ValueError('Common factory reserve m+7 required')
        N=m+6;cp=f.current.current.profile
        geom=geometry(N,cp.Rmin,cp.Rmax);factor=angular.native_factor(N,cp.Rmax,geom['Umin'])
        fields=f.fields(N)
        wave=means.covariance.defect(f.state,N,target_rate=Q(2))
        edge=full.refresh.maps.update.signed.maps.edge_absorption(wave['pole'],f.state.target.width)
        budget=f.current.initial.budget;base=budget.fixed(N,budget.pole)
        compact=budget.compact(0,budget.pole)
        if set(compact)!={'KT','KA','I'}:raise ValueError('Original compact/time floor contract required')
        a={k:factor*fields[k] for k in ('x','u','z')}
        a.update({k:factor*edge*wave[k] for k in ('W','Z')})
        a.update({k:factor*convert(base[source]) for k,source in (('b0','b'),('V','V'),('G','G'),('c','c'))})
        a.update(e=self.angular.native_coefficient(m+4,Umin=geom['Umin']),
                 r=EP(geom['r']),g=EP(geom['g']),
                 T=9600*6**m*factorial(m+4)*convert(compact['I'])*EP({1:1}))
        return dict(data=a,chart='factorial-through-native-U',source_order=N,angular_order=m+4,
                    factory_reserve=m+7,Umin=geom['Umin'],Umax=geom['Umax'])
    def coefficient(self,m):return polynomial(m,self.inputs(m)['data'])

def record():
    return dict(status='same-hull-native-nonlinear-radial-composition',
                radial_shell_source=199996,factory_reserve=200003,profile_reserve=200030,
                temporal_raw_constant=9600,temporal_extra_jets=4,
                independent_final_field_or_temporal_callback=False,
                actual_continuum_inequalities_certified=False,unconditional_residual_or_blowup=False)
