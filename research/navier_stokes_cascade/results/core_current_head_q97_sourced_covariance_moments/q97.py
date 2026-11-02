"""Changed complex covariance and fixed-mean wave moment increments."""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q96_sourced_response/q96.py'
s=importlib.util.spec_from_file_location('q97_response',p);response=importlib.util.module_from_spec(s);s.loader.exec_module(response)
full=response.full;EP=response.EP;cv=response.cv;prod=response.prod

def covariance_envelope(old,new):
    cross=prod(old,new);square=prod(new,new)
    return dict(cross=[2*x for x in cross],square=square,
                change=[2*x+y for x,y in zip(cross,square)])

def moment_envelope(D,a,b,edge,c_rho):
    a,b=Q(a),Q(b)
    if not 0<a<b or len(D)<3 or len(c_rho)!=len(D)-1:raise ValueError('Positive full hull and two extra covariance jets required')
    D=list(map(cv,D));c=list(map(cv,c_rho));L=(b-a)*cv(edge)
    P=[2*Q(1,a)*L*D[j]+L*D[j+1] for j in range(len(D)-1)]
    Jt=[L*b*b*D[j] for j in range(len(D)-1)]
    Jz=[2*L*b*D[j]+Q(1,2)*L*b*b*D[j+1] for j in range(len(D)-1)]
    pressure_core=[L*b*D[j]+Q(1,2)*L*b*b*D[j+1] for j in range(len(D)-1)]
    cmP=prod(c,P);pressure_moment=[x+y for x,y in zip(pressure_core,cmP)]
    axial=[Jz[j+1]+cmP[j+1] for j in range(len(D)-2)]
    return dict(P=P,Jtheta=Jt,Jz=Jz,pressure_moment=pressure_moment,whole_axial=axial)

class SourcedCovariance:
    def __init__(self,sourced_response):
        if not isinstance(sourced_response,response.SourcedResponse):raise TypeError('Same-state Q96 complete response required')
        self.response=sourced_response;self.state=sourced_response.state
    def covariance(self,N,*args,**kwargs):
        full.order(N);old=self.response.inputs.source.waves
        if old.state is not self.state:raise ValueError('One identical incoming first-real state required')
        if 'harmonic' in kwargs:raise ValueError('Both original +/-1 and +/-2 response pairs are required')
        new=[self.response.build(N,*args,harmonic=m,**kwargs) for m in (1,2)]
        A=[factorial(j)*old.through(j)['U'] for j in range(N+1)]
        V=[2*sum((v['complete'][j] for v in new),EP()) for j in range(N+1)]
        poles=[old.through(j)['pole'] for j in range(N+1)]+[v['pole'] for v in new]
        q=max(poles);out=covariance_envelope(A,V)
        return dict(out,old=A,new=V,order=N,wave_pole=q,covariance_pole=2*q,
                    mean_factory_reserve=N+3,source_state=self.state,
                    stage='post_first_sourced_wave_fixed_means',
                    actual_native_contracts_certified=False,means_updated=False)
    def moments(self,n,bundle):
        full.order(n)
        if bundle['source_state'] is not self.state or bundle['order']<n+2:
            raise ValueError('Same changed covariance through n+2 required for the whole axial product')
        cp=self.state.target.source.profile
        edge=full.refresh.maps.update.signed.maps.edge_absorption(bundle['covariance_pole'],self.state.target.width)
        c=[243*factorial(j)*cv(cp.c_rho(j)) for j in range(n+2)]
        out=moment_envelope(bundle['change'][:n+3],cp.Rmin,cp.Rmax,edge,c)
        out.update(covariance_pole=bundle['covariance_pole'],output_order=n,
                   full_pointwise_pressure_bound_supplied=False,means_updated=False)
        return out

def record():
    return dict(status='first-sourced-complex-covariance-and-global-moment-increments',
        covariance_order='3/2-kappa',new_square_order='2-2*kappa',
        pointwise_radial_source_order='3/2-2*kappa',Haar_moment_order='3/2-kappa',
        whole_axial_change_order='5/2-kappa',both_complex_cross_orientations=True,
        old_covariance_and_pressure_memory_retained=True,mean_velocities_held_fixed=True,
        full_pointwise_pressure_or_fresh_means_supplied=False,
        actual_native_contracts_certified=False,complete_residual_or_blowup=False)
