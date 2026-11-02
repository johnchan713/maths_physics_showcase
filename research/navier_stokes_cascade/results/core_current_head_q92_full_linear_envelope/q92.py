"""Complete regular linear envelope; lower weighted norms remain hypotheses."""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q91_weighted_mean_wave/q91.py'
s=importlib.util.spec_from_file_location('q92_meanwave',p);meanwave=importlib.util.module_from_spec(s);s.loader.exec_module(meanwave)
full=meanwave.full;EP=full.EP;convert=full.convert
KEYS={'b','G','g','r','normal','Eik','bR','Bz','c','K','T','R','Pi','psi'}

def envelope(n,data,harmonic=1):
    full.order(n)
    if set(data)!=KEYS:raise ValueError('All regular-wave, pressure and full-background slots required')
    if type(harmonic) is not int or harmonic==0:raise ValueError('Frozen nonzero integer harmonic required')
    a={k:convert(v) for k,v in data.items()};h=abs(harmonic);d=n+2
    b,G,g,r,N,e,bR,Bz,c,K,T,R,pi,psi=(a[k] for k in ('b','G','g','r','normal','Eik','bR','Bz','c','K','T','R','Pi','psi'))
    groups=dict(slow=d*(1+b*g+G),phase=2*h*e,base=bR+Bz+b*r,
        viscosity=2*d*d*g*g+d*d+r*d*g+r*r,
        mixed=4*h*N*d*(g+1),divergence=2*h*N*(d*(g+1)+r),angular=4*h*N*r)
    E=sum(groups.values(),EP());principal=2*c*d+K+4*h*h*N*N
    return dict(F=psi*(E*T+d*(g+1)*pi+(principal+E)*R),groups=groups,
                pressure=psi*d*(g+1)*pi,principal_on_curl=psi*principal*R,
                input_order=n+2,root_or_target_order=n+3,actual_inputs_certified=False)

def gains(k=Q(1,100000)):
    k=Q(k)
    if not 0<k<Q(1,4):raise ValueError('Retained kappa range required')
    return dict(slow=1-k,phase=Q(1,2),base=Q(1),viscosity=1-2*k,
                mixed=Q(1,2)-k,divergence=Q(1,2)-k,angular=Q(1,2))

def record():
    return dict(status='explicit-complete-regular-linear-source-envelope',
        input_order_shift=2,complete_curl_root_order_shift=3,second_graph_multiplier_differentiated=True,
        pressure_and_all_viscous_groups_retained=True,regular_and_cut_wave_slots_identified=False,
        independent_Clin_callback_removed_conditionally=True,weighted_wave_inputs_certified=False,
        full_actual_source_residual_or_blowup=False)
