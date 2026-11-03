"""Scalar-input analytic disks and real Gevrey binding of the original loop."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('q130_previous',HERE.parent/'core_current_head_q129_moving_inverse_gevrey/q129.py')
previous=importlib.util.module_from_spec(s);s.loader.exec_module(previous)
cutoff=previous.previous;exact=previous.exact
IDENTITY='same-six-scalar-loop-inputs-with-full-pressure-and-augmented-variance'
def scalar_jets(state,recipe,Ecap,afloor):
 if state.get('stage')!='PREMOD-state' or recipe.get('stage')!='real-Gevrey3-cutoff-root':raise ValueError('Q128 complete PREMOD state and same cutoff/root recipe required')
 A,B,P,rho=map(exact,(state['a'],state['b'],state['p'],state['rho']))
 e,afloor=exact(Ecap),exact(afloor)
 if min(A,e,afloor,rho)<=0 or min(B,P)<0:raise ValueError('Actual consistent scalar caps and positive original floor required')
 ra=recipe['a_radius']
 if ra!=min(rho/2,rho*afloor/(4*A)):raise ValueError('Same original a floor/radius required')
 L=max(Q(1),1/rho,1/ra,recipe['amplitude_base'],recipe['mu_base'])
 Ms=recipe['amplitude_prefactor'];Mmu=max(recipe['real_mu_value_cap'],recipe['mu_prefactor'])
 Mt=2*B/afloor;Mv=A+2*B*B/afloor+A*Ms*Ms
 N=A+Mt+P+Mmu+e+Mv
 return dict(stage='six-real-scalar-Gevrey3',prefactors=dict(a=A,ts=Mt,p=P,mu=Mmu,E=e,v=Mv),
  N=N,L=3*L,actual_PREMOD_inputs_certified=False,slow_analytic_extension_claimed=False)
def parameter_disk(G,Eexp,identity):
 G,Eexp=exact(G),exact(Eexp)
 if identity!=IDENTITY or G<2 or Eexp<1:raise ValueError('Same original six real caps and outward exp(G) bound required')
 E2=9*Eexp**2
 r0=min(Q(1,4),1/(8*G),1/(64*G*Eexp))
 T=2*G+24*G*Eexp*E2
 r=min(r0/2,r0/(112*T*(G+1)))
 W=4*G*G*(1+T*T);F=8*G**3*(T+1)
 return dict(stage='analytic-scalar-parameter-disk',G=G,Eexp=Eexp,E2=E2,r0=r0,T=T,r=r,W=W,F=F,
  real_caps_and_outward_exponential_bound_are_premises=True,actual_profile_certified=False)
def bind(jets,disk):
 if jets.get('stage')!='six-real-scalar-Gevrey3' or disk.get('stage')!='analytic-scalar-parameter-disk':raise ValueError('Same scalar jets and analytic parameter disk required')
 density_floor=1/(8*disk['G']**2)
 B=jets['L']*(1+jets['N']/disk['r']);C=1/disk['r']
 inv=previous.inverse(disk['W'],B,C,density_floor,previous.IDENTITY)
 primitive=previous.compose(inv,disk['F'],B,C,previous.IDENTITY)
 return dict(stage='bound-original-loop-conditional',W=disk['W'],F=disk['F'],B=B,C=C,
  inverse=inv,primitive=primitive,density_floor=density_floor,actual_original_state_caps_floors_and_real_G_caps_certified=False,
  actual_final_swirl_repair_cone_or_scale_certified=False,slow_analytic_extension_claimed=False)
def record():
 return dict(status='six-scalar-real-Gevrey-to-original-mixed-loop-binding',
  scalar_parameters='a, ts=-b/a, p=p2, mu, E, v=a+b^2/a+a*s^2',
  density_floor='1/(8*G^2) from a>=1/G, 0<v<=G and 2*pi<8; not G^-3 from G>=2 alone',
  mixed='|d_eta^i d_theta^j (w,fA,fB)|/(i!j!) <= (W,F,F)*B^i*C^j*(i!)^2',
  pressure='Same Q128 p2 with complete incoming Pi and Pi_eta',
  exponential='Eexp>=exp(G) is a premise; E2=9*Eexp^2>=exp(2*(G+1))',
  scope='Removable signed loop, fixed phase origin and original augmented-variance means; Q129 inverse/primitive bounds apply',
  actual_state_caps_swirl_final_repair_scale_pulse_residual_or_blowup_certified=False,
  slow_analytic_extension_claimed=False)
