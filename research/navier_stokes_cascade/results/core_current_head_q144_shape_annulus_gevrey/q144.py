"""Mixed Gevrey-3 shape/annulus recurrence and the real shear floor."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
 s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
previous=load('q144_previous','core_current_head_q143_post_taper_gevrey/q143.py')
annular=load('q144_annular','core_current_head_q127_annular_radius_axis/q127.py')
old=load('q144_q70','core_current_head_q70_shape_annulus_mixed/q70.py')
mixed=previous.mixed;exact=mixed.exact
IDENTITY='same-Q143-endpoint-Q70-shape-and-Q127-earlier-annular-root'
REAL_ROOT=Q(572,10**16);WIDTH=Q(2,25)

def floor_certificate(B0,J6,root_cap=REAL_ROOT,width=WIDTH):
 """Use real first derivatives, not high-order seminorms, in a=1-2(log E)_y."""
 B0,J6,root_cap,width=map(exact,(B0,J6,root_cap,width))
 if B0<1 or not 0<=J6<=B0/100 or not 0<=root_cap<=REAL_ROOT or width!=WIDTH:
  raise ValueError('Same Q70 shape hierarchy, real accepted root and width required')
 T=2000*B0
 shape=Q(4,5)-128*(J6+6)/T       # ||sigma'||<=64
 q0=3*root_cap*64
 q1=3*root_cap*old.bump(1)       # three relative E bumps
 if q0>=1:raise ValueError('Original relative swirl nonvanishing required')
 earlier=Q(4,5)-2*q1/(1-q0)
 floor=min(shape,Q(4,5),earlier)
 if floor<=0:raise ValueError('Positive real shape/annular shear gate failed')
 return dict(stage='Q70-Q127-real-shear-floor',shape=shape,restoration=Q(4,5),
  earlier_annulus=earlier,floor=floor,T=T,relative_value=q0,relative_slope=q1,
  used_real_low_order_root_cap=True,high_order_complex_root_cap_used_for_floor=False)

def shape_annulus(post,H,Z,B0,J6,L0,Ecap,invEcap,R,P,length,identity):
 """Close fixed mixed envelopes from Xi through the earlier annular join."""
 vals=tuple(map(exact,(H,Z,B0,J6,L0,Ecap,invEcap,R,P,length)))
 H,Z,B0,J6,L0,Ecap,invEcap,R,P,length=vals
 if identity!=IDENTITY or post.get('stage')!='Q71-post-taper-value-based-mixed-Gevrey3':
  raise ValueError('Same Q143 endpoint and Q70/Q127 construction required')
 h,t,Xi,Lambda,old_length,Benv=post['reference_sources']['parameters']
 if not 0<H<=1 or Z<0 or B0<1 or not 0<=J6<=B0/100 or L0<J6 or min(Ecap,invEcap,P,length)<=0 or R<Xi or P<16 or Ecap<4*P or invEcap<20/P:
  raise ValueError('Original positive shape, strip, target and annular caps required')
 if length.denominator!=1 or length<2000*B0+10 or R/Xi>2**length.numerator:
  raise ValueError('Certified logarithmic shape/annular length required')
 endpoint=post['state']
 if Ecap<endpoint['E']['P'] or invEcap<endpoint['invE']['P']:
  raise ValueError('Shape value caps must contain the actual Q143 endpoint')
 root=annular.radius(H,Z,annular.PREMOD_ID)
 floor=floor_certificate(B0,J6)

 # Shape: differentiate the explicit logarithm at its actual value.  The
 # endpoint angular envelope and d_y log E driver enter the same recurrence.
 Be0=max(endpoint['E']['B'],endpoint['invE']['B'],endpoint['U']['B'],1/H)
 Etrace=mixed.env(endpoint['E']['P'],1,Be0,'constant-in-y-Gevrey3')
 Utrace=mixed.env(endpoint['U']['P'],1,Be0,'constant-in-y-Gevrey3')
 radial_step=mixed.env(162,max(Q(1),3260/(2000*B0)),1)
 ell=mixed.env(L0,1,Be0,'constant-in-y-Gevrey3')
 logf=mixed.env(6,1,max(Q(1),1/H),'constant-in-y-Gevrey3')
 dlog=mixed.add(mixed.env(Q(1,10)),mixed.mul(mixed.dy(radial_step),mixed.add(ell,logf)))
 zero=mixed.env(0,1,dlog['B'],'constant-in-y-Gevrey3')
 log_increment=mixed.integrate(dlog,zero,length)
 ratio_cap=max(Q(1),Ecap*endpoint['invE']['P'])
 Eshape=mixed.mul(Etrace,mixed.exp_comp(log_increment,ratio_cap))
 eta=mixed.env(8,1,max(Q(1),1/H))
 restore=mixed.env(162,3260,1)
 Ushape=mixed.add(Utrace,eta,mixed.mul(mixed.add(Utrace,eta),restore))

 # Earlier join: Q127 supplies a genuine complex radius for the same
 # eta-dependent root.  Q70 supplies the fixed real radial bumps.
 coeff=mixed.env(annular.D+annular.A,1,max(Q(1),1/root['sigma']))
 bump=mixed.env(162*3260,max(Q(1),8*3260/WIDTH),1)
 correction=mixed.scale(mixed.mul(coeff,bump),3)
 idealE=mixed.env(Ecap,1,max(Q(1),1/H))
 Eann=mixed.mul(idealE,mixed.add(mixed.env(1),correction))
 Uann=mixed.add(eta,mixed.scale(mixed.mul(coeff,bump),2))
 E=mixed.add(Eshape,Eann);U=mixed.add(Ushape,Uann)
 invE=mixed.reciprocal(E,1/invEcap)

 # Convert E=sqrt(2X)F and reintegrate all five ACTUAL Q143 rows.  The
 # crude rational sqrt caps  sqrt(2R)<=2R and 1/sqrt(2X)<=1 suffice.
 F=mixed.mul(E,mixed.env(1));invF=mixed.mul(invE,mixed.env(2*R))
 field=dict(stage='Q72-value-based-field-recurrence',F=F,invF=invF,U=U,
  exp_value_cap=Ecap,order_dependent_norm_exponentiated=False)
 T=max(v['P'] for v in endpoint['moments'].values());Bt=max(v['B'] for v in endpoint['moments'].values())
 traces=mixed.trace_rows(T,Bt);geo=mixed.geometry(R,1/Xi,2*R,1)
 state=mixed.state(field,traces,h,geo,length,mixed.IDENTITY)
 return dict(stage='conditional-Q70-Q127-shape-annulus-mixed-Gevrey3',field=field,state=state,
  reference_root=root,floor_certificate=floor,a_floor=floor['floor'],H=state['H'],Ay=state['Ay'],Be=state['Be'],
  shape_log_increment=log_increment,shape_multiplier_value_cap=ratio_cap,
  incoming_actual_moments=traces,full_p2_Pi_and_Pi_eta_retained=True,
  all_five_actual_incoming_rows_reintegrated=True,selected_order_zero_complex_and_loop_gates_certified=False)

def gate_request(recipe):
 if recipe.get('stage')!='conditional-Q70-Q127-shape-annulus-mixed-Gevrey3':
  raise ValueError('Same shape/annular recurrence required')
 s=recipe['state']
 return dict(stage='Q144-regional-gate-request',H=recipe['H'],Ay=recipe['Ay'],Be=recipe['Be'],
  afloor=recipe['a_floor'],p2cap=s['p']['P'],acap=s['a']['P'],bcap=s['b']['P'],
  missing=('same-profile complex value-cap binding','G>=max(2,abs(mu),abs(p2),abs(mu*p2))',
           'Eexp>=exp(G)','V0<=unchanged N'),complete=False)

def record():
 return dict(status='shape-annulus-mixed-Gevrey3-recurrence-and-positive-shear-floor',
  correction='Q70 real first-derivative bounds give a shape floor 1298/3125 and the accepted earlier root gives a larger annular floor',
  mixed='Q143 endpoint plus Q70 explicit shape/restoration and Q127 same-root radius give fixed y/eta Gevrey-3 envelopes',
  memory='all five actual Q143 moment rows are reintegrated; p2 retains Pi and Pi_eta',
  floor='a>=min(1298/3125,4/5,144709422735784/195312499997855)=1298/3125',
  blocking_inequality='bind the same-profile complex value caps, G>=max(2,|mu|,|p2|,|mu*p2|), Eexp>=exp(G), and V0<=N before Q140-Q138',
  actual_complete_profile_pulse_covariance_residual_or_blowup_certified=False,historical_frozen_records_modified=False)
