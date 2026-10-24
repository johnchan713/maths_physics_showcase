"""First-real fresh-mean scalar protocol and deferred, full-request obligations."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 sp=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
native=load('real_refresh_native','core_native_haar_primitive_join/bounds.py')
post=load('real_refresh_post','core_post_signed_mean_closure/budgets.py')
EP,order,exact=native.EP,native.order,native.exact

def rates(kappa=Q(1,100000)):
 k=exact(kappa)
 if not 0<k<=Q(1,100):raise ValueError('Retained positive small-kappa regime required')
 B=1-k;alpha=2-3*k;D=alpha-k;nu=3-6*k
 p1=post.strict_integrations(D,nu,k);p0=post.strict_integrations(B,nu,k)
 groups=dict(post_signed_Haar=nu,fast_difference=D+1-2*k,slow_difference=alpha+1-2*k,A1=D+p1*k)
 return dict(B=B,alpha=alpha,D=D,nu=nu,fast=D,slow=alpha,
  fast_scalar_change=min(B+D,2*D),postfast_scalar=min(Q(2),alpha,B+D,2*D),
  postslow_scalar=min(alpha+1,2*alpha,alpha+2),whole_postslow_axial=alpha+2,
  tangential_groups=groups,conditional_mean_order=min(groups.values()),
  radial_compensator=alpha+1,A0=B+p0*k,A1_integrations=p1,A0_integrations=p0,
  strict_fast_difference_margin=D+1-2*k-nu,
  complete_wave_mean_feedback_order=min(D,alpha),
  actual_mean_block_or_full_residual_certified=False)

def requests(n,*,old_pole,target_allowance,retained_wave_pole,kappa=Q(1,100000)):
 order(n);p,d,qold=map(exact,(old_pole,target_allowance,retained_wave_pole))
 if min(p,d,qold)<0:raise ValueError('Nonnegative whole-request pole bounds required')
 r=rates(kappa);trace=post.dependency_trace(n,r['A1_integrations'],r['A0_integrations'])
 N=trace['maxima']['new_covariance'];q=max(qold,6*p,6*N+9,6*N+d+6)
 return dict(trace=trace,q=q,operator_domains=dict(KT1=2*q,KA1=2*q,Kp=(2*q,4*q),Kmom=(2*q,4*q)),
  target_allowance=d,old_pole=p,actual_input_bounds_through_entire_nested_request_graph_required=True,
  low_order_parent_inputs_sufficient=False,actual_high_jets_evaluated=False)

def current_postfast(n,pressure_result,old_full,fast_increment,moment):
 """pressure_result already includes OLD scalar + COMPLETE signed change.

 Add only actual fast_change; adding wave_change again double counts it.
 Mean inputs through n+1, including the actual axial gd-A1 gd field.
 """
 order(n)
 if set(pressure_result['updated_scalars'])!={'P','Jtheta','Jz'}:raise ValueError('All current post-signed scalar defects required')
 f=post.scalar.fast_change(n,old_full,fast_increment,moment)
 return {k:EP(pressure_result['updated_scalars'][k])+f[k] for k in ('P','Jtheta','Jz')}

def slow_response(n,current,current_next,fixed,fixed_next):
 return post.preserved.response(n,current,current_next,fixed,fixed_next)

def final_tangential(n,haar,A1,old_full,fast_increment,slow_increment,fixed):
 """Conditional inherited full differences: all fields through n+3."""
 return post.final_row_envelope(n,haar,A1,old_full,fast_increment,slow_increment,fixed)
