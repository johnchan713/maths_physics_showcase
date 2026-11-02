"""Current incoming Haar means and the sourced postslow scalar reset."""
from pathlib import Path
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q102_sourced_postfast_slow/q102.py'
s=importlib.util.spec_from_file_location('q103_slow',p);slow=importlib.util.module_from_spec(s);s.loader.exec_module(slow)
full=slow.full;EP=slow.EP;cv=slow.cv

def haar_sum(initial,previous):
    if set(initial)!={'X','U','Z'} or set(previous)!={'x','u','z'}:
        raise ValueError('Both complete slow histories required')
    return {k:243*(cv(initial[j])+cv(previous[k])) for k,j in (('x','X'),('u','U'),('z','Z'))}

def reset(n,old_haar,increment,base_b,moment):
    full.order(n)
    return full.refresh.post.preserved.haar_slow_reset(n,old_haar,increment,base_b,moment)

class CurrentPostslow:
    def __init__(self,fresh):
        if not isinstance(fresh,slow.PostfastSlow):raise TypeError('Q102 same-state fresh slow norm map required')
        self.fresh=fresh;self.state=fresh.state;self.current=fresh.rows.current
        b=self.current.initial.budget;cp=self.current.current.profile
        a,z=b.hull['Rminus'],b.hull['Rplus']
        if not cp.Rmin<=a<z<=cp.Rmax:raise ValueError('Whole reserved potential hull must fit current norm hull')
        self.moment=(z-a)*(1/a+1+z+z*z)
    def check(self,n,bundle):
        full.order(n)
        if self.fresh.means.state is not self.state or bundle['source_state'] is not self.state:
            raise ValueError('Incoming Haar and wave state must agree')
        if bundle['stage']!='post_first_sourced_wave_fixed_means' or bundle['means_updated']:
            raise ValueError('Original wave-only bundle required for the derived fresh means')
        if n+9>self.fresh.rows.maximum or bundle['order']<n+11:
            raise ValueError('Postslow reset needs rows n+9 and covariance n+11')
    def old_haar(self,n):
        full.order(n)
        if n>self.fresh.means.postslow.max_output:raise ValueError('Outside old complete slow-history reserve')
        initial=self.current.initial.slow_fields(n,False)
        previous=self.current.slow(n,False)
        return haar_sum(initial,previous)
    def scalars(self,n,bundle):
        self.check(n,bundle);N=n+2
        old=self.old_haar(N);inc=self.fresh.slow(N,bundle,False)
        if inc['source_state'] is not self.state or inc['order']<N or inc['weighted'] is not False:
            raise ValueError('Same-state unweighted complete fresh slow field required')
        b=self.current.initial.budget
        base=243*cv(b.fixed(N,2*b.pole)['b'])
        out=reset(n,old,inc['fields'],base,self.moment)
        return dict(scalars={k:out[k] for k in ('P','Jtheta','Jz')},old_haar=old,
                    fresh_slow=inc['fields'],order=n,scalar_order='5/2-kappa',
                    source_state=self.state,norm_stage='after_current_fresh_slow',
                    actual_native_contracts_certified=False,means_updated=False)
    def axial(self,n,bundle):
        full.order(n)
        a=self.scalars(n+1,bundle)['scalars'];cp=self.current.current.profile
        return slow.whole_axial(n,a,243*cv(cp.c_rho(n+1)))
    def compensator(self,n,bundle):
        full.order(n)
        a=self.scalars(n,bundle)['scalars'];cp=self.current.current.profile
        return 243*cv(cp.weighted_rho(n))*a['P']

def record():
    return dict(status='current-haar-postslow-reset-and-pressure-scalars',
                incoming_haar='initial_slow_plus_previous_first_real_slow',
                full_fast_correlations_retained_in_Q102=True,
                scalar_and_compensator_order='5/2-kappa',whole_axial_order='7/2-kappa',
                moving_pressure_product_retained=True,original_whole_hull_G0_required=True,
                actual_native_contracts_certified=False,full_pointwise_pressure_or_radial_residual_supplied=False,
                unconditional_pulse_residual_or_blowup=False)
