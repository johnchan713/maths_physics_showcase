"""Full postfast correlations and same-profile fresh slow norm composition."""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q101_final_rows_fast_join/q101.py'
s=importlib.util.spec_from_file_location('q102_rows',p);rows=importlib.util.module_from_spec(s);s.loader.exec_module(rows)
full=rows.full;EP=rows.EP;cv=rows.cv
edge=full.refresh.maps.update.signed.maps.edge_absorption
KEYS=('P','Jtheta','Jz')

def correlation_join(n,before,old,increment,moment):
    full.order(n)
    if set(before)!=set(KEYS):raise ValueError('All current scalar channels required')
    delta=full.refresh.post.scalar.fast_change(n,old,increment,moment)
    return dict(scalars={k:cv(before[k])+cv(delta[k]) for k in KEYS},fast_change=delta)

def whole_axial(n,scalar_next,c_rho_next):
    full.order(n)
    return (n+1)*(cv(scalar_next['Jz'])+cv(c_rho_next)*cv(scalar_next['P']))

class PostfastSlow:
    def __init__(self,current_rows):
        if not isinstance(current_rows,rows.FinalRows):raise TypeError('Q101 derived full total-row factory required')
        self.rows=current_rows;self.state=current_rows.state
        self.means=current_rows.means
        self.changed=current_rows.changes.scalars
    def postfast(self,n,bundle):
        full.order(n)
        if self.means.state is not self.state or bundle['source_state'] is not self.state:
            raise ValueError('One identical incoming final state required')
        if bundle['stage']!='post_first_sourced_wave_fixed_means' or bundle['means_updated']:
            raise ValueError('Wave-only changed covariance before fresh means required')
        if n+6>self.rows.maximum or bundle['order']<n+8:
            raise ValueError('Postfast scalars need rows n+6 and covariance n+8')
        N=n+1;fast=self.rows.fast(N,bundle);p=self.rows.common_pole(bundle)
        if fast['order']<N or fast['pole']!=p or fast['source_state'] is not self.state:
            raise ValueError('Complete current fast field and its exact domain required')
        old=self.means.fields(N)
        inc={k:243*cv(edge(p,self.state.target.width))*cv(fast[k]) for k in ('x','u','z')}
        cp=self.rows.current.current.profile;a,b=cp.Rmin,cp.Rmax
        M=(b-a)*(1/a+1+b+b*b)
        changed=self.changed.bound(n,bundle)
        out=correlation_join(n,changed['changed'],old,inc,M)
        out.update(before_fast=changed['changed'],old_full=old,fast_unweighted=inc,
                   order=n,pole=p,source_state=self.state,scalar_order='3/2-kappa',
                   fast_change_order='5/2-3*kappa',means_updated=False,
                   actual_native_contracts_certified=False)
        return out
    def slow(self,n,bundle,weighted=False):
        full.order(n)
        if type(weighted) is not bool:raise TypeError('Choose the original weighted or unweighted slow norm')
        a=self.postfast(n,bundle)['scalars'];an=self.postfast(n+1,bundle)['scalars']
        b=self.rows.current.initial.budget
        k={key:cv(v) for key,v in b.slow(n).items()}
        kn={key:cv(v) for key,v in b.slow(n+1).items()}
        if set(k)!={'H','Kp','Kz','Kg','Kpsi'} or set(kn)!=set(k):
            raise ValueError('Same complete five-map coefficients required')
        fields=full.refresh.post.preserved.response(n,a,an,k,kn)
        if weighted:fields={key:cv(b.floor)*v for key,v in fields.items()}
        return dict(fields=fields,order=n,weighted=weighted,source_state=self.state,
                    tangential_order='3/2-kappa',radial_order='5/2-kappa',
                    actual_native_contracts_certified=False,means_updated=False)
    def axial(self,n,bundle):
        full.order(n)
        a=self.postfast(n+1,bundle)['scalars']
        cp=self.rows.current.current.profile
        return whole_axial(n,a,243*cv(cp.c_rho(n+1)))

def record():
    return dict(status='sourced-postfast-correlations-and-fresh-slow-map',
                scalar_order='3/2-kappa',fast_correlation_order='5/2-3*kappa',
                slow_tangential_order='3/2-kappa',slow_radial_order='5/2-kappa',
                postfast_axial_scalar_row_order='5/2-kappa',
                full_old_means_and_axial_fast_velocity_retained=True,
                wave_change_added_once=True,whole_moving_pressure_product_retained=True,
                actual_native_contracts_certified=False,postslow_reset_supplied=False,
                unconditional_pulse_residual_or_blowup=False)
