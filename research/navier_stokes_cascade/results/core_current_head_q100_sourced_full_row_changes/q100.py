"""Full sourced-wave momentum row increments, not a stale incoming row."""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q99_changed_scalar_rows/q99.py'
s=importlib.util.spec_from_file_location('q100_scalars',p);scalars=importlib.util.module_from_spec(s);s.loader.exec_module(scalars)
pressure=scalars.pressure;EP=scalars.EP;cv=scalars.cv;full=scalars.full

def envelope(n,D,a,graph,pressure_next):
    full.order(n);a=Q(a)
    if a<=0 or len(D)<n+2:raise ValueError('Positive same hull and covariance through n+1 required')
    F=pressure.through(D,n);Fn=pressure.through(D,n+1)
    r=sum((a**(-j-1) for j in range(n+1)),Q(0))
    dr=(n+1)*(1+2*cv(graph))*Fn
    dz=(n+1)*Fn;dp=(n+1)*cv(pressure_next)
    theta=dr+2*r*F+dz;axial=dr+r*F+dz+dp
    return dict(radial_transport=dr,angular_connection=2*r*F,axial_connection=r*F,
                covariance_axial_gradient=dz,pressure_axial_gradient=dp,
                Etheta=theta,Ez=axial,Qtheta=2*theta,Qz=2*axial)

class SourcedRowChanges:
    def __init__(self,changed_scalars):
        if not isinstance(changed_scalars,scalars.ChangedScalars):raise TypeError('Q99 same-state changed scalar adapter required')
        self.scalars=changed_scalars;self.pressure=changed_scalars.pressure
        self.state=changed_scalars.state
    def bound(self,n,bundle):
        full.order(n)
        # Q98 now needs pressure through n+1, not only the integrated moment.
        pm=self.pressure.bound(n+1,bundle)
        b=self.scalars.means.current.initial.budget
        g=cv(b._data(n+1,b.pole)['graph'])
        cp=self.pressure.maps.profile
        out=envelope(n,bundle['change'],cp.Rmin,g,pm['pressure'])
        out.update(source_state=self.state,covariance_pole=bundle['covariance_pole'],
                   pressure_order_requested=n+1,covariance_order_requested=n+2,
                   increment_order='3/2-2*kappa',incoming_final_rows_retained=True,
                   total_changed_row_bound_supplied=False,means_updated=False,
                   actual_native_contracts_certified=False)
        return out
    def fresh_fast(self,*args,**kwargs):
        raise NotImplementedError('Incoming FINAL nonzero-Haar rows and their compact axial remainder must be composed before a fresh fast inverse')

def record():
    return dict(status='current-sourced-full-row-increment-envelope',
                increment_order='3/2-2*kappa',pressure_request='n+1',covariance_request='n+2',
                full_pressure_and_connections_retained=True,projection_norm_factor=2,
                incoming_final_rows_required=True,stale_prefast_row_reuse_rejected=True,
                total_changed_row_or_fresh_mean_supplied=False,
                actual_native_contracts_certified=False,complete_residual_or_blowup=False)
