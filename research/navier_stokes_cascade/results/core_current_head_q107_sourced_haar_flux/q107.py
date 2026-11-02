"""Current fixed-mean sourced Haar flux; full old rows remain additive."""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q100_sourced_full_row_changes/q100.py'
s=importlib.util.spec_from_file_location('q107_changes',p);changes=importlib.util.module_from_spec(s);s.loader.exec_module(changes)
pressure=changes.pressure;full=changes.full;EP=changes.EP;cv=changes.cv

def envelope(n,D,a,pressure_next):
    full.order(n);a=Q(a)
    if a<=0 or len(D)<n+2:raise ValueError('Positive common hull and covariance through n+1 required')
    F=pressure.through(D,n);Fn=pressure.through(D,n+1)
    r=sum((a**(-j-1) for j in range(n+1)),Q(0));d=n+1
    # Explicit epsilon in slow gradients is dropped only in coefficient arithmetic.
    return dict(theta=d*Fn+2*r*F+d*Fn,
                axial=d*Fn+r*F+d*(Fn+cv(pressure_next)),
                radial_flux=d*Fn+2*r*F,slow_covariance=d*Fn,
                slow_pressure=d*cv(pressure_next))

class HaarChanges:
    def __init__(self,sourced_changes):
        if not isinstance(sourced_changes,changes.SourcedRowChanges):raise TypeError('Q100 current sourced changes required')
        self.changes=sourced_changes;self.pressure=sourced_changes.pressure;self.state=sourced_changes.state
        if not isinstance(self.pressure,pressure.ShiftedPressure) or self.pressure.state is not self.state:
            raise ValueError('Same Q98 complete shifted pressure required')
    def bound(self,n,bundle):
        full.order(n)
        if self.pressure.state is not self.state or self.changes.state is not self.state:raise ValueError('Current sourced pressure state changed')
        pm=self.pressure.bound(n+1,bundle) # original stage, pole, trace and order guards
        out=envelope(n,bundle['change'],self.pressure.maps.profile.Rmin,pm['pressure'])
        out.update(source_state=self.state,covariance_pole=bundle['covariance_pole'],order=n,
                   covariance_request=n+2,pressure_request=n+1,
                   epsilon_order='3/2-kappa',graph_multiplier_required=False,
                   incoming_haar_rows_still_required=True,means_updated=False,
                   actual_native_contracts_certified=False,full_cone_or_residual_certified=False)
        return out

def rates(kappa=Q(1,100000)):
    k=Q(kappa)
    if not 0<k<Q(1,4):raise ValueError('Original 0<kappa<1/4 required')
    return dict(changed_haar=Q(3,2)-k,full_changed_row=Q(3,2)-2*k,
                gain=k,old_conditional_mean_target=2-2*k,
                remaining_target_gap=Q(1,2)-k)

def record():
    return dict(status='current-sourced-haar-flux-and-compensated-pressure',
                changed_haar_order='3/2-kappa',full_changed_row_order_unchanged='3/2-2*kappa',
                local_pressure_and_moving_rho_retained=True,
                haar_cutoff_term_cancels_only_after_unit_mass_and_trace_identity=True,
                incoming_haar_or_target_cancellation_certified=False,
                actual_native_contracts_certified=False,full_cone_pulse_residual_or_blowup=False)
