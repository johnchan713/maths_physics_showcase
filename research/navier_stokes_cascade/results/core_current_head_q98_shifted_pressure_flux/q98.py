"""Conditional current covariance pressure via the whole-hull flux identity."""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q97_sourced_covariance_moments/q97.py'
s=importlib.util.spec_from_file_location('q98_covariance',p);covariance=importlib.util.module_from_spec(s);s.loader.exec_module(covariance)
full=covariance.full;EP=covariance.EP;cv=covariance.cv

def through(D,n):
    full.order(n)
    if len(D)<=n:raise ValueError('All raw mixed covariance jets required')
    return sum((Q(5**j,factorial(j))*cv(D[j]) for j in range(n+1)),EP())

def envelope(n,D,a,b,edge,rho,KT,KA):
    """Five-variable factorial-through norms, at the SAME covariance pole.

    Conditional on weighted derivative and smooth whole-hull trace contracts.
    rho is the weighted normalized through-n norm of the original compensator.
    """
    full.order(n);a,b=Q(a),Q(b)
    if not 0<a<b or len(D)<n+2:raise ValueError('Positive full hull and covariance through n+1 required')
    F=through(D,n);Z=(n+1)*through(D,n+1)
    radius=sum((a**(-j-1) for j in range(n+1)),Q(0))
    H=2*radius*F+Z
    P=(b-a)*cv(edge)*H
    compensated=H+cv(rho)*P
    return dict(local_flux=F,H=H,P=P,compensated=compensated,
                pressure=F+cv(KT)*compensated,
                A0=cv(KA)*compensated,
                radial_residual=cv(rho)*P+cv(KA)*compensated)

class ShiftedPressure:
    def __init__(self,changed,*,c0,c1,shape,time_inverse,maximum,pole):
        if not isinstance(changed,covariance.SourcedCovariance):raise TypeError('Current Q97 covariance factory required')
        self.changed=changed;self.state=changed.state
        scalars=changed.response.inputs.source.means.current
        self.maps=scalars.current
        old=scalars.fast.rows.pressure
        self.operator=full.compact(self.maps,old.d,c0,c1,shape,time_inverse,maximum,pole,0)
    def bound(self,n,bundle):
        full.order(n)
        if bundle['source_state'] is not self.state or bundle['stage']!='post_first_sourced_wave_fixed_means' or bundle['means_updated']:
            raise ValueError('Identical fixed-mean changed-wave stage required')
        if bundle['order']<n+1:raise ValueError('Current complete covariance through n+1 required')
        p=bundle['covariance_pole'];op=self.operator(n,p)
        cp=self.maps.profile
        edge=full.refresh.maps.update.signed.maps.edge_absorption(p,self.state.target.width)
        out=envelope(n,bundle['change'],cp.Rmin,cp.Rmax,edge,cv(cp.weighted_rho(n)),op['KT'],op['KA'])
        out.update(output_order=n,covariance_pole=p,source_state=self.state,
                   epsilon_order='3/2-kappa',whole_hull_trace_contract_required=True,
                   actual_native_contracts_certified=False,means_updated=False)
        return out

def record():
    return dict(status='current-shifted-pressure-flux-bound',pressure_increment_order='3/2-kappa',
        direct_source_order='3/2-2*kappa',covariance_request='n+1',
        same_moving_T0_and_A0=True,local_pressure_flux_retained=True,
        pressure_memory_added_not_reset=True,whole_hull_mixed_trace_contract_required=True,
        actual_native_contracts_certified=False,fresh_means_or_complete_residual_supplied=False,
        complete_residual_or_blowup=False)
