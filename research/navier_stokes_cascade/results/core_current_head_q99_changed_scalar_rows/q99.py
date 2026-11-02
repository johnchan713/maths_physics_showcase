"""Postslow plus sourced-wave scalar moments; no fresh-mean assertion."""
from pathlib import Path
from math import factorial
from fractions import Fraction as Q
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q98_shifted_pressure_flux/q98.py'
s=importlib.util.spec_from_file_location('q99_pressure',p);pressure=importlib.util.module_from_spec(s);s.loader.exec_module(pressure)
EP=pressure.EP;cv=pressure.cv;full=pressure.full
KEYS=('P','Jtheta','Jz')

def join(n,old,old_axial,moments):
    full.order(n)
    if set(old)!=set(KEYS):raise ValueError('Complete incoming POSTSlow scalars required')
    if any(len(moments[k])<n+2 for k in KEYS) or len(moments['whole_axial'])<n+1:
        raise ValueError('Moment derivatives through n+1 and whole axial through n required')
    delta={k:pressure.through(moments[k],n+1) for k in KEYS}
    changed={k:cv(old[k])+delta[k] for k in KEYS}
    # Bound integrated row changes after factoring their one explicit epsilon.
    theta=sum((Q(5**j,factorial(j))*moments['Jtheta'][j+1] for j in range(n+1)),EP())
    axial=pressure.through(moments['whole_axial'],n)
    return dict(old={k:cv(old[k]) for k in KEYS},delta=delta,changed=changed,
                integrated_theta=(n+1)*cv(old['Jtheta'])+theta,
                integrated_axial=cv(old_axial)+axial,
                delta_integrated_theta=theta,delta_integrated_axial=axial)

class ChangedScalars:
    def __init__(self,shifted_pressure):
        if not isinstance(shifted_pressure,pressure.ShiftedPressure):raise TypeError('Q98 same-state full shifted pressure required')
        self.pressure=shifted_pressure;self.changed=shifted_pressure.changed
        self.means=self.changed.response.inputs.source.means
        self.state=self.changed.state
    def bound(self,n,bundle):
        full.order(n)
        # Q97 checks n+2 covariance reserve; Q98 also enforces exact stage/pole.
        moment=self.changed.moments(n,bundle)
        pm=self.pressure.bound(n,bundle)
        post=self.means.postslow
        old=post.scalars(n+1)
        out=join(n,old,post.axial(n),moment)
        out.update(pressure_increment=pm,source_state=self.state,
                   scalar_order='3/2-kappa',integrated_row_order='5/2-kappa',
                   covariance_pole=bundle['covariance_pole'],scalar_output_order=n+1,
                   fresh_mean_source_rows_supplied=False,means_updated=False,
                   actual_native_contracts_certified=False)
        return out

def record():
    return dict(status='postslow-sourced-scalar-and-integrated-row-join',
                old_scalar_stage='completed-first-real-postslow',wave_change_added_once=True,
                scalar_order='3/2-kappa',integrated_row_order='5/2-kappa',
                whole_pressure_compensator_derivative_retained=True,
                full_nonzero_Haar_rows_or_fresh_means_supplied=False,
                actual_native_contracts_certified=False,complete_residual_or_blowup=False)
