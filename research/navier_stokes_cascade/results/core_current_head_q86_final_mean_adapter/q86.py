"""Final first-real mean coefficient adapter; no outgoing fine-row callback."""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util
_path=Path(__file__).resolve().parent.parent/'core_current_head_q85_covariance_defect_adapter/q85.py'
_spec=importlib.util.spec_from_file_location('q86_covariance',_path)
covariance=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(covariance)
full=covariance.full

def ledger(N):
    full.order(N)
    return dict(final_fields=N,incoming_rows=N+10,covariance_state=N+11,
                signed_target=N+14,profile_covariance=N+27,
                full_nested_DAG_expanded=False)

class FinalMeans:
    def __init__(self,postslow):
        if not isinstance(postslow,full.post.Postslow):raise TypeError('Corrected Q47 Postslow object required')
        current=postslow.current;rows=current.fast.rows
        if not isinstance(rows,full.FullRows):raise TypeError('Full shifted-pressure Q47 rows required, not the old Haar norm')
        if rows.current is not current.current or rows.initial is not current.initial:
            raise ValueError('One identical current pressure and initial mean state required')
        self.postslow=postslow;self.current=current;self.state=current.current.state
        if self.state.target.root is not self.state.primary.primary.root:raise ValueError('Original root identity lost')
    def fields(self,N):
        full.order(N)
        if N>self.postslow.max_output:raise ValueError('Outside full factory reserve')
        old=self.current.old_full(N);fast=self.current.fast_unweighted(N);slow=self.current.slow(N,False)
        if any(set(a)!={'x','u','z'} for a in (old,fast,slow)):raise ValueError('All three complete velocity components required')
        return {k:243*(old[k]+fast[k]+slow[k]) for k in ('x','u','z')}

def record():
    k=Q(1,100000);B=1-k;D=2-4*k;a=2-3*k
    return dict(status='corrected-full-pressure-final-mean-adapter',reserve=ledger(200002),
                tangential_orders=[str(v) for v in (B,D,a)],radial_orders=[str(v+1) for v in (B,D,a)],
                fine_angular_slot_supplied=False,pressure_or_axial_A1_removed=False,
                actual_field_inequalities_certified=False,unconditional_residual_or_blowup=False)
