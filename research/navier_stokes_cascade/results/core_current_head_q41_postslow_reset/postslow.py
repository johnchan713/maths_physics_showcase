"""First-real same-state Haar slow reset and whole axial pressure moment."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('q41_fast',HERE.parent/'core_current_head_q40_current_fast_rows/current_fast.py')
fast=importlib.util.module_from_spec(sp);sp.loader.exec_module(fast)
refresh=fast.refresh
EP,order,convert=fast.EP,fast.order,fast.convert

def reset(n,old_haar,slow,b,moment):
    """Exact reserved-hull linear solve is a premise; outputs are remainder bounds."""
    return refresh.post.preserved.haar_slow_reset(n,old_haar,slow,b,moment)

class Postslow:
    def __init__(self,scalars,max_output):
        if not isinstance(scalars,refresh.CurrentScalars):raise TypeError('Current Q39 scalar recipe required')
        if not isinstance(scalars.fast,fast.FastFields):raise TypeError('Derived Q40 fast recipe required')
        if scalars.fast.rows.current is not scalars.current:raise ValueError('One identical current pressure state required')
        order(max_output)
        if scalars.max_output<max_output+3:raise ValueError('Whole axial moment needs slow fields through m+3')
        if scalars.fast.max_output<max_output+5:raise ValueError('Whole reset needs new fast fields through m+5')
        self.current,self.max_output=scalars,max_output
        h=scalars.initial.budget.hull;a,b=h['Rminus'],h['Rplus']
        self.moment=(b-a)*(1/a+1+b+b*b)
    def check(self,n):
        order(n)
        if n>self.max_output+1:raise ValueError('Outside postslow scalar ledger')
    def old_haar(self,n):
        order(n);self.check(max(0,n-2))
        s=self.current.initial.slow_fields(n,False)
        return {k:convert(s[j]) for k,j in [('x','X'),('u','U'),('z','Z')]}
    @lru_cache(None)
    def scalars(self,n):
        self.check(n);N=n+2;b=self.current.initial.budget
        out=reset(n,self.old_haar(N),self.current.slow(N,False),convert(b.fixed(N,2*b.pole)['b']),self.moment)
        return {k:out[k] for k in ('P','Jtheta','Jz')}
    def axial(self,n):
        order(n)
        if n>self.max_output:raise ValueError('Outside whole pressure-moment ledger')
        a=self.scalars(n+1);c=self.current.current.profile.c_rho(n+1)
        return (n+1)*(a['Jz']+c*a['P'])

def build(current,initial_A1,compact,max_output):
    order(max_output)
    rows=fast.CurrentRows(current,initial_A1,max_output+10)
    fields=fast.FastFields(rows,compact,max_output+5)
    scalars=refresh.CurrentScalars(current,fields,max_output+3)
    return Postslow(scalars,max_output)

def record():
    alpha=2-Q(3,100000)
    return dict(result='same-state-postslow-scalar-and-whole-pressure-reset',
        old_Haar_is_retained_initial_slow=True,full_old_means_still_used_in_fast_update=True,
        current_fast_slow_same_state_required=True,whole_reserved_potential_hull=True,
        postslow_order=str(min(alpha+1,2*alpha,alpha+2)),whole_axial_order=str(alpha+2),
        Q37_pressure_ceiling='m+11',Q36_target_ceiling='m+14',profile_covariance_ceiling='m+27',
        exact_zero_Haar_G0_and_five_map_moments_required=True,
        actual_high_jets_full_rows_source_residual_or_blowup_certified=False)
