"""Current nonzero-Haar row decomposition and complete fast norm adapter."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('q40_refresh',HERE.parent/'core_current_head_q39_current_scalar_refresh/refresh_scalars.py')
refresh=importlib.util.module_from_spec(sp);sp.loader.exec_module(refresh)
EP,order,convert=refresh.EP,refresh.order,refresh.convert

class CurrentRows:
    """Initial exact A1 tail and same-domain fixed inputs remain conditional."""
    def __init__(self,current,initial_A1,max_row):
        if not isinstance(current,refresh.maps.NewPoleMaps):raise TypeError('Current Q38 pressure state required')
        order(max_row)
        if current.state.max_output<max_row+1:raise ValueError('Rows need covariance through n+3')
        if not callable(initial_A1):raise TypeError('Original initial A1-tail norm recipe required')
        self.current,self.initial_A1,self.max_row=current,initial_A1,max_row
        self.initial=current.target.source.haar.initial
        N=max_row+3;p=self.initial.budget.pole;d=current.target.pole
        self.q=max(6*p,6*N+9,6*N+d+6)
        self.pole=2*self.q
    @lru_cache(None)
    def rows(self,n):
        order(n)
        if n>self.max_row:raise ValueError('Outside declared current-row ledger')
        N=n+3;b=self.initial.budget
        f={k:convert(v) for k,v in b._fast(N).items()}
        s0=self.initial.slow_fields(N,True)
        s={k:convert(s0[j]) for k,j in [('x','X'),('u','U'),('z','Z')]}
        fixed={k:convert(v) for k,v in b._data(N,b.pole).items()}
        cp=self.current.profile;a,z=cp.Rmin,cp.Rmax
        fixed.update(Kp=self.current.recipe(N,self.pole)['Kp'],rho=cp.weighted_rho(N),
            Kmom=EP((z-a)*(1/a+1+z+z*z)*refresh.maps.update.signed.maps.edge_absorption(self.pole,cp.hull['t1'])))
        zero={k:EP() for k in ('x','u','z')}
        fast=refresh.post.fast.difference(n,zero,f,fixed)
        slow=refresh.post.fast.difference(n,f,s,fixed)
        cov=self.current.state.covariance(N)
        if cov['pole']>self.pole:raise ValueError('Complete covariance exceeds common row pole')
        wave=refresh.post.fast.covariance_rows(n,cov['complete'],fixed)
        tail=self.initial_A1(n,self.pole)
        if set(tail)!={'order','pole','Ez'}:raise ValueError('Original axial A1 tail and coverage required')
        order(tail['order']);order(tail['pole'])
        if tail['order']<n or tail['pole']>self.pole:raise ValueError('Initial A1 estimate does not embed in this row domain')
        return dict(Etheta=2*(fast['Etheta']+slow['Etheta']+wave['Etheta']),
            Ez=2*(fast['Ez']+slow['Ez']+wave['Ez']+convert(tail['Ez'])),
            old_fast=fast,old_slow=slow,signed_wave=wave,initial_A1=convert(tail['Ez']),pole=self.pole)

def fast_recipe(n,rows,compact,pole):
    """Original inverse, including physical S and actual axial 1+KA1."""
    order(n);order(pole)
    a=compact(n,pole);b=compact(n+1,pole)
    for v,j in ((a,n),(b,n+1)):
        if set(v)!={'order','pole','KT','KA','I'}:raise ValueError('Full compact/time inverse contract required')
        order(v['order']);order(v['pole'])
        if v['order']<j or v['pole']<pole:raise ValueError('Compact inverse domain too short')
    v=refresh.post.fast_envelope(n,rows,convert(b['KT']),convert(a['KA']),EP({1:1})*convert(a['I']))
    return dict(v,order=n,pole=pole)

class FastFields:
    def __init__(self,rows,compact,max_output):
        if not isinstance(rows,CurrentRows):raise TypeError('Derived current-row decomposition required')
        order(max_output)
        if rows.max_row<max_output+5:raise ValueError('Fast inverse needs current rows through n+5')
        self.rows,self.compact,self.max_output=rows,compact,max_output
    @lru_cache(None)
    def __call__(self,n):
        order(n)
        if n>self.max_output:raise ValueError('Outside new fast-field ledger')
        return fast_recipe(n,self.rows.rows,self.compact,self.rows.pole)

def record():
    return dict(result='current-rows-to-complete-fast-norm-composition',
        initial_A1_retained=True,full_old_fast_slow_differences=True,
        full_signed_covariance_and_pressure=True,nonzero_Haar_factor=2,
        fast_axial_A1_velocity_retained=True,time_inverse_S_retained=True,
        Q39_slow_output_m_requires_Q37_pressure='m+8',Q36_target='m+11',Q36_profile_covariance='m+24',
        original_A1_compact_time_and_background_contracts_required=True,
        actual_fast_source_postslow_or_full_residual_certified=False,blowup=False)
