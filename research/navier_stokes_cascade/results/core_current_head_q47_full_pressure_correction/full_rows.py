"""Correct Q40 full torus pressure rows; Q38 remains the Haar map."""
from functools import lru_cache
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('q47_moving',HERE.parent/'core_current_head_q46_moving_cutoff/moving.py')
moving=importlib.util.module_from_spec(s);s.loader.exec_module(moving)
post=moving.radial.provider.postslow
fast=post.fast
refresh,EP,order,convert=fast.refresh,fast.EP,fast.order,fast.convert

def compact(current,d,c0,c1,shape,time_inverse,maximum,pole,e):
    cp=current.profile
    return moving.build(d=d,Rmin=cp.Rmin,Rmax=cp.Rmax,Xa=cp.hull['Xa'],Xb=cp.hull['Xb'],
        c0=c0,c1=c1,width=current.target.width,shape=shape,time_inverse=time_inverse,
        max_order=maximum,max_pole=pole,e=e)

class FullRows(fast.CurrentRows):
    def __init__(self,current,initial_A1,max_row,*,d,c0,c1,shape,time_inverse):
        super().__init__(current,initial_A1,max_row)
        self.pressure=compact(current,d,c0,c1,shape,time_inverse,max_row+3,self.pole,0)
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
        fixed.update(Kp=self.pressure(N,self.pole)['KT'],rho=cp.weighted_rho(N),
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


def build(current,initial_A1,m,*,d,c0,c1,shape,time_inverse):
    order(m);shape=tuple(shape)
    rows=FullRows(current,initial_A1,m+10,d=d,c0=c0,c1=c1,shape=shape,time_inverse=time_inverse)
    axial=compact(current,d,c0,c1,shape,time_inverse,m+6,rows.pole,1)
    fields=fast.FastFields(rows,axial,m+5)
    return post.Postslow(refresh.CurrentScalars(current,fields,m+3),m)

def record():
    return dict(result='full-torus-pressure-correction',full_pressure_channel=0,axial_channel=1,
                full_row_pressure_request='n+3',combined_cutoff_shape_ceiling='m+14',
                original_A1_and_actual_source_required=True,complete_residual_or_blowup_certified=False)
