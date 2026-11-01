"""Support-aware adapter to the historical finalized radial tuple recipe."""
from pathlib import Path
import importlib.util
from fractions import Fraction as Q
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('q80_radial',HERE.parent/'core_background_radial_reconstruction/budgets.py')
radial=importlib.util.module_from_spec(spec);spec.loader.exec_module(radial)

def clipped(h,n,m,r,H,V,Frows,Urows,*,zero_total):
    if zero_total is not True:raise ValueError('Exact zero total axial moment is required')
    if type(n) is not int or n<1:raise ValueError('Only positive-order coefficients have this support')
    h,r,H,V=map(radial.exact,(h,r,H,V))
    if not 0<r<=H or V<=0:raise ValueError('Positive ordered chart and support required')
    # Validate all original order and row requirements, even in the zero branch.
    full=radial.reconstruct(h,n,m,2*H*H,r/2,2*H,Frows,Urows)
    a=r*r/8;b=2*H*H
    if a>=V:return dict(bound=Q(0),native_upper=None,exterior=True)
    end=min(b,V)
    result=radial.reconstruct(h,n,m,end,r/2,2*H,Frows,Urows)
    assert result['through'][m]<=full['through'][m]
    return dict(bound=result['through'][m],native_upper=end,exterior=False)

def record():
    f=(Q(7),Q(11));u=(Q(5),Q(17));h=Q(1,200)
    full=radial.reconstruct(h,1,1,200,Q(1,2),20,f,u)['through'][1]
    cut=clipped(h,1,1,1,10,3,f,u,zero_total=True)
    return {'status':'support-saturated-positive-coefficient-tuple',
      'fixture_full':str(full),'fixture_clipped':str(cut['bound']),
      'fixture_head':str(max(1,9856*(1+2*h)*cut['bound'])),
      'original_band_or_schedule_changed':False,'leading_tuple_clipped':False,
      'pressure_tail_removed':False,'unconditional_pulse_or_residual_or_blowup':False}
