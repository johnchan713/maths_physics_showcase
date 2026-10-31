"""Q29 weighted raw compact constants composed with Q40's finite interface."""
from pathlib import Path
import importlib.util
from functools import lru_cache
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('q42_postslow',HERE.parent/'core_current_head_q41_postslow_reset/postslow.py')
postslow=importlib.util.module_from_spec(s);s.loader.exec_module(postslow)
fast=postslow.fast
EP,order,convert=fast.EP,fast.order,fast.convert
KEYS={'f','g','k','chi','one_minus_chi','chi_U','jacobian','length','edge','kernel','interior'}

def constants(n,data):
    order(n)
    if set(data)!=KEYS:raise ValueError('All Q29 raw multiplier and weighted geometry bounds required')
    v={k:convert(x) for k,x in data.items()}
    a=2**(3*n)*v['f']*v['jacobian']
    return dict(KT=a*v['g']*(v['chi']+v['one_minus_chi'])*v['kernel'],
                KA=a*v['k']*v['chi_U']*v['length']*v['edge']*v['interior'])

class CompactProvider:
    """Lower raw-U, same-profile Q29 hypotheses remain analytical contracts.

    No monotonicity in pole is assumed. Each request is evaluated at its
    exact pole. `I` bounds the original time inverse divided by S; Q40
    inserts S. Input fields must satisfy Q29's weighted raw derivative
    hypothesis; this adapter is not a native/chart norm conversion.
    """
    def __init__(self,multipliers,time_inverse,max_order,max_pole):
        order(max_order);order(max_pole)
        if not callable(multipliers) or not callable(time_inverse):raise TypeError('Lower recipes required')
        self.multipliers,self.time_inverse=multipliers,time_inverse
        self.max_order,self.max_pole=max_order,max_pole
    @lru_cache(None)
    def __call__(self,n,pole):
        order(n);order(pole)
        if n>self.max_order or pole>self.max_pole:raise ValueError('Outside declared compact archive')
        d=self.multipliers(n,pole)
        if set(d)!={'order','pole','cutoff_order','bounds'}:raise ValueError('Declared raw input and cutoff coverage required')
        for k in ('order','pole','cutoff_order'):order(d[k])
        if d['order']<n or d['pole']!=pole or d['cutoff_order']<n+1:raise ValueError('Raw jets, exact pole or cutoff jets insufficient')
        c=constants(n,d['bounds'])
        return dict(c,I=convert(self.time_inverse(n,pole)),order=n,pole=pole)

def record():
    return dict(result='explicit-current-compact-provider',fast_output_m_requires_compact='m+1',
        fast_output_m_requires_cutoff='m+2',Q41_output_m_requires_compact='m+6',
        Q41_output_m_requires_cutoff='m+7',exact_pole_requested=True,
        raw_input_chart_conversion_and_actual_profile_hypotheses_required=True,
        initial_A1_tail_full_source_residual_or_blowup_certified=False)
