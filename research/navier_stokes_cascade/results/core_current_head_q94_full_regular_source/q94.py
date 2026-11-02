"""First-real regular source join; actual native contracts remain conditional."""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q93_current_weighted_wave/q93.py'
s=importlib.util.spec_from_file_location('q94_waves',p);waves=importlib.util.module_from_spec(s);s.loader.exec_module(waves)
linear=waves.linear;full=waves.full;EP=waves.EP;convert=waves.convert

def wave_wave(n,g,r,U,ratio_bound=1):
    full.order(n)
    if isinstance(ratio_bound,(float,bool)) or Q(ratio_bound)<0:raise ValueError('Exact nonnegative input-harmonic ratio bound required')
    g,r,U=map(convert,(g,r,U))
    return (1+Q(ratio_bound))*((n+1)*(g+1)+r)*U*U

class CurrentSource:
    def __init__(self,wave_inputs):
        if not isinstance(wave_inputs,waves.CurrentWaves):raise TypeError('Current original-root Q93 wave composition required')
        self.waves=wave_inputs;self.means=wave_inputs.means;self.state=wave_inputs.state
    def coefficient(self,n,linear_fixed):
        full.order(n)
        if set(linear_fixed)!={'Eik','bR','Bz','psi'}:raise ValueError('Only the four remaining full-R linear coefficient bounds required')
        N=n+2
        if N>self.means.postslow.max_output:raise ValueError('Common final-mean factory reserve n+2 required')
        w=self.waves.through(N)
        if w['wave_state'] is not self.state or self.means.state is not self.state:raise ValueError('One identical first-real source state required')
        fields=self.means.fields(N);M=sum(fields.values(),EP())
        budget=self.means.current.initial.budget;a=budget._data(N,budget.pole)
        g,r=convert(a['graph']),convert(a['r'])
        matrix=sum((Q(5**j,factorial(j))*self.waves.jet('matrix',j) for j in range(N+1)),EP())
        factor=2*convert(self.state.primary.chart.transfer(N))*16*188**N*factorial(N)
        K=factor*waves.rescale(matrix)
        data={k:convert(a[k]) for k in ('b','G','c')}
        data.update(g=g,r=r,K=K,normal=w['normal'],**{k:w[k] for k in ('T','R','Pi')})
        data.update({k:convert(v) for k,v in linear_fixed.items()})
        lin=linear.envelope(n,data)['F']
        ww=wave_wave(n,g,r,w['U'])
        mw=linear.meanwave.envelope(n,M,w['U'],g,w['normal'],r)
        return dict(F=lin+ww+mw,linear=lin,wave_wave=ww,mean_wave=mw,
                    source_order=n,input_order=N,common_pole=2*w['pole'],wave_pole=w['pole'],
                    mean_factory_reserve=N,chart='factorial-through-R',stage='first_real',
                    actual_native_contracts_certified=False,complete_residual_improvement=False)

def record():
    return dict(status='same-state-complete-first-real-regular-source-join',
                source_order='1-kappa',source_pole='2*p_(j+2)',
                input_order_shift=2,mean_factory_reserve_shift=2,
                current_target_reserve_shift=16,profile_reserve_shift=29,
                output_harmonic_division=False,all_ordered_wave_products_retained=True,
                native_provider_inequalities_certified=False,full_residual_or_blowup=False)
