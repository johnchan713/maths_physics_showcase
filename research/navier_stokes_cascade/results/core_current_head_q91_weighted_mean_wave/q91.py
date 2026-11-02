"""Weighted mean-wave source envelope from corrected unweighted means."""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q86_final_mean_adapter/q86.py'
s=importlib.util.spec_from_file_location('q91_means',p);means=importlib.util.module_from_spec(s);s.loader.exec_module(means)
full=means.full;EP=full.EP;convert=full.convert

def envelope(n,M,U,g,N,r,harmonic=1):
    full.order(n)
    if type(harmonic) is not int or harmonic==0:raise ValueError('Fixed nonzero integer harmonic required')
    M,U,g,N,r=map(convert,(M,U,g,N,r))
    return 3*(2*(n+1)*(g+1)+2*abs(harmonic)*N+2*r)*M*U

class MeanWave:
    def __init__(self,postslow):self.means=means.FinalMeans(postslow)
    def coefficient(self,n,*,wave,normal,wave_state,harmonic=1):
        full.order(n);f=self.means
        if wave_state is not f.state:raise ValueError('Weighted wave must refer to the identical first-real state')
        if means.FinalMeans(f.postslow).state is not f.state:raise ValueError('Changed current pressure state')
        N=n+1
        fields=f.fields(N);M=sum(fields.values(),EP())
        budget=f.current.initial.budget;data=budget._data(N,budget.pole)
        return dict(coefficient=envelope(n,M,wave,data['graph'],normal,data['r'],harmonic),
                    source_order=n,field_and_wave_order=N,factory_reserve=N,
                    chart='factorial-through-R',mean_weight_required=0,
                    wave_weight_required='sqrt(zeta)*P',actual_inputs_certified=False)

def record():
    return dict(status='weighted-first-real-mean-wave-from-unweighted-means',
                mean_component_norm_callback_removed=True,extra_field_jets=1,
                mean_radial_edge_weight_required=0,source_radial_weight='1/2',
                epsilon_order='1-kappa',weighted_wave_input_still_required=True,
                actual_full_source_or_complete_residual_certified=False,unconditional_blowup=False)
