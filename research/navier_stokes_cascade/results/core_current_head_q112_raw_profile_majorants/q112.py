"""Closed raw e=0 multiplier and moving-cutoff composition majorants."""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q111_closed_chart_majorant/q111.py'
s=importlib.util.spec_from_file_location('q112_charts',p);charts=importlib.util.module_from_spec(s);s.loader.exec_module(charts)
order=charts.order

def power_recipe(N,a,b,role):
    order(N);a,b=Q(a),Q(b)
    if not 0<a<=b or role not in ('f0','k0'):raise ValueError('Positive original hull and raw e=0 profile required')
    u=min(Q(1),a**3);H=max(Q(1),a**-2,b*b)
    return dict(role=role,scalar=3*H,factorial_argument=N,geometric_base=3/u,power=N,
                norm='raw-through-U',materialized=False)

def cutoff_recipe(N,a,Xb):
    order(N);a,Xb=Q(a),Q(Xb)
    if min(a,Xb)<=0:raise ValueError('Original positive radial and physical X hull required')
    u=min(Q(1),a**3)
    return dict(scalar=Q(1),factorial_argument=N,geometric_base=512*(1+Xb)/u,power=N,
                shape_sum_order=N,shape_sum_required=True,
                norm='raw-through-(U,Z,T)',slow_scale_motion_retained=True,materialized=False)

def materialize(recipe,*,shape_sum=None,maximum=40):
    order(maximum);n=recipe['factorial_argument'];order(n);order(recipe['power'])
    if max(n,recipe['power'])>maximum:raise ValueError('Keep giant profile factors symbolic')
    value=Q(recipe['scalar'])*factorial(n)*Q(recipe['geometric_base'])**recipe['power']
    if recipe.get('shape_sum_required'):
        if shape_sum is None or Q(shape_sum)<0:raise ValueError('Original normalized shape-jet sum is required')
        value*=Q(shape_sum)
    return value

class KernelProfiles:
    def __init__(self,bridge):
        if not isinstance(bridge,charts.ChartBridge):raise TypeError('Q111 identical-state chart bridge required')
        self.bridge=bridge;self.state=bridge.state
    def plan(self,n,*,kappa=Q(1,100000),target=None):
        out=self.bridge.plan(n,kappa=kappa,target=target)
        if out['source_state'] is not self.state:raise ValueError('Current profile/source state changed')
        a,b=(out['hull'][k] for k in ('Rmin','Rmax'));J=out['request']['raw_source_order']
        cp=self.bridge.radial.radial.current.current.profile
        return dict(charts=out,source_f0=power_recipe(J,a,b,'f0'),output_k0=power_recipe(n,a,b,'k0'),
                    cutoff=cutoff_recipe(n+1,a,cp.hull['Xb']),raw_product_two_power=n,
                    source_factorial_multiplicity_including_f0=3,channel=0,
                    source_state=self.state,analytic_power_profiles_derived=True,
                    actual_shape_jets_and_source_certified=False,actual_A0_gain_certified=False)
    def recipe(self,n,bundle,*,geometry,kappa=Q(1,100000),target=None):
        profiles=self.plan(n,kappa=kappa,target=target)
        old=self.bridge.recipe(n,bundle,geometry=geometry,kappa=kappa,target=target)
        return dict(profile_majorants=profiles,charted_tail=old,source_state=self.state,
                    actual_shape_source_and_full_pole_proof_still_required=True,
                    full_cone_or_residual_certified=False)

def record():
    return dict(status='closed-raw-radial-power-and-moving-cutoff-majorants',
                f0_and_k0='3*max(1,Rmin^-2,Rmax^2)*N!*(3/Umin)^N',
                moving_chi='N!*[512*(1+Xb)/Umin]^N*sum(shape[0:N+1])',
                original_shape_and_slow_motion_retained=True,source_factorial_multiplicity_including_f0=3,
                giant_profile_factors_materialized=False,actual_shape_source_or_A0_gain_certified=False,
                full_cone_pulse_residual_or_blowup=False)
