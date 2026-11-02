"""Closed all-fixed-order chart constants; no giant source jets are evaluated."""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q110_radial_tail_request/q110.py'
s=importlib.util.spec_from_file_location('q111_radial',p);radial=importlib.util.module_from_spec(s);s.loader.exec_module(radial)

def order(n):
    if type(n) is not int or n<0:raise ValueError('Nonnegative integer derivative order required')

def hull(a,b,d):
    a,b,d=map(Q,(a,b,d))
    if not 0<a<=b or not Q(1,3)<d<Q(5,2):raise ValueError('Original positive hull and radial exponent range required')
    return dict(Rmin=a,Rmax=b,d=d,Umin=min(Q(1),a**3),Umax=max(Q(1),b**3))

def forward_recipe(N,H,Umin):
    order(N);H,Umin=Q(H),Q(Umin)
    if min(H,Umin)<=0:raise ValueError('Same-hull positive coordinate bounds required')
    return dict(scalar=486,factorial_argument=N,geometric_base=16*(1+H)/min(Q(1),Umin),power=N,
                input='factorial-through-R',output='factorial-through-U',materialized=False)

def reverse_recipe(N,a,b):
    order(N);a,b=Q(a),Q(b)
    if not 0<a<=b:raise ValueError('Positive original R hull required')
    B=max(Q(1),b*b,1/a)
    return dict(scalar=2,factorial_argument=N,geometric_base=16*(1+B)/min(Q(1),a),power=N,
                input='raw-U',output='factorial-through-R',materialized=False)

def materialize(recipe,*,maximum=40):
    """Moderate diagnostic only; the analytical recipe has no order cutoff."""
    order(maximum);n=recipe['factorial_argument'];order(n)
    order(recipe['power'])
    if n>maximum or recipe['power']>maximum:raise ValueError('Keep the high-order factorial/power symbolic')
    return Q(recipe['scalar'])*factorial(n)*Q(recipe['geometric_base'])**recipe['power']

class ChartBridge:
    def __init__(self,radial_audit):
        if not isinstance(radial_audit,radial.RadialAudit):raise TypeError('Q110 identical-state radial audit required')
        self.radial=radial_audit;self.state=radial_audit.state
    def plan(self,n,*,kappa=Q(1,100000),target=None):
        self.radial.check()
        if self.radial.state is not self.state:raise ValueError('Current physical state changed')
        c=self.radial.radial.current;cp=c.current.profile
        h=hull(cp.Rmin,cp.Rmax,c.fast.rows.pressure.d)
        req=radial.request(n,kappa=kappa,target=target);J=req['raw_source_order']
        return dict(request=req,hull=h,source_chart=forward_recipe(J,h['Rmax'],h['Umin']),
                    output_chart=reverse_recipe(n,h['Rmin'],h['Rmax']),
                    additional_raw_source_factorial=J,source_factorial_multiplicity=2,
                    source_state=self.state,chart_majorants_derived=True,
                    actual_source_or_kernel_profiles_certified=False,radial_tail_certified=False)
    def recipe(self,n,bundle,*,geometry,kappa=Q(1,100000),target=None):
        charts=self.plan(n,kappa=kappa,target=target)
        out=self.radial.recipe(n,bundle,geometry=geometry,kappa=kappa,target=target)
        if out['source_state'] is not self.state:raise ValueError('Same current source/geometry required')
        return dict(raw_tail=out,charts=charts,source_state=self.state,
                    original_source_unweighting_and_output_floor_retained=True,
                    raw_kernel_profiles_and_high_jet_bounds_still_required=True,
                    full_cone_or_actual_pulse_certified=False)

def record():
    return dict(status='closed-same-hull-chart-majorants-for-current-radial-request',
                forward='486*N!*[16*(1+Rmax)/min(1,Umin)]^N',
                reverse='2*N!*[16*(1+max(1,Rmax^2,1/Rmin))/min(1,Rmin)]^N',
                source_jet_order=100005,source_factorial_multiplicity=2,
                giant_factorials_or_source_jets_materialized=False,
                actual_source_kernel_or_radial_tail_certified=False,
                full_cone_pulse_residual_or_blowup=False)
