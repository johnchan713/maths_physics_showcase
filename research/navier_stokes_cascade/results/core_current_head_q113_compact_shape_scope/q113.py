"""Compact-shape quantifiers and role binding; no replacement cutoff choice."""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q112_raw_profile_majorants/q112.py'
s=importlib.util.spec_from_file_location('q113_profiles',p);profiles=importlib.util.module_from_spec(s);s.loader.exec_module(profiles)
order=profiles.order
ROLES={
 'compact-X':dict(argument='X=R^2/(2s)',source='core_current_head_q30_fixed_support_join/README.md'),
 'borel':dict(argument='a_n*q',source='core_borel_cutoff_certificate/README.md'),
 'variance':dict(argument='8*(v-2)/delta-1',source='core_current_head_q66_cutoff_amplitude/README.md')}

def require_compact(role):
    if role!='compact-X':raise ValueError('Borel/variance formulas do not identify the original compact-X shape')
    return dict(ROLES[role])

def transition_comparison(c0,c1,index):
    """Comparison class only, never a selected actual cutoff."""
    c0,c1=Q(c0),Q(c1);order(index)
    if not c0<c1 or index<1:raise ValueError('Ordered fixed collars and positive comparison index required')
    gap=c1-c0;width=gap/(4*index);center=(c0+c1)/2
    return dict(center=center,half_width=width,transition_length=2*width,
                derivative_lower_bound=1/(2*width),plateaus_retained=True,
                chosen_actual_shape=False)

class CompactShapeScope:
    def __init__(self,kernels):
        if not isinstance(kernels,profiles.KernelProfiles):raise TypeError('Q112 identical-state raw profiles required')
        self.kernels=kernels;self.state=kernels.state
    def plan(self,n,*,role='compact-X',kappa=Q(1,100000),target=None):
        lineage=require_compact(role);out=self.kernels.plan(n,kappa=kappa,target=target)
        if out['source_state'] is not self.state:raise ValueError('Current compact profile state changed')
        N=out['cutoff']['shape_sum_order']
        return dict(profile_majorants=out,lineage=lineage,shape_order=N,
                    shape_constant='sum_{j=0}^N sup_[c0,c1] |hat_chi^(j)|/j!',
                    conditional_on_original_fixed_Cinfinity_choice=True,
                    qualitative_fixed_order_band_uniformity=True,
                    uniform_over_arbitrary_cutoff_choices=False,
                    numerical_shape_value_supplied=False,numerical_shape_bound_verified=False,
                    source_state=self.state,actual_source_or_A0_gain_certified=False)
    def recipe(self,n,bundle,*,geometry,role='compact-X',kappa=Q(1,100000),target=None):
        scope=self.plan(n,role=role,kappa=kappa,target=target)
        out=self.kernels.recipe(n,bundle,geometry=geometry,kappa=kappa,target=target)
        return dict(scope=scope,kernel_recipe=out,source_state=self.state,
                    raw_kernel_profile_finiteness_conditional_on_fixed_choice=True,
                    actual_high_jet_source_still_required=True,full_cone_or_residual_certified=False)

def record():
    return dict(status='compact-shape-identity-and-fixed-order-uniformity-audit',
                qualitative_profile_finiteness_already_follows_from_original_fixed_smooth_choice=True,
                finite_constants_depend_on_selected_shape_and_order=True,
                geometry_alone_gives_no_numeric_derivative_bound=True,
                borel_or_variance_substitution_rejected=True,
                actual_numeric_shape_or_high_jet_source_certified=False,
                original_cutoff_unchanged=True,full_cone_pulse_residual_or_blowup=False)
