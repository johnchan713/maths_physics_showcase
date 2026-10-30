"""Same-family initial mean/projected-source assembly; all bounds conditional."""
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
    sp=importlib.util.spec_from_file_location(name,HERE.parent/path)
    m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
norm=load('q35_norm','core_current_head_q34_haar_norm_correction/haar_norms.py')
haar=load('q35_haar','core_symbolic_initial_haar_join/bounds.py')
fresh=haar.sharp.fresh
chart=fresh.seed.chart
roots=chart.bridge.primary.curl.entrance.inverse.roots

def copy_root(root):
    """Exact norm-polynomial rebinding only; preserve every original value datum."""
    callbacks={k:(lambda j,f=getattr(root,k):roots.P(f(j).c))
               for k in ('H','target','wave','curl','loghalf')}
    return roots.RootBudget(**callbacks,K0=root.K0,beta=root.beta,R=root.R,
                           mass=root.mass,width=root.width,trace=root.trace)

def build(root,stress,*,Rmin,Rmax,fixed,compact,slow,hull,fast_hull,
          max_output,c_rho,trace=None):
    """root must carry the SAME original U-chart canonical H/target/W/curl bounds.

    Fixed profiles, whole hulls, full pressure and operators remain inputs.
    Current means, scalar defects and higher-stress norm callbacks do not.
    """
    haar.order(max_output)
    class TracedPrimary(chart.ChartedPrimary):
        def covariance_R(self,n):
            if trace is not None:trace.add(('covariance',n))
            return super().covariance_R(n)
    primary=TracedPrimary(chart.bridge.WeightedPrimary(copy_root(root)),chart.ReverseChart(Rmin,Rmax))
    pole=12*(max_output+13)+18
    budget=fresh.FreshMeans(primary,fixed,compact,slow,hull,max_output+1,pole)
    initial=haar.sharp.SharpInitial(budget,max_output,fast_hull)
    def higher(n,p):
        if 3*n+3>p:raise ValueError('Higher-stress pole exceeds this common ledger')
        if trace is not None:trace.add(('higher_stress',n))
        b=norm.common.normal.transfer(n,stress(n),Rmin,Rmax)
        return haar.EP(norm.count(n)*b['higher_stress_raw'])
    return haar.InitialHaar(initial,max_output,higher,c_rho)

def record():
    return dict(result='conditional-sharp-initial-projected-source-composition',
        covariance_ceiling='m+13',common_pole='12*(m+13)+18',
        initial_mean_and_scalar_callbacks_removed=True,
        corrected_R_norms_used=True,whole_moving_pressure_retained=True,
        two_projector_channels_distinct=True,old_means_reset=False,
        fixed_operator_and_canonical_premises_certified=False,
        signed_primitive_root_or_full_source_certified=False,
        complete_residual_or_blowup=False)
