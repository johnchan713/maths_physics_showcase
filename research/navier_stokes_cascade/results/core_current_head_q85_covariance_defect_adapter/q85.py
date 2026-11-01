"""Same-root covariance coefficient join, conditional on Q51's physical rate."""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util
_here=Path(__file__).resolve().parent
_spec=importlib.util.spec_from_file_location('q85_full_rows',_here.parent/'core_current_head_q47_full_pressure_correction/full_rows.py')
full=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(full)
state_module=full.refresh.maps.update

def defect(state,n,*,target_rate):
    if not isinstance(state,state_module.SignedCovariance):raise TypeError('Original Q37 same-root state required')
    full.order(n)
    if target_rate!=Q(2) or isinstance(target_rate,(bool,float)):
        raise ValueError('Q51 order-two norm premise required; do not regrade the old ledger by assertion')
    if state.target.root is not state.primary.primary.root:raise ValueError('Original signed and primary root must coincide')
    if n>state.max_output:raise ValueError('Conservative covariance/pressure reserve exceeded')
    v=state.coefficients(n);cov=state.covariance(n)
    transfer=state_module.EP(state.primary.chart.transfer(n))
    conversion=243 # sum_(j<=n) 5^j/j! <= exp(5) < 3^5
    curl=conversion*transfer*2**n*v['B']*v['B']
    principal=conversion*transfer*2**n*v['A']*v['A']
    signed=conversion*cov['complete']
    Z=curl+signed
    return dict(W=principal+Z,Z=Z,oldcurl=curl,signed=signed,raw_to_factorial=conversion,pole=cov['pole'],order=n,
                target_order_required=n+3,physical_rate_premise='Q51 on the same actual target')

def record():
    k=Q(1,100000)
    return dict(status='same-root-full-covariance-defect-adapter',
                defect_term_rates=[str(v) for v in (2,3,2-2*k,3-2*k,4-2*k)],
                Q84_covariance_request=200002,conservative_target_reserve=200005,
                old_target_rate_suffices=False,physical_wave_rebuilt=False,
                actual_field_inequalities_certified=False,unconditional_residual_or_blowup=False)
