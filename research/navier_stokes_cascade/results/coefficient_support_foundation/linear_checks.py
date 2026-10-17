"""Actual linear-map quadrature diagnostics on manufactured, non-small data."""
from functools import lru_cache
import mpmath as mp
import sys
from support_bounds import HERE,load

folder=HERE.parent/'annular_first_correction'
_original_search_path=sys.path[:]
sys.path.insert(0,str(folder))
try:
    extension=load('coefficient_support_actual_extension',folder/'extension.py')
    fixtures=load('coefficient_support_actual_fixture_integrands',folder/'fixtures.py')
finally:sys.path[:]=_original_search_path

@lru_cache(None)
def numerical_review():
    rows=[]
    with mp.workdps(80):
        direct_rule=fixtures.Rule(32,8)
        for lam in ('.001','.01'):
            repair=extension.Repair(lam,order=24,panels=6)
            radius=mp.mpf(3);e_star=mp.mpf(7)
            def incoming(t):return mp.matrix([1+t+t*t,2-t/3,3+t*t*t,4-t+t*t,5+t/2])
            def amplitude(t):return e_star/(1+t*t)
            def root(t):return extension.solve_patch(incoming(t),radius,amplitude(t),repair)
            for eta in map(mp.mpf,('-.5','0','.4','1')):
                c=root(eta);changes=fixtures.direct_patch_changes(c,repair,radius,amplitude(eta),direct_rule)
                closed=incoming(eta)+changes
                norm=max(1,max(map(abs,c)))
                error=max(map(abs,extension.normalized_target(closed,radius,amplitude(eta),repair.lam)))/norm
                target=extension.normalized_target(incoming(eta),radius,amplitude(eta),repair.lam)
                nonlinear_gap=max(map(abs,repair.value(c)+target))/norm
                # Differentiate the ordinary inputs and e(eta)^-1 BEFORE solving.
                second=mp.matrix([mp.diff(lambda t:root(t)[i],eta,2) for i in range(5)])
                dz=mp.matrix([mp.diff(lambda t:extension.normalized_target(incoming(t),radius,amplitude(t),repair.lam)[i],eta,2) for i in range(5)])
                jet_error=max(map(abs,repair.B*second+dz))/norm
                frozen=mp.matrix([mp.diff(lambda t:extension.normalized_target(incoming(t),radius,amplitude(eta),repair.lam)[i],eta,2) for i in range(5)])
                freeze_gap=max(map(abs,dz-frozen))
                rows.append(dict(lambda_value=lam,eta=str(eta),normalized_closure_per_coefficient=float(error),
                    second_jet_equation_per_coefficient=float(jet_error),nonlinear_repair_gap_per_coefficient=float(nonlinear_gap),
                    frozen_amplitude_second_jet_gap=float(freeze_gap),coefficient_max=float(norm)))
    return dict(rows=rows,checks=dict(
        independent_ordinary_integrands_close_all_five_moments=all(r['normalized_closure_per_coefficient']<1e-12 for r in rows),
        same_linear_inverse_controls_second_eta_jets=all(r['second_jet_equation_per_coefficient']<1e-60 for r in rows),
        nonlinear_leading_repair_rejected_for_positive_order=all(r['nonlinear_repair_gap_per_coefficient']>1 for r in rows),
        freezing_actual_amplitude_derivatives_rejected=all(r['frozen_amplitude_second_jet_gap']>1 for r in rows),
        no_small_target_condition=all(r['coefficient_max']>1 for r in rows),
        manufactured_data_not_actual_selected_profile_or_uniform_proof=True))
