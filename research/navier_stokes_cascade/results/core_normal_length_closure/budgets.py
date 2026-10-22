"""Original radius-free normal floor from the retained length budget."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


frame=load('normal_length_frame',HERE.parent/'frame_column_foundation/frame_budget.py')
selector=load('normal_length_selector',HERE.parent/'core_primary_tail_recheck/proof_checks.py')
exact=selector.exact


def frame_inputs(p):
    if set(p)!=set(frame.inputs()):raise ValueError('Complete ORIGINAL primitive bounds required; no diagnostic defaults')
    return frame.inputs(**p)


def normal_length_certificate(p):
    """Algebraic implication for valid original bounds, before an ODE solve.

    Physical consistency follows from lambda0^2<=2*F0*|g0|. A set of
    positive manufactured numbers need not be physically consistent.
    """
    p=frame_inputs(p);b=frame.budgets(p);sr=b['Sref'];Lric=b['Riccati_acceptance']
    consistency=p['lam_min']**2<=2*p['Fmax']*p['gmax']
    margins=dict(operator_lower=b['CE']-2*p['lam_min']*b['Da'],
                 normal_length=Lric-16*(1+sr)*b['Da'],
                 radial_tilt_length=Lric-8*(1+sr)*b['Ds'])
    return dict(original_budgets=b,physical_reference_consistency=consistency,proof_margins=margins,
                length_rule_floor=max(Q(25),Lric),normal_threshold=b['normal_acceptance'],
                gates=dict(physical_reference_consistency=consistency,
                           operator_and_length_implications=consistency and all(v>=0 for v in margins.values()),
                           original_Riccati_and_25_dominate_original_normal_threshold=consistency and max(Q(25),Lric)>=b['normal_acceptance']),
                inverse_normal_square=4/p['bmin']**2,
                floor_constant_is_original_Bs_lower_not_phase_beta=True,
                source_budgets_selected_or_changed=False)


def tangential_error(p,S,L,epsilon,epsilon2,k_inverse):
    """Five complete component-sum terms, using the ORIGINAL pulse length."""
    p=frame_inputs(p);S,L,e,e2,ki=map(exact,(S,L,epsilon,epsilon2,k_inverse))
    if S<1 or not 1<=L<=S or min(e,e2,ki)<0:raise ValueError('Original 1<=L<=S and nonnegative slow errors required')
    b=frame.primitive(p);r=p['Rmin'];B=p['bmax'];u=p['u'];g=p['gmin']
    terms=dict(angular_box=B*p['Dbox']/(r*S**3),
               angular_tilt=B*u*p['Rmax']/(r*g*L),angular_rounding=ki/r,
               axial_tilt=B*u/(g*L),complete_axial_transport=e*L*b['Zc'])
    slow=S*S*(e+e2+ki)<=1
    return dict(terms=terms,error_upper=sum(terms.values(),Q(0)),Dt_over_L=b['Dt']/L,
                original_slow_gate=slow,implication_valid=not slow or sum(terms.values(),Q(0))<=b['Dt']/L)


def original_selector_data(p):
    if set(p)!=selector.KEYS:raise ValueError('All original selector inputs required')
    p={k:exact(v) for k,v in p.items()}
    if min(p.values())<=0 or not 0<p['h']<Q(1,100) or p['epsilon']>1 or p['qstar']>1 or p['radius']>Q(1,6057072) or p['slow_cost']<1 or p['slow_tolerance']>1:raise ValueError('Original construction and selector ranges required')
    tilt=selector.inherited.tilt_certificate(p['epsilon'])
    gamma=p['omega']**2*tilt['chosen_d0']/(9*tilt['kappa'])
    B0=min(p['rb']/8,p['omega']**2*tilt['chosen_d0']/(16*p['lpsi']),gamma/(128*p['cpert']))
    m=selector.inherited.ceil_log2(p['cost']/B0)
    slow=Q(3*factorial(8)*4**8)*p['slow_cost']/(p['slow_tolerance']*p['h']**8)
    return p,gamma,m,slow


def verify_existing(p,selector_inputs,existing_ell):
    """Verify one supplied old band and its original length SEMANTICS.

    No invocation of select_tail, no new band or length is returned.
    Success is a conditional implication from valid actual input bounds,
    not certification that the supplied rational fixtures bound a fluid.
    """
    p=frame_inputs(p);data=normal_length_certificate(p)
    sp,gamma,m,slow=original_selector_data(selector_inputs)
    checks=selector.band_checks(sp,existing_ell,gamma,m,slow)
    gates=dict(data['gates'])
    gates.update({name:bool(v) for name,v in checks.items()})
    gates.update(original_length_has_its_declared_Riccati_and_25_meaning=sp['length']>=data['length_rule_floor'],
                 original_geometric_length_upper_is_below_S=12*sp['radius']<1)
    return dict(gates=gates,conditional_normal_floor_implication_passes=all(gates.values()),
                existing_ell=existing_ell,existing_length_budget=sp['length'],
                original_budget=data,actual_pulse_length_lower=2*sp['radius']*existing_ell**2,
                inverse_normal_square=data['inverse_normal_square'],
                normal_nonvanishing_on_rule_defined_accepted_tail=all(gates.values()),
                actual_input_enclosures_or_arbitrary_historical_band_certified=False,
                stronger_phase_beta_bound_certified=False,full_propagator_or_residual_certified=False,
                bands_lengths_supports_or_representatives_selected_or_changed=False)
