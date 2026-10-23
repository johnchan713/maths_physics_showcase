"""Deferred seam-trace requirements; no canonical trace enclosure is supplied."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('trace_axial_parent',HERE.parent/'core_axial_potential_pressure_gain/bounds.py')
parent=importlib.util.module_from_spec(spec);spec.loader.exec_module(parent)
spec=importlib.util.spec_from_file_location('trace_original_geometry',HERE.parent/'auxiliary_geometry_audit/geometry.py')
geometry=importlib.util.module_from_spec(spec);spec.loader.exec_module(geometry)

def exact(x):
    if isinstance(x,bool) or not isinstance(x,(int,Q)):raise TypeError('Use exact integers or fractions')
    return Q(x)

def order(n):
    if isinstance(n,bool) or not isinstance(n,int) or n<0:raise ValueError('Nonnegative integer order required')
    return n

def transported_modes(coefficients,relative_cover):
    """Finite exact fixture on ORIGINAL relative overlap covers <=3.

    Pi_i f(Y)=f(J^i Y) sends k to (J^i)^T k. Factor the common Pi_i
    before comparing; never build the huge band-dependent absolute cover.
    """
    order(relative_cover)
    if relative_cover>geometry.DELTA_MAX:raise ValueError('Use the original relative overlap cover <=3')
    A=geometry.matrix_power(relative_cover);out={}
    for k,value in coefficients.items():
        if not isinstance(k,tuple) or len(k)!=2 or any(type(j) is not int for j in k):raise ValueError('Integer torus mode pair required')
        ell=tuple(sum(A[j][i]*k[j] for j in range(2)) for i in range(2))
        out[ell]=out.get(ell,Q(0))+exact(value)
    return {k:v for k,v in out.items() if v}

def boundary_recipe(p):
    order(p)
    if p==0:raise ValueError('Boundary integration formula requires p>=1')
    return dict(integrations=p,trace_first=0,trace_last=p-1,
        term='exp(i omega s) (-1)^r Delta_r(s)/(i omega)^(r+1)',
        Delta='left minus right of weighted Fourier coefficient, including exterior zero pieces',
        bulk='(-1)^p (i omega)^(-p) sum integral F_piece^(p) exp(i omega U)',
        zero_trace_through_for_M_minus_p=p-2 if p>=2 else None,
        last_trace_may_be_bounded_at_bulk_order=True,materialized=False)

def seam_orders(source_order,p,trace_rates,*,kappa=Q(1,100000),target=None):
    """Fixed-mode epsilon arithmetic only; None means an EXACT zero trace.

    A finite mapping r->extra epsilon order does not supply all other traces.
    Require a complete mapping for this finite diagnostic, not giant p.
    """
    a=exact(source_order);k=exact(kappa);order(p)
    if p<1 or not 0<k<Q(1,4):raise ValueError('Positive p and original kappa interval required')
    if set(trace_rates)!=set(range(p)):raise ValueError('Declare every trace r=0,...,p-1, or use the deferred recipe')
    leaks={r:a+exact(b)+(r+1)*k for r,b in trace_rates.items() if b is not None}
    floor=min([a+p*k]+list(leaks.values()))
    t=None if target is None else exact(target)
    return dict(bulk=a+p*k,leaks=leaks,guaranteed_power=floor,
        strict_target_pass=None if t is None else floor>t,
        exact_phase_cancellation_assumed=False,actual_trace_inputs_certified=False)

def A1_trace_request(*,physical_power=None,h=Q(1,200),kappa=Q(1,100000)):
    """Zeroth-output physical route only. All higher-output requests stay old."""
    a=parent.A1_request(0,physical_power=physical_power,h=h,kappa=kappa)
    p=a['integration_order']
    return dict(stage='same postwave',output_order=0,source='A1 Dz Q Phi',
        source_order=Q(2),target=a['normalized_target'],integrations=p,
        direct_bound_without_nonstationary_IBP=p==0,
        matching_Dz_weighted_Phi_through=p-2 if p>=2 else None,
        sufficient_matching_weighted_Phi_through_for_moving_seam=p-1 if p>=1 else None,
        original_raw_Phi_order=a['raw_Phi_order'],
        boundary_recipe=boundary_recipe(p) if p else None,
        higher_output_parent_graph_retained=True,
        actual_closed_cover_trace_identities_verified=False)

def certificate():
    return dict(parent_commit='c9269c3412d4284b4c0b72ad5dd6663427e285ad',
        projection_defect='A_e D_e h=-R^-e chi_R sum internal [r^e h]_right-minus-left on the characteristic, with zero exterior endpoint values',
        primitive_defect='T_e D_e h=h-R^-e(sum_(s<R) jump_s-chi sum_all jump_s), in fixed characteristic coordinates',
        pressure='Q pm=Q K1+T0 Q r1 plus the primitive defect of Q K1 if classical pieces fail to glue; no change to the smooth global identity',
        axial='A1 Q Ez=A1 Dz Q Phi plus A1 D1 Q Omega if Omega fails its value trace contract; Omega=Lambda+2G Psi',
        moving_seam='partial_Z Delta_r F=Delta_r(partial_Z F)+s_Z Delta_(r+1) F',
        normalized=A1_trace_request(),physical_illustration=A1_trace_request(physical_power=0),
        physical_scope='zeroth output only; sufficient separate trace norm test, exact phase cancellation left open; not a summed all-jet canonical bound',
        unchanged='full pm, moving rho P, K0/K1, Pi1/Pi2, six covariance entries, two viscosities, G_R, cutoff commutators and full older norm/domain graph',
        actual_input_validation_complete_regular_decay_or_blowup_proved=False)

def certify_actual_traces():
    raise NotImplementedError('Missing canonical whole closed-cover mixed trace identities for Q Omega, Q K1 and weighted partial_Z Q Phi, with original chart/cover transport and all requested norm constants')
