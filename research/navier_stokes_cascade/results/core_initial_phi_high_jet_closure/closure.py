"""Finite-order induction for the corrected first-real post-signed Phi input."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "a1_request_graph", HERE.parent/"core_initial_a1_request_graph/requests.py")
request_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(request_module)


PRIMITIVES = (
    "analytic_background_and_Borel_coefficients",
    "same_band_phase_normal_and_value_cone_gaps",
    "canonical_homogeneous_ODE_and_full_path_cutoffs",
    "fixed_profile_rho_and_rho_squared_projectors",
    "fixed_compact_temporal_and_slow_inverse_maps",
)


DEPENDENCIES = {
    "primary_root_and_complete_curl": PRIMITIVES[:3],
    "primary_covariance": ("primary_root_and_complete_curl",),
    "initial_fast_means": ("primary_covariance", PRIMITIVES[4]),
    "initial_slow_means": ("initial_fast_means", PRIMITIVES[4]),
    "initial_pressure_and_scalars": (
        "primary_covariance", "initial_fast_means", "initial_slow_means", PRIMITIVES[3], PRIMITIVES[4]),
    "initial_Haar_rows": (
        "primary_covariance", "initial_fast_means", "initial_slow_means", "initial_pressure_and_scalars"),
    "first_signed_target": ("initial_Haar_rows", "primary_root_and_complete_curl", PRIMITIVES[3]),
    "first_signed_covariance": ("first_signed_target", "primary_root_and_complete_curl"),
    "postsigned_pressure_and_rows": (
        "primary_covariance", "first_signed_covariance", "initial_fast_means", "initial_slow_means",
        "initial_pressure_and_scalars", PRIMITIVES[3], PRIMITIVES[4]),
    "Qcirc_Psi": ("postsigned_pressure_and_rows", "initial_fast_means", "initial_slow_means", PRIMITIVES[4]),
    "Qcirc_Phi": (
        "Qcirc_Psi", "first_signed_covariance", "postsigned_pressure_and_rows",
        "initial_fast_means", "initial_slow_means"),
}


def topological_order():
    done, order = set(PRIMITIVES), []
    remaining = dict(DEPENDENCIES)
    while remaining:
        ready = sorted(name for name, deps in remaining.items() if set(deps) <= done)
        if not ready:
            raise ValueError("The retained first-real dependency graph is cyclic or incomplete")
        for name in ready:
            order.append(name)
            done.add(name)
            del remaining[name]
    return tuple(order)


def fixed_order_induction(n=0):
    graph = request_module.request_graph(n)
    order = topological_order()
    J = graph["largest_immediate_order"]
    nodes = {name: dict(max_order=J, finite_positive_polynomial=True,
                        band_uniform_if_primitive_contract_holds=True)
             for name in order}
    return dict(
        output_order=n,
        request=graph,
        primitive_contract=PRIMITIVES,
        topological_order=order,
        nodes=nodes,
        requested_order_is_finite=isinstance(J, int) and J >= 0,
        numerical_materialization_required=False,
        desired_A1_output_used_as_input=False,
        primitive_contract_independently_certified=False,
    )


def rate_closure(kappa=Q(1, 100000)):
    k = Q(kappa)
    if not Q(0) < k <= Q(1, 100):
        raise ValueError("Use the retained positive small-kappa regime")
    B = Q(1)-k
    groups = dict(Qcirc_Psi=Q(1), gamma_square=2*B,
                  complete_Wzz=Q(1), Qcirc_complete_pressure=Q(1))
    phi = min(groups.values())
    dz_phi = phi+1
    p = request_module.request_graph(0, k)["integrations"]
    tail = dz_phi+p*k
    nu = Q(3)-6*k
    return dict(B=B, Phi_groups=groups, Qcirc_Phi=phi, Dz_Qcirc_Phi=dz_phi,
                integrations=p, A1_tail=tail, required=nu,
                strict_margin=tail-nu, strict=tail > nu)


def closure_statement(n=0):
    induction = fixed_order_induction(n)
    graph = induction["request"]
    return dict(
        conditional_Phi_inequality=(
            f"N_{graph['raw_Phi_order']}(Qcirc Phi_postsigned) "
            f"<= epsilon P_{graph['raw_Phi_order']}(S)"),
        conditional_DzPhi_inequality=(
            f"N_{graph['raw_Dz_Phi_order']}(Dz Qcirc Phi_postsigned) "
            f"<= epsilon^2 Q_{graph['raw_Dz_Phi_order']}(S)"),
        high_jet_is_not_an_additional_independent_premise=True,
        implication_requires_every_named_primitive_contract=True,
        actual_primitive_contracts_certified_here=False,
        full_pressure_memory_retained=True,
        whole_support_and_seam_closure_reused=True,
        A0_or_complete_residual_closed=False,
    )
