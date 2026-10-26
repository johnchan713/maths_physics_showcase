"""Conditional first-real radial A0 and compact physical-tail closure."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


phi = load("initial_phi_closure", HERE.parent/"core_initial_phi_high_jet_closure/closure.py")
mean = load("initial_a0_mean_graph", HERE.parent/"core_post_signed_mean_closure/budgets.py")


PRIMITIVES = phi.PRIMITIVES
DEPENDENCIES = dict(phi.DEPENDENCIES)
DEPENDENCIES.update({
    "first_real_fast_means": ("Qcirc_Phi", "postsigned_pressure_and_rows", PRIMITIVES[4]),
    "first_real_postfast_scalars": (
        "first_real_fast_means", "first_signed_covariance", "initial_pressure_and_scalars"),
    "first_real_slow_means": ("first_real_postfast_scalars", PRIMITIVES[4]),
    "final_mean_fields": (
        "initial_fast_means", "initial_slow_means", "first_real_fast_means", "first_real_slow_means"),
    "Qcirc_gr_final": ("final_mean_fields", "primary_covariance", "first_signed_covariance"),
    "Pnew_and_radial_compensator": ("Qcirc_gr_final", PRIMITIVES[3], PRIMITIVES[4]),
    "A0_Qcirc_gr_final": ("Qcirc_gr_final", "Pnew_and_radial_compensator", PRIMITIVES[4]),
})


def _order(n):
    if not isinstance(n, int) or isinstance(n, bool) or n < 0:
        raise ValueError("A nonnegative integer order is required")
    return n


def _kappa(kappa):
    k = Q(kappa)
    if not Q(0) < k <= Q(1, 100):
        raise ValueError("Use the retained positive small-kappa regime")
    return k


def topological_order():
    done, order = set(PRIMITIVES), []
    remaining = dict(DEPENDENCIES)
    while remaining:
        ready = sorted(name for name, deps in remaining.items() if set(deps) <= done)
        if not ready:
            raise ValueError("The retained first-real A0 dependency graph is cyclic or incomplete")
        for name in ready:
            order.append(name)
            done.add(name)
            del remaining[name]
    return tuple(order)


def normalized_request(n=0, kappa=Q(1, 100000)):
    n, k = _order(n), _kappa(kappa)
    B, nu = Q(1)-k, Q(3)-6*k
    p0 = (nu-B)//k+1
    p1 = (nu-Q(2))//k+1
    trace = mean.dependency_trace(n, p1, p0)
    return dict(
        output_order=n, kappa=k, source_order=B, target_order=nu,
        A0_integrations=p0, A1_integrations=p1,
        A0_source_order=trace["A0_source_order"],
        A1_DzPhi_source_order=trace["A1_source_order"],
        A1_Phi_source_order=trace["A1_source_order"]+1,
        achieved_order=B+p0*k, strict_margin=B+p0*k-nu,
        trace=trace, largest_same_prefix_order=max(trace["maxima"].values()),
    )


def physical_request(m=0, N=0, h=Q(1, 200), kappa=Q(1, 100000)):
    m, k, h, N = _order(m), _kappa(kappa), Q(h), Q(N)
    if not Q(0) < h < Q(1, 100) or N < 0:
        raise ValueError("Use fixed 0<h<1/100 and a nonnegative physical power")
    B = Q(1)-k
    Lambda = Q(3, 2)+2*h
    target = (N+Lambda+2*m)/h
    p0 = max(1, (target-B)//k+1)
    p1 = max(1, (target-Q(2))//k+1)
    trace = mean.dependency_trace(m, p1, p0)
    E0 = h*(B+p0*k)-Lambda-2*m
    E1 = h*(Q(2)+p1*k)-Lambda-2*m
    return dict(
        physical_jet=m, requested_q_power=N, h=h, kappa=k,
        normalization=Lambda, normalized_target=target,
        A0_integrations=p0, A1_integrations=p1,
        A0_source_order=trace["A0_source_order"],
        A1_DzPhi_source_order=trace["A1_source_order"],
        A1_Phi_source_order=trace["A1_source_order"]+1,
        A0_physical_exponent=E0, A1_physical_exponent=E1,
        A0_strict_margin=E0-N, A1_strict_margin=E1-N,
        trace=trace, largest_same_prefix_order=max(trace["maxima"].values()),
    )


def fixed_order_induction(*, physical=False, n=0, N=0, h=Q(1, 200)):
    request = physical_request(n, N, h) if physical else normalized_request(n)
    J = request["largest_same_prefix_order"]
    order = topological_order()
    nodes = {name: dict(max_order=J, finite_positive_polynomial=True,
                        band_uniform_if_primitive_contract_holds=True)
             for name in order}
    return dict(
        request=request, primitive_contract=PRIMITIVES, topological_order=order,
        nodes=nodes, requested_order_is_finite=isinstance(J, int) and J >= 0,
        numerical_materialization_required=False,
        desired_A0_output_used_as_input=False,
        primitive_contract_independently_certified=False,
    )


def radial_rate_closure(kappa=Q(1, 100000)):
    r = normalized_request(0, kappa)
    B, k = r["source_order"], r["kappa"]
    compensator = Q(3)-3*k
    groups = dict(
        time_on_radial_mean=B+1,
        radial_covariance_flux=Q(1)-k,
        axial_flux_at_least=Q(1),
        linear_centripetal_mean=B,
        quadratic_centripetal_and_covariance=min(2*B, Q(1)),
        vector_viscosity=B+2-2*k,
    )
    return dict(
        radial_source_groups=groups, Qcirc_gr_final=min(groups.values()),
        integrations=r["A0_integrations"],
        A0_tail=B+r["A0_integrations"]*k,
        radial_compensator=compensator, required=r["target_order"],
        full_radial_row=min(B+r["A0_integrations"]*k, compensator),
        strict_margin=r["strict_margin"], strict=r["strict_margin"] > 0,
    )


def physical_interface(m=0, N=0, h=Q(1, 200)):
    r = physical_request(m, N, h)
    return dict(
        compact_vector="(-A0 Qcirc gr_final, 0, +A1 Dz Qcirc Phi_postsigned)",
        conditional_bound=f"|Dphysical^{m} Tcompact| <= C_{{{m},{N}}} q^{N}",
        each_fixed_jet_and_power_only=True,
        A0_margin=r["A0_strict_margin"], A1_margin=r["A1_strict_margin"],
        same_accepted_band_and_profiles=True,
        Gaussian_and_base_exterior_contracts_reused_not_reproved=True,
        complete_regular_residual_improvement_proved=False,
        later_cutoff_interactions_and_diagonal_proved=False,
    )


def closure_statement():
    r = normalized_request()
    return dict(
        conditional_radial_source_inequality=(
            f"N_{r['A0_source_order']}(Qcirc gr_final) "
            f"<= epsilon^(1-kappa) P_{r['A0_source_order']}(S+K_all)"),
        conditional_A0_tail=(
            "N_0(A0 Qcirc gr_final) <= epsilon^(3-5 kappa) P(S+K_all)"),
        high_jet_is_not_an_additional_independent_premise=True,
        implication_requires_every_named_primitive_contract=True,
        actual_primitive_contracts_certified_here=False,
        complete_pressure_memory_and_unit_mass_compensation_retained=True,
        whole_support_endpoint_and_seam_closure_reused=True,
        compact_physical_tail_interface_closed_conditionally=True,
        complete_residual_or_blowup_closed=False,
    )
