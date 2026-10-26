"""First-sourced same-state compact-tail dependency and exponent ledger."""
from fractions import Fraction as Q


PRIMITIVES = (
    "analytic_background_and_Borel_coefficients",
    "same_band_phase_normal_root_and_value_cone_gaps",
    "canonical_forced_ODE_and_full_path_cutoffs",
    "fixed_rho_projectors_and_complete_pressure_maps",
    "fixed_compact_temporal_and_reserved_slow_inverse_maps",
    "whole_support_endpoint_and_seam_extensions",
)

DEPENDENCIES = {
    "first_real_state": PRIMITIVES,
    "first_real_regular_source": ("first_real_state", PRIMITIVES[1], PRIMITIVES[5]),
    "first_sourced_wave": ("first_real_regular_source", PRIMITIVES[1], PRIMITIVES[2]),
    "complete_covariance_change": ("first_real_state", "first_sourced_wave"),
    "postwave_pressure_moments": ("complete_covariance_change", PRIMITIVES[3]),
    "postwave_Qcirc_rows": ("first_real_state", "complete_covariance_change", "postwave_pressure_moments"),
    "postwave_A1_response": ("postwave_Qcirc_rows", PRIMITIVES[4], PRIMITIVES[5]),
    "fresh_fast_means": ("postwave_Qcirc_rows", "postwave_A1_response", PRIMITIVES[4]),
    "postfast_pressure_scalars": ("fresh_fast_means", "complete_covariance_change", PRIMITIVES[3]),
    "fresh_slow_means": ("postfast_pressure_scalars", PRIMITIVES[4]),
    "final_first_sourced_state": ("first_real_state", "first_sourced_wave", "fresh_fast_means", "fresh_slow_means"),
    "Qcirc_gr_S": ("final_first_sourced_state", "complete_covariance_change", PRIMITIVES[3]),
    "A0_Qcirc_gr_S": ("Qcirc_gr_S", PRIMITIVES[4], PRIMITIVES[5]),
}


def _order(value, name):
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _parameters(h=Q(1, 200), kappa=Q(1, 100000)):
    h, k = Q(h), Q(kappa)
    if not Q(0) < h < Q(1, 100):
        raise ValueError("use the retained 0<h<1/100 regime")
    if not Q(0) < k < Q(1, 4):
        raise ValueError("use the retained 0<kappa<1/4 regime")
    return h, k


def topological_order():
    done, out, remaining = set(PRIMITIVES), [], dict(DEPENDENCIES)
    while remaining:
        ready = sorted(name for name, deps in remaining.items() if set(deps) <= done)
        if not ready:
            raise ValueError("the first-sourced compact dependency graph is cyclic or incomplete")
        for name in ready:
            out.append(name)
            done.add(name)
            del remaining[name]
    return tuple(out)


def normalized_request(n=0, kappa=Q(1, 100000)):
    n = _order(n, "n")
    _, k = _parameters(Q(1, 200), kappa)
    B, C, D = Q(1)-k, Q(3, 2)-k, Q(3, 2)-2*k
    p0 = max(1, (C-B)//k+1)
    p1 = max(1, (C-D)//k+1)
    return request_ledger(n, 0, p0, p1, B, D, C)


def physical_request(m=0, N=0, h=Q(1, 200), kappa=Q(1, 100000)):
    m, N = _order(m, "m"), _order(N, "N")
    h, k = _parameters(h, kappa)
    B, D = Q(1)-k, Q(3, 2)-2*k
    Lambda = Q(3, 2)+2*h
    target = Q(N+2*m)+Lambda
    target /= h
    p0 = max(1, (target-B)//k+1)
    p1 = max(1, (target-D)//k+1)
    out = request_ledger(m, N, p0, p1, B, D, target)
    out.update(h=h, normalization=Lambda,
               A0_physical_exponent=h*(B+p0*k)-Lambda-2*m,
               A1_physical_exponent=h*(D+p1*k)-Lambda-2*m)
    return out


def request_ledger(base_order, requested_power, p0, p1, B, D, target):
    J0, J1 = base_order+2*p0+3, base_order+2*p1+3
    final_field = J0+3
    direct_A1_first_real = J1+200004
    radial_nested_A1_first_real = final_field+200013
    radial_sourced_wave_first_real = final_field+28
    largest = max(direct_A1_first_real, radial_nested_A1_first_real,
                  radial_sourced_wave_first_real)
    return dict(
        output_order=base_order, requested_q_power=requested_power,
        source_orders=dict(A0=B, A1=D), target_order=target,
        integrations=dict(A0=p0, A1=p1), raw_sources=dict(A0=J0, A1=J1),
        final_field_order=final_field,
        direct_A1_first_real_ceiling=direct_A1_first_real,
        radial_nested_A1_first_real_ceiling=radial_nested_A1_first_real,
        radial_sourced_wave_first_real_ceiling=radial_sourced_wave_first_real,
        largest_first_real_ceiling=largest,
        radial_branch_dominates=largest == radial_nested_A1_first_real,
    )


def fixed_order_induction(*, physical=False, n=0, N=0,
                          h=Q(1, 200), kappa=Q(1, 100000)):
    request = (physical_request(n, N, h, kappa) if physical
               else normalized_request(n, kappa))
    order = topological_order()
    return dict(
        request=request, primitives=PRIMITIVES, topological_order=order,
        requested_ceiling_is_finite=(
            isinstance(request["largest_first_real_ceiling"], int)
            and request["largest_first_real_ceiling"] >= 0),
        desired_A0_output_used_earlier=False,
        postwave_A1_uses_pre_fast_row=True,
        every_node_is_finite_positive_polynomial_given_primitives=True,
        numerical_high_jet_materialization_required=False,
        primitive_contracts_independently_certified=False,
    )


def residual_interface(m=0, N=0, h=Q(1, 200), kappa=Q(1, 100000), eta=1):
    r = physical_request(m, N, h, kappa)
    h, k, eta = Q(h), Q(kappa), Q(eta)
    if eta <= 0:
        raise ValueError("a positive polynomial reserve is required")
    D, Lambda = Q(3, 2)-2*k, r["normalization"]
    regular_exponent = h*D-Lambda-2*m
    return dict(
        compact_vector="(-A0 Qcirc gr_S,0,+A1 Qcirc Ez_postwave)",
        conditional_compact_bound=f"|Dphysical^{m} Tcompact_S| <= C_{{{m},{N}}} q^{N}",
        compact_flat_for_each_fixed_m_N=(
            r["A0_physical_exponent"] > N and r["A1_physical_exponent"] > N),
        regular_order=D, regular_physical_exponent=regular_exponent,
        sufficient_later_order_threshold=(Q(N)+Lambda+2*m+eta)/h,
        current_regular_decay_certified=regular_exponent-eta >= N,
        Gaussian_and_base_exterior_flatness_reused_conditionally=True,
        complete_residual_decay_proved=False,
    )
