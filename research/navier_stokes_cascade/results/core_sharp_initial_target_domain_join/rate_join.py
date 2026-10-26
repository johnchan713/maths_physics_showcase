"""Exact rate controls for the sharper initial target/domain join."""
from fractions import Fraction as Q


NODES = (
    "weighted_initial_haar_rows",
    "same_choice_projection",
    "two_edge_primitive",
    "original_root_signed_map",
    "first_real_covariance",
    "compact_slow_domain",
    "strict_first_real_final_row",
)
EDGES = (
    ("weighted_initial_haar_rows", "same_choice_projection"),
    ("same_choice_projection", "two_edge_primitive"),
    ("two_edge_primitive", "original_root_signed_map"),
    ("original_root_signed_map", "first_real_covariance"),
    ("first_real_covariance", "strict_first_real_final_row"),
    ("compact_slow_domain", "strict_first_real_final_row"),
)


def topological_order(nodes=NODES, edges=EDGES):
    incoming = {node: 0 for node in nodes}
    outgoing = {node: [] for node in nodes}
    for left, right in edges:
        if left not in incoming or right not in incoming:
            raise ValueError("Unknown dependency node")
        incoming[right] += 1
        outgoing[left].append(right)
    ready = [node for node in nodes if incoming[node] == 0]
    result = []
    while ready:
        node = ready.pop(0)
        result.append(node)
        for right in outgoing[node]:
            incoming[right] -= 1
            if incoming[right] == 0:
                ready.append(right)
    if len(result) != len(nodes):
        raise ValueError("Dependency cycle")
    return result


def rates(kappa=Q(1, 100000)):
    k = Q(kappa)
    if not 0 < k < Q(1, 4):
        raise ValueError("Need 0<kappa<1/4")
    alpha = 2 - 3 * k
    stale_source = alpha - k
    nu = stale_source + 1 - 2 * k
    beta = 2 - 2 * k
    projected_source = min(beta, Q(3))
    pointwise_response = projected_source - k
    viscosity = pointwise_response + 1 - 2 * k
    return {
        "kappa": k,
        "alpha": alpha,
        "stale_D": stale_source,
        "nu": nu,
        "beta": beta,
        "projected_source": projected_source,
        "pointwise_response": pointwise_response,
        "radial_viscosity": viscosity,
        "strict_margin": viscosity - nu,
    }


def covariance_rates(kappa=Q(1, 100000)):
    row = rates(kappa)
    k, beta = row["kappa"], row["beta"]
    signed_principal = beta - Q(1, 2)
    signed_curl = beta - k
    return {
        "signed_principal": signed_principal,
        "signed_curl": signed_curl,
        "old_new_cross": beta + 1 - 2 * k,
        "signed_square": 2 * beta - 1,
        "curl_square": 2 * beta - 2 * k,
        "real_pair_remainder": min(beta + 1 - 2 * k, 2 * beta - 1, 2 * beta - 2 * k),
    }


def final_groups(kappa=Q(1, 100000)):
    row = rates(kappa)
    k, alpha, source = row["kappa"], row["alpha"], row["pointwise_response"]
    source_groups = {
        "original_fast_difference": alpha,
        "original_slow_difference": row["beta"],
        "signed_covariance_response": source,
        "corrected_original_A1_complement": Q(2),
    }
    fast_groups = {
        "slow_time": source + 1,
        "base_radial_flux": source + 1 - k,
        "radial_viscosity": source + 1 - 2 * k,
        "old_mean_transport": (1 - k) + source + 1 - k,
        "changed_axial_pressure": source + 1,
        "new_mean_square": 2 * source + 1 - k,
    }
    final = {
        "postsigned_Haar": covariance_rates(kappa)["real_pair_remainder"],
        "fresh_fast_difference": min(fast_groups.values()),
        "retained_A1_raw_bound": 3 - 5 * k,
        "sharp_slow_difference": 3 - 3 * k,
    }
    return {
        "source_groups": source_groups,
        "complete_fast_source": min(source_groups.values()),
        "fast_groups": fast_groups,
        "final_groups": final,
        "conditional_final_order": min(final.values()),
        "strict_margin": min(final.values()) - row["nu"],
    }


def review():
    order = topological_order()
    row = rates()
    covariance = covariance_rates()
    final = final_groups()
    checks = {
        "dependency_route_is_acyclic": set(order) == set(NODES),
        "Haar_rate_exceeds_stale_alpha_by_kappa": row["beta"] - row["alpha"] == row["kappa"],
        "moment_bumps_do_not_limit_beta": row["projected_source"] == row["beta"],
        "linear_maps_preserve_projected_epsilon_order": row["pointwise_response"] == row["alpha"],
        "signed_principal_rate": covariance["signed_principal"] == Q(3, 2) - 2 * row["kappa"],
        "real_pair_remainder_rate": covariance["real_pair_remainder"] == 3 - 4 * row["kappa"],
        "complete_fast_source_rate": final["complete_fast_source"] == row["alpha"],
        "radial_viscosity_rate_is_improved": final["fast_groups"]["radial_viscosity"] == 3 - 5 * row["kappa"],
        "conditional_final_order": final["conditional_final_order"] == 3 - 5 * row["kappa"],
        "strict_margin_is_kappa": final["strict_margin"] == row["kappa"] > 0,
        "stale_zero_margin_is_corrected": row["radial_viscosity"] - row["nu"] == row["kappa"],
    }
    return {
        "checks": checks,
        "rates": {key: str(value) for key, value in row.items()},
        "covariance_rates": {key: str(value) for key, value in covariance.items()},
        "final_groups": {
            key: ({name: str(value) for name, value in values.items()} if isinstance(values, dict) else str(values))
            for key, values in final.items()
        },
        "scope": "exact rational rate propagation; actual same-state high-jet and support enclosures remain separate continuum obligations",
    }
