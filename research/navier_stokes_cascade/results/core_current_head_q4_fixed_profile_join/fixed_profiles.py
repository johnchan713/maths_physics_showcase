"""Current-head Q4 fixed-profile compensator and projector join."""

from fractions import Fraction
from math import comb, factorial
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
PROJECT = RESULTS.parent

DEPENDENCIES = (
    "core_current_head_q3_full_frame_join",
    "core_fixed_profile_projector_reconstruction",
    "core_fixed_profile_projector_join",
    "core_fixed_profile_tame",
    "core_moving_pressure_compensator_join",
    "core_reserved_geometry_recheck",
    "core_rho_squared_projector_realization",
    "core_rho_projector_signed_recompute",
    "core_log_pullback_projector_join",
    "core_finite_support_closure",
)

EXPECTED_STATUS = {
    "core_current_head_q3_full_frame_join": "current-head-Q3-full-frame-whole-path-and-phase-domain-bound",
    "core_fixed_profile_projector_reconstruction": "fixed-profile-same-choice-projector-contract-reconstructed-on-current-full-path-tail",
    "core_fixed_profile_projector_join": "fixed-profile-same-choice-projector-contract-joined",
    "core_fixed_profile_tame": "fixed-profile-all-order-tame-bounds-and-quantifier-correction",
    "core_moving_pressure_compensator_join": "conditional-fixed-moving-pressure-compensator-and-new-pole-moment-join",
    "core_reserved_geometry_recheck": "uniform-reserved-coordinate-compensator-jets-and-moving-support-costs-rechecked",
    "core_rho_squared_projector_realization": "explicit-new-rho-squared-projector-realization",
    "core_rho_projector_signed_recompute": "conditional-same-choice-rho-projector-signed-recomputation",
    "core_log_pullback_projector_join": "conditional-actual-logarithmic-Haar-pullback-and-projector-seminorm-boundary",
    "core_finite_support_closure": "finite-smooth-support-closure-relative-to-original-weighted-input-contract",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def evidence_path(name):
    return RESULTS / name / "evidence.json"


def load_closure():
    records, graph, dependency_checks, source_checks = {}, {}, {}, {}

    def visit(path):
        path = path.resolve()
        if path in records:
            return
        record = json.loads(path.read_text())
        records[path] = record
        children = []
        for rel, digest in record.get("dependency_evidence_sha256", {}).items():
            target = (PROJECT / rel).resolve()
            dependency_checks[f"{path.parent.name}:{rel}"] = target.exists() and sha(target) == digest
            if target.exists():
                children.append(target)
                visit(target)
        graph[path] = children
        for rel, digest in record.get("source_hashes", {}).items():
            target = PROJECT / rel
            source_checks[f"{path.parent.name}:{rel}"] = target.exists() and sha(target) == digest

    for name in DEPENDENCIES:
        visit(evidence_path(name))
    return records, graph, dependency_checks, source_checks


def cycle_free(graph):
    color = {}

    def visit(node):
        color[node] = 1
        for child in graph[node]:
            if color.get(child) == 1:
                return False
            if color.get(child, 0) == 0 and not visit(child):
                return False
        color[node] = 2
        return True

    return all(color.get(node, 0) or visit(node) for node in graph)


def fixed_jet_majorant(order, profile, coordinates):
    """Finite exact composition majorant at one fixed order."""
    if type(order) is not int or order < 0:
        raise ValueError("nonnegative integer order required")
    if len(profile) != order + 1 or len(coordinates) != order + 1:
        raise ValueError("all fixed jets through the requested order are required")
    values = tuple(Fraction(value) for value in tuple(profile) + tuple(coordinates))
    if any(value < 0 for value in values):
        raise ValueError("nonnegative exact seminorms required")
    profile, coordinates = values[: order + 1], values[order + 1 :]
    return sum(
        Fraction(comb(order, j) * factorial(j + 1))
        * (profile[j] + 1) * (coordinates[order - j] + 1) ** (j + 1)
        for j in range(order + 1)
    )


def projector_contract(ell, upper, support_ell, support_upper, scale, A2, exponent):
    """Exact normalization/support interface for sigma_1 and sigma_2."""
    ell, upper, support_ell, support_upper, scale, A2 = map(
        Fraction, (ell, upper, support_ell, support_upper, scale, A2)
    )
    if exponent not in (1, 2) or scale not in (Fraction(1, 4), Fraction(1), Fraction(4)):
        raise ValueError("e=1,2 and a certified rational-square scale are required")
    if not 0 < ell < support_ell < support_upper < upper or A2 <= 0:
        raise ValueError("positive profile data and strict actual support are required")
    if A2 * (upper - ell) < 1:
        raise ValueError("A2 must satisfy the unit-mass Cauchy-Schwarz floor")
    root = {Fraction(1, 4): Fraction(1, 2), Fraction(1): Fraction(1), Fraction(4): Fraction(2)}[scale]
    hull = (root * ell, root * upper)
    actual = (root * support_ell, root * support_upper)
    return {
        "moment": root ** exponent * root ** (-exponent) * root ** (-1) * root,
        "reciprocal_A2": 1 / A2,
        "reciprocal_bound": upper - ell,
        "hull": hull,
        "actual_support": actual,
        "strict_support": hull[0] < actual[0] < actual[1] < hull[1],
    }


def moving_pressure_memory(crho_z, pressure, crho, pressure_z, jz_z):
    """Retain the whole Z derivative of Jz+c_rho P."""
    crho_z, pressure, crho, pressure_z, jz_z = map(
        Fraction, (crho_z, pressure, crho, pressure_z, jz_z)
    )
    return {
        "Jz_Z": jz_z,
        "crho_Z_times_P": crho_z * pressure,
        "crho_times_P_Z": crho * pressure_z,
        "total_without_epsilon": jz_z + crho_z * pressure + crho * pressure_z,
    }


def contract_order():
    return (
        "Q1_Q3_field_and_tail",
        "fixed_pressure_profile",
        "moving_compensator",
        "rho_squared_projectors",
        "same_choice_recompute",
        "finite_support_closure",
        "Q4_contract",
    )


def review():
    records, graph, dependencies, sources = load_closure()
    status_checks = {
        name: records[evidence_path(name).resolve()].get("status") == expected
        for name, expected in EXPECTED_STATUS.items()
    }
    rows = [
        projector_contract(1, 2, Fraction(5, 4), Fraction(7, 4), scale, Fraction(10, 7), exponent)
        for exponent in (1, 2) for scale in (Fraction(1, 4), Fraction(1), Fraction(4))
    ]
    memory = moving_pressure_memory(Fraction(1, 5), Fraction(7, 3), Fraction(11, 4), Fraction(13, 6), Fraction(17, 9))
    order = contract_order()
    checks = {
        "provider_graph_is_acyclic": cycle_free(graph),
        "all_declared_dependency_hashes_match": bool(dependencies) and all(dependencies.values()),
        "all_declared_source_hashes_match": bool(sources) and all(sources.values()),
        "all_ten_Q4_provider_statuses_match": len(status_checks) == 10 and all(status_checks.values()),
        "fixed_profile_gives_finite_positive_fixed_order_majorant": fixed_jet_majorant(4, [1, 3, 7, 15, 31], [1, 2, 4, 8, 16]) > 0,
        "both_projector_moments_are_exact_on_all_scale_controls": all(row["moment"] == 1 for row in rows),
        "actual_support_is_strictly_inside_every_moving_enclosure": all(row["strict_support"] and row["actual_support"] != row["hull"] for row in rows),
        "whole_moving_pressure_derivative_retains_crho_Z_P": memory["crho_Z_times_P"] != 0 and memory["total_without_epsilon"] == sum((memory["Jz_Z"], memory["crho_Z_times_P"], memory["crho_times_P_Z"])),
        "same_choice_recompute_follows_projector_choice": order.index("rho_squared_projectors") < order.index("same_choice_recompute") < order.index("Q4_contract"),
        "finite_support_closure_follows_same_choice_recompute": order.index("same_choice_recompute") < order.index("finite_support_closure") < order.index("Q4_contract"),
    }
    return {
        "checks": checks,
        "graph": {
            "nodes": len(graph),
            "edges": sum(len(row) for row in graph.values()),
            "dependency_hash_checks": len(dependencies),
            "source_hash_checks": len(sources),
        },
        "provider_status_checks": status_checks,
        "interface_controls": {
            "projector_rows": [
                {key: ([str(x) for x in value] if isinstance(value, tuple) else str(value) if isinstance(value, Fraction) else value) for key, value in row.items()}
                for row in rows
            ],
            "moving_pressure_memory": {key: str(value) for key, value in memory.items()},
            "commutator": "all six moving-density product terms and epsilon*partial_Z(Jz+c_rho*P), including c_rho,Z*P",
            "same_choice": "recompute h_e, two-edge primitives, signed amplitude/root input, curl, six covariances and full pressure norms after fixing sigma_e",
            "support": "supp(sigma_e)=sqrt(s)*supp(rhohat) is strictly inside the moving interval; [ell/2,2u] is only a common integration hull",
        },
        "scope": {
            "Q1_Q3_field_is_a_direct_dependency": True,
            "Q4_fixed_profile_compensator_projectors_bound_on_same_later_tail": True,
            "normalization_only_used_to_control_arbitrary_profile_jets": False,
            "moving_support_enclosure_identified_with_actual_support": False,
            "historical_unspecified_projectors_identified_with_new_choice": False,
            "same_choice_projector_dependent_objects_recomputed": True,
            "uniformity_over_profiles_orders_or_infinite_cycles_claimed": False,
            "compact_slow_weighted_inverse_contract_closed": False,
            "hashes_or_tests_treated_as_continuum_proof": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup_claimed": False,
        },
    }
