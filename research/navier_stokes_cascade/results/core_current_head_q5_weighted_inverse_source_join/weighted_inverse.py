"""Current-head Q5 weighted inverse and log-source join."""

from fractions import Fraction
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
PROJECT = RESULTS.parent

DEPENDENCIES = (
    "core_current_head_q4_fixed_profile_join",
    "core_compact_slow_inverse_reconstruction",
    "core_compact_slow_inverse_domain_join",
    "core_direct_compact_slow_provider_reconstruction",
    "core_current_log_source_reconciliation",
    "core_log_renormalized_actual_source_induction",
    "core_fast_inverse_budget",
    "core_slow_inverse_budget",
    "core_initial_haar_row_budget",
    "core_native_haar_primitive_join",
)

EXPECTED_STATUS = {
    "core_current_head_q4_fixed_profile_join": "current-head-Q4-fixed-profile-compensator-and-same-choice-projectors-bound",
    "core_compact_slow_inverse_reconstruction": "compact-slow-weighted-inverse-contract-reconstructed-on-current-same-choice-tail",
    "core_compact_slow_inverse_domain_join": "compact-slow-weighted-domain-and-support-contract-joined",
    "core_direct_compact_slow_provider_reconstruction": "direct-compact-slow-operators-proved-with-actual-source-supremum-frontier",
    "core_current_log_source_reconciliation": "current-five-contract-log-renormalized-actual-source-induction-reconciled-conditionally",
    "core_log_renormalized_actual_source_induction": "log-renormalized-actual-source-tail-induction-proved-with-residual-replay-frontier",
    "core_fast_inverse_budget": "explicit-fast-inverse-actual-axial-repair-and-current-residual-jet-budget",
    "core_slow_inverse_budget": "explicit-fixed-profile-slow-inverse-and-coupled-scalar-substitution",
    "core_initial_haar_row_budget": "initial-complete-weighted-Haar-row-budget-and-sharper-average-rate",
    "core_native_haar_primitive_join": "conditional-native-Haar-pressure-and-target-primitive-join",
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


def two_edge_constant(inner, outer, interior):
    """Symbolic three-region compact-kernel coefficient interface."""
    inner, outer, interior = map(Fraction, (inner, outer, interior))
    if min(inner, outer, interior) <= 0:
        raise ValueError("all two-edge branches need positive finite coefficients")
    return {"inner": inner, "outer": outer, "interior": interior, "sum": inner + outer + interior}


def request_graph(output_order, raw_A0, raw_A1):
    """Finite derivative requests, including the log route's extra integration."""
    if any(type(value) is not int or value < 0 for value in (output_order, raw_A0, raw_A1)):
        raise ValueError("nonnegative integer requests required")
    return {
        "fast_tangential_or_axial": output_order + 4,
        "fast_radial_potential": output_order + 5,
        "prefast_source": output_order + 8,
        "whole_axial_pressure_moment": output_order + 9,
        "A0_with_extra_integration": raw_A0 + 2,
        "A1_with_extra_integration": raw_A1 + 2,
    }


def log_reserve(h, kappa, degree):
    """Exact exponent allocation for q^(h*kappa)(1+|log eps|)^P."""
    h, kappa = map(Fraction, (h, kappa))
    if h <= 0 or h > 1 or kappa <= 0 or type(degree) is not int or degree < 0:
        raise ValueError("0<h<=1, kappa>0 and finite nonnegative log degree required")
    gain = h * kappa
    spent = gain / 2
    return {
        "gain": gain,
        "spent_on_logs": spent,
        "retained": gain - spent,
        "constant_exponent": Fraction(degree * degree, 1) / (2 * gain),
    }


def pure_power_equality(kappa):
    """Old positive-majorant equality; it is not the log-domain frontier."""
    kappa = Fraction(kappa)
    if not 0 < kappa < Fraction(1, 10):
        raise ValueError("small positive kappa required")
    B = 1 - kappa
    D = 2 - 4 * kappa
    nu = 3 - 6 * kappa
    viscosity = D + 1 - 2 * kappa
    return {
        "B": B,
        "D": D,
        "nu": nu,
        "radial_viscosity": viscosity,
        "equality": viscosity == nu,
        "source_margin": Fraction(0),
        "viscosity_margin": viscosity - nu,
    }


def fixed_stage_domain(stage, jet, target, log_degree):
    """Freeze the correct order of finite-stage source quantifiers."""
    if any(type(value) is not int or value < 0 for value in (stage, jet, target, log_degree)):
        raise ValueError("fixed finite nonnegative indices are required")
    return {
        "stage_fixed_first": True,
        "jet_fixed_first": True,
        "target_fixed_first": True,
        "finite_log_degree": log_degree,
        "stage_uniform_constant_required": False,
        "derivative_uniform_constant_required": False,
    }


def review():
    records, graph, dependencies, sources = load_closure()
    status_checks = {
        name: records[evidence_path(name).resolve()].get("status") == expected
        for name, expected in EXPECTED_STATUS.items()
    }
    kernel = two_edge_constant(Fraction(4, 49), Fraction(401, 100), Fraction(1, 4900))
    requests = request_graph(0, 399995, 200001)
    reserve = log_reserve(Fraction(1, 200), Fraction(1, 100000), 37)
    equality = pure_power_equality(Fraction(1, 100000))
    domain = fixed_stage_domain(7, 13, 19, 37)
    checks = {
        "provider_graph_is_acyclic": cycle_free(graph),
        "all_declared_dependency_hashes_match": bool(dependencies) and all(dependencies.values()),
        "all_declared_source_hashes_match": bool(sources) and all(sources.values()),
        "all_ten_Q5_provider_statuses_match": len(status_checks) == 10 and all(status_checks.values()),
        "both_edge_and_interior_kernel_branches_are_positive": min(kernel["inner"], kernel["outer"], kernel["interior"]) > 0,
        "fast_and_pressure_derivative_losses_are_distinct": requests["fast_tangential_or_axial"] == 4 and requests["fast_radial_potential"] == 5 and requests["whole_axial_pressure_moment"] == 9,
        "extra_compact_integration_costs_two_A0_and_A1_derivatives": requests["A0_with_extra_integration"] == 399997 and requests["A1_with_extra_integration"] == 200003,
        "half_log_reserve_is_strictly_positive": reserve["gain"] == reserve["spent_on_logs"] + reserve["retained"] and reserve["retained"] > 0,
        "old_pure_power_radial_viscosity_has_exact_zero_margin": equality["equality"] and equality["source_margin"] == 0 and equality["viscosity_margin"] == 0,
        "fixed_stage_quantifiers_require_no_uniform_order_constant": domain["stage_fixed_first"] and domain["finite_log_degree"] == 37 and not domain["stage_uniform_constant_required"] and not domain["derivative_uniform_constant_required"],
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
            "two_edge_kernel": {key: str(value) for key, value in kernel.items()},
            "request_graph": requests,
            "log_reserve": {key: str(value) for key, value in reserve.items()},
            "pure_power_equality": {key: str(value) if isinstance(value, Fraction) else value for key, value in equality.items()},
            "fixed_stage_domain": domain,
            "weighted_domain": "epsilon^b*(1+|log epsilon|)^P*zeta*delta^-p with finite P after fixing stage, jet and target",
            "fast_map": "retain y+Ma*v_r, u=-c^-1 Nt^-1 Ftheta, gd=-c^-1 Nt^-1 Fz, gamma=gd-A1*gd, beta=-epsilon*partial_Z T1*gd",
            "slow_map": "fixed five-profile response on the complete reserved potential hull with finite fixed lambda^-1 conditioning",
            "pressure": "reconstruct after fast and slow maps; retain Jz+c_rho*P, c_rho,Z*P, axial velocity, six covariances and Pi1/Pi2",
        },
        "scope": {
            "Q1_Q4_object_is_a_direct_dependency": True,
            "Q5_weighted_inverse_and_log_source_contract_bound_on_same_later_tail": True,
            "raw_pure_power_supremum_promoted_as_exact_current_source_domain": False,
            "old_zero_margin_equality_erased": False,
            "physical_characteristic_shift_removed": False,
            "actual_A1_axial_velocity_removed": False,
            "finite_stage_log_source_induction_instantiated": True,
            "six_group_complete_residual_replayed": False,
            "hashes_or_tests_treated_as_continuum_proof": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup_claimed": False,
        },
    }
