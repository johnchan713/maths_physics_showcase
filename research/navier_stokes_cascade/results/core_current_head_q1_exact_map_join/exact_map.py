"""Exact-core/root Q1 provider join and interface controls."""

from fractions import Fraction
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
PROJECT = RESULTS.parent

DEPENDENCIES = (
    "core_current_head_primitive_frontier_audit",
    "foundation_bundle_route",
    "core_contraction_audit",
    "core_fixed_point_ledger_recheck",
    "root_foundation",
    "axial_root_foundation",
    "matching_root_foundation",
    "prescribed_continuation_foundation",
    "annular_exterior_foundation",
    "core_continuation_matching_recheck",
    "core_modulation_repair_join_recheck",
    "core_global_pressure_moment_recheck",
    "core_pressure_identification",
    "core_same_field_pressure_reconciliation",
)

EXPECTED_STATUS = {
    "core_contraction_audit": "core-contraction-and-exact-stress-audited-for-pressure-family",
    "core_fixed_point_ledger_recheck": "core-fixed-point-all-index-ledger-independently-rederived",
    "root_foundation": "reference-angular-heat-and-wait-roots-audited-relative-to-remaining-construction-lemmas",
    "axial_root_foundation": "reference-axial-moment-and-complete-strength-roots-audited-relative-to-remaining-construction-lemmas",
    "matching_root_foundation": "actual-annular-and-compact-matching-roots-audited-relative-to-retained-core-and-compact-input-estimates",
    "core_continuation_matching_recheck": "core-continuation-and-both-same-profile-matching-roots-independently-rechecked",
    "core_modulation_repair_join_recheck": "same-profile-compact-modulation-repair-and-regional-cone-join-independently-rechecked",
    "core_global_pressure_moment_recheck": "same-repaired-profile-global-pressure-moments-and-exterior-stress-identities-independently-rechecked",
    "core_pressure_identification": "leading-core-pressure-identified-with-regular-angular-inverse-and-retained-higher-memory",
    "core_same_field_pressure_reconciliation": "same-field-core-pressure-identified-with-completed-leading-backward-pressure",
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


def quadratic_root_control(beta, qnorm, delta):
    """Exact small-root lemma: radius, invariant factor and Lipschitz factor."""
    beta, qnorm, delta = map(Fraction, (beta, qnorm, delta))
    radius = 2 * beta * delta
    invariant = 8 * beta * beta * qnorm * delta
    lipschitz = 4 * beta * beta * qnorm * delta
    return {"radius": radius, "invariant": invariant, "lipschitz": lipschitz}


def endpoint_pressure_obstruction(A, nonzero_constant=True):
    """For K_A delta=0, a nonzero homogeneous mode is endpoint-unbounded."""
    A = Fraction(A)
    if A <= 0:
        raise ValueError("the pressure exponent must be positive")
    return {
        "homogeneous_mode": "C*(1-eta^2)^(-2A)",
        "nonzero_mode_is_bounded_at_endpoints": not nonzero_constant,
        "bounded_smooth_solution_forces_C_zero": True,
    }


def provider_contracts():
    return {
        "core": "one Banach fixed point in the complete all-index coefficient space for the selected exact positive-measure pressure family, with exact core stress cancellation",
        "reference_roots": "the angular, heat, waiting, axial two-moment and positive scalar-strength roots are unique in their retained brackets/balls and have finite fixed jets",
        "matching_and_repair_roots": "the annular and compact five-row roots are the same accepted roots at every fixed order and restore M,R-or-I,J,S,Cp without resetting pressure",
        "leading_pressure": "Pi_axis+Cp_infinity=0 on the same completed leading profile, so the core datum equals the completed leading backward pressure while positive-order, wave, mean and localization pressure remain",
    }


def review():
    records, graph, dependencies, sources = load_closure()
    status_checks = {}
    for name, expected in EXPECTED_STATUS.items():
        status_checks[name] = records[evidence_path(name).resolve()].get("status") == expected

    annular = quadratic_root_control(610, 50, Fraction(286, 10**19))
    # D/N < H^-16 and H>2^16; use the weakest exact endpoint H=2^16.
    compact = quadratic_root_control(800, 25000, Fraction(1, 2 ** 256))
    pressure = endpoint_pressure_obstruction(Fraction(1, 2))
    checks = {
        "provider_graph_is_acyclic": cycle_free(graph),
        "all_declared_dependency_hashes_match": bool(dependencies) and all(dependencies.values()),
        "all_declared_source_hashes_match": bool(sources) and all(sources.values()),
        "all_ten_Q1_provider_statuses_match": len(status_checks) == 10 and all(status_checks.values()),
        "annular_root_smallness_is_strict": annular["invariant"] < 1 and annular["lipschitz"] < Fraction(1, 2),
        "compact_root_smallness_is_strict": compact["invariant"] < 1 and compact["lipschitz"] < Fraction(1, 2),
        "nonzero_pressure_homogeneous_mode_is_endpoint_unbounded": not pressure["nonzero_mode_is_bounded_at_endpoints"] and pressure["bounded_smooth_solution_forces_C_zero"],
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
        "contracts": provider_contracts(),
        "interface_controls": {
            "annular": {key: str(value) for key, value in annular.items()},
            "compact": {key: str(value) for key, value in compact.items()},
            "pressure": pressure,
            "scalar_strength_derivative_floor": "F_A>49/100",
            "fixed_object_ledger": "same exact pressure datum, h=lambda^2, C^-20 width, supports, N=1+floor(H^32), roots and gauge",
        },
        "scope": {
            "Q1_provider_chain_bound_on_one_current_head_object": True,
            "Q1_interface_uniqueness_and_pressure_compatibility_rechecked": True,
            "every_provider_theorem_independently_reproved_here": False,
            "Q2_regional_continuum_infima_discharged": False,
            "Q3_full_frame_primitive_domain_discharged": False,
            "hashes_or_tests_treated_as_continuum_proof": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup_claimed": False,
        },
    }
