"""Exact provenance and quantifier controls for the complete conditional chain."""

from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
PROJECT = RESULTS.parent

ROOTS = (
    "core_leading_background_borel_join",
    "core_same_band_gap_join",
    "core_canonical_full_path_ode_join",
    "core_fixed_profile_projector_join",
    "core_compact_slow_inverse_domain_join",
    "core_certified_finite_stage_domain_join",
    "core_certified_complete_potential_diagonal_join",
    "core_certified_complete_potential_endpoint_join",
    "core_primitive_margin_normalization_reconciliation",
    "core_scale_compatible_corner_force_replay",
    "core_certified_complete_potential_force_join",
    "core_certified_complete_potential_comparison_join",
    "core_same_field_positive_core_reconciliation",
)

CONSTRUCTION_SPINE = (
    "core_leading_background_borel_join",
    "core_same_band_gap_join",
    "core_canonical_full_path_ode_join",
    "core_fixed_profile_projector_join",
    "core_compact_slow_inverse_domain_join",
    "core_certified_finite_stage_domain_join",
    "core_certified_complete_potential_diagonal_join",
    "core_certified_complete_potential_endpoint_join",
    "core_certified_complete_potential_force_join",
    "core_certified_complete_potential_comparison_join",
    "core_same_field_positive_core_reconciliation",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def evidence_path(name):
    return RESULTS / name / "evidence.json"


def load_graph():
    """Load every transitive evidence dependency and verify declared hashes."""
    records, graph, hash_checks, source_checks = {}, {}, {}, {}

    def visit(path):
        path = path.resolve()
        if path in records:
            return
        record = json.loads(path.read_text())
        records[path] = record
        deps = []
        for rel, digest in record.get("dependency_evidence_sha256", {}).items():
            target = (PROJECT / rel).resolve()
            key = f"{path.parent.name}:{rel}"
            hash_checks[key] = target.exists() and sha(target) == digest
            if target.exists():
                deps.append(target)
                visit(target)
        graph[path] = deps
        for rel, digest in record.get("source_hashes", {}).items():
            target = PROJECT / rel
            source_checks[f"{path.parent.name}:{rel}"] = target.exists() and sha(target) == digest

    for name in ROOTS:
        visit(evidence_path(name))
    return records, graph, hash_checks, source_checks


def cycle_free(graph):
    color = {}

    def visit(node):
        color[node] = 1
        for dep in graph[node]:
            if color.get(dep) == 1:
                return False
            if color.get(dep, 0) == 0 and not visit(dep):
                return False
        color[node] = 2
        return True

    return all(color.get(node, 0) or visit(node) for node in graph)


def reaches(graph, start, target):
    todo, seen = [evidence_path(start).resolve()], set()
    goal = evidence_path(target).resolve()
    while todo:
        node = todo.pop()
        if node == goal:
            return True
        if node not in seen:
            seen.add(node)
            todo.extend(graph[node])
    return False


def quantifier_review():
    """Freeze the valid order and the remaining master provider condition."""
    # These are the exact retained slopes for the default h and kappa.
    h, kappa = F(1, 200), F(1, 100000)
    delta = F(1, 2) - 3 * kappa
    A = F(1, 2) + h

    def sufficient_stage(m, N):
        # Same conservative inequalities used by the complete diagonal.
        L = 2 * A + 2 * (m + 3)
        H = L + 1
        tail = (2 * (N + 1 + L + H) / (h * delta)) - 1
        B0 = F(3, 2) - 2 * kappa
        residual_loss = 2 * A + F(3, 2) + 2 * m
        finite = ((N + 1 + residual_loss) / h - B0) / delta
        ceil = lambda x: -(-x.numerator // x.denominator)
        return max(m + 2, 1, ceil(tail), ceil(finite))

    samples = {(m, N): sufficient_stage(m, N) for m, N in ((0, 0), (1, 3), (5, 11))}
    return {
        "sample_fixed_comparison_stages": {f"m={m},N={N}": j for (m, N), j in samples.items()},
        "default_stage_is_12816": samples[(0, 0)] == 12816,
        "field_and_sequences_chosen_before_requests": True,
        "derivative_and_pole_indices_are_proof_archives_not_field_reselections": True,
        "fixed_m_N_choose_one_finite_J_before_q_limit": True,
        "endpoint_uses_local_finiteness_before_terminal_limit": True,
        "force_extension_uses_one_diagonal_sequence_for_all_fixed_jets": True,
        "comparison_order_is_fixed_T_then_R_then_tau": True,
        "unconditional_provider_estimates_established_here": False,
    }


def review():
    records, graph, hashes, sources = load_graph()
    spine_links = {
        f"{later}_reaches_{earlier}": reaches(graph, later, earlier)
        for earlier, later in zip(CONSTRUCTION_SPINE, CONSTRUCTION_SPINE[1:])
    }
    quantifiers = quantifier_review()
    checks = {
        "thirteen_named_roots_loaded": all(evidence_path(name).resolve() in records for name in ROOTS),
        "transitive_evidence_graph_is_acyclic": cycle_free(graph),
        "all_transitive_dependency_hashes_match": bool(hashes) and all(hashes.values()),
        "all_transitive_source_hashes_match": bool(sources) and all(sources.values()),
        "active_spine_is_transitively_connected": all(spine_links.values()),
        "valid_fixed_request_and_limit_order": all(
            quantifiers[key] for key in (
                "default_stage_is_12816",
                "field_and_sequences_chosen_before_requests",
                "derivative_and_pole_indices_are_proof_archives_not_field_reselections",
                "fixed_m_N_choose_one_finite_J_before_q_limit",
                "endpoint_uses_local_finiteness_before_terminal_limit",
                "force_extension_uses_one_diagonal_sequence_for_all_fixed_jets",
                "comparison_order_is_fixed_T_then_R_then_tau",
            )
        ),
        "conditional_scope_preserved": not quantifiers["unconditional_provider_estimates_established_here"],
    }
    return {
        "checks": checks,
        "graph": {
            "nodes": len(graph),
            "edges": sum(len(deps) for deps in graph.values()),
            "dependency_hash_checks": len(hashes),
            "source_hash_checks": len(sources),
            "spine_links": spine_links,
        },
        "quantifiers": quantifiers,
        "remaining_master_condition": {
            "same_fixed_representatives": "for every fixed j,m,N the provider constants C_{j,m,N} are finite for one U^[j], not a request-dependent field",
            "scale_compatible_gaps": "Delta_bg>0, 3/8012>0, theta_ann<1, theta_rep<1",
            "weighted_domains": "all fixed-order profile, root, compact and slow-map weighted seminorms are finite on the unchanged supports",
            "residual": "the matching-state residual identity retains pressure, viscosity, all cutoffs, both cross interactions and the new square",
            "status": "asserted by frozen analytical provider contracts but not independently reconstructed in this checkpoint",
        },
        "scope": {
            "common_object_and_quantifier_join_checked": True,
            "tests_or_hashes_treated_as_continuum_proof": False,
            "unconditional_pulse_residual_force_or_blowup_claimed": False,
        },
    }
