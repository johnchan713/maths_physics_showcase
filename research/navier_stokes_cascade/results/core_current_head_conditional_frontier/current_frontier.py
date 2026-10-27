"""Current-head dependency closure and live continuum frontier."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
PROJECT = RESULTS.parent
START = "core_log_reserved_pressure_comparison_replay"

FRONTIER = (
    "core_leading_background_borel_join",
    "core_same_band_gap_join",
    "core_canonical_full_path_ode_join",
    "core_fixed_profile_projector_join",
    "core_compact_slow_inverse_domain_join",
)

DISCHARGE_CHAIN = (
    "core_complete_dependency_quantifier_audit",
    "core_direct_compact_slow_provider_reconstruction",
    "core_log_renormalized_actual_source_induction",
    "core_log_reserved_complete_residual_replay",
    "core_log_reserved_cartesian_corner_replay",
    "core_log_reserved_endpoint_compatibility_join",
    "core_log_reserved_localization_force_replay",
    "core_log_reserved_preterminal_energy_replay",
    "core_log_reserved_pressure_comparison_replay",
)

def evidence_path(name):
    return RESULTS / name / "evidence.json"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def load_graph(start=START):
    """Load the transitive evidence DAG and verify every declared binding."""
    records, graph, dependency_checks, source_checks = {}, {}, {}, {}
    def visit(path):
        path = path.resolve()
        if path in records:
            return
        record = json.loads(path.read_text())
        records[path] = record
        dependencies = []
        for rel, digest in record.get("dependency_evidence_sha256", {}).items():
            target = (PROJECT / rel).resolve()
            key = f"{path.parent.name}:{rel}"
            dependency_checks[key] = target.exists() and sha(target) == digest
            if target.exists():
                dependencies.append(target)
                visit(target)
        graph[path] = dependencies
        for rel, digest in record.get("source_hashes", {}).items():
            target = PROJECT / rel
            source_checks[f"{path.parent.name}:{rel}"] = target.exists() and sha(target) == digest
    visit(evidence_path(start))
    return records, graph, dependency_checks, source_checks

def cycle_free(graph):
    color = {}
    def visit(node):
        color[node] = 1
        for dependency in graph[node]:
            if color.get(dependency) == 1:
                return False
            if color.get(dependency, 0) == 0 and not visit(dependency):
                return False
        color[node] = 2
        return True
    return all(color.get(node, 0) or visit(node) for node in graph)

def reaches(graph, later, earlier):
    todo = [evidence_path(later).resolve()]
    target = evidence_path(earlier).resolve()
    seen = set()
    while todo:
        node = todo.pop()
        if node == target:
            return True
        if node not in seen:
            seen.add(node)
            todo.extend(graph.get(node, ()))
    return False

def frontier_contracts():
    """State the five same-object propositions still imported downstream."""
    return {
        "P1_leading_and_borel": "for the same completed leading profile, every fixed coefficient seminorm C_bg(n,m) is finite before the real Borel scale b_n is selected, and the selected locally finite sum retains full pressure",
        "P2_same_band_gaps": "on one derivative-independent accepted tail, Delta_bg>0, the phase normal has a positive floor, and the value, spectral and factored cone gaps hold for the P1 field",
        "P3_full_path_ode": "the canonical full-path propagator has finite fixed-jet bounds on the P2 tail while retaining n dot f in pressure, moving normal, varying growth and every cutoff derivative",
        "P4_fixed_profiles": "one fixed compensator and the same-choice projectors have finite fixed-order seminorms, exact moments and common reserved supports, with all projector-dependent pressure objects recomputed",
        "P5_weighted_maps": "the shifted two-edge compact maps, zero-Haar inverse and fixed slow response preserve zeta times delta^-p weighted fixed-order domains and complete pressure memory on the P1--P4 field",
    }

def review():
    records, graph, dependencies, sources = load_graph()
    chain_links = {f"{later}_reaches_{earlier}": reaches(graph, later, earlier) for earlier, later in zip(DISCHARGE_CHAIN, DISCHARGE_CHAIN[1:])}
    frontier_links = {f"current_head_reaches_{name}": reaches(graph, START, name) for name in FRONTIER}
    checks = {
        "current_head_graph_is_acyclic": cycle_free(graph),
        "all_dependency_hashes_match": bool(dependencies) and all(dependencies.values()),
        "all_source_hashes_match": bool(sources) and all(sources.values()),
        "discharge_chain_is_transitively_connected": all(chain_links.values()),
        "five_primitive_contracts_are_in_current_closure": all(frontier_links.values()),
        "old_provider_frontier_is_followed_by_log_source_and_residual_replays": reaches(graph, "core_log_reserved_complete_residual_replay", "core_complete_dependency_quantifier_audit") and reaches(graph, "core_log_reserved_complete_residual_replay", "core_log_renormalized_actual_source_induction"),
    }
    return {
        "checks": checks,
        "graph": {"nodes": len(graph), "edges": sum(len(row) for row in graph.values()), "dependency_hash_checks": len(dependencies), "source_hash_checks": len(sources), "chain_links": chain_links, "frontier_links": frontier_links},
        "discharged_relative_to_frontier": {
            "pure_power_source_frontier": "corrected: finite logarithmic degree is allowed and absorbed by one additional compact integration",
            "compact_and_slow_operator_bounds": "reconstructed directly",
            "actual_source_induction": "proved in the log-weighted domain for every fixed stage",
            "complete_residual": "replayed with all six groups and separate logarithmic reserves",
            "terminal_corner_and_endpoint": "replayed on the same fixed-gauge field",
            "localization_and_force": "all twelve localization terms and compact future force replayed",
            "energy": "exact preterminal identity and separate endpoint bounds replayed",
            "pressure_and_comparison": "force divergence retained and expanding-ball remainders replayed",
        },
        "live_frontier": frontier_contracts(),
        "quantifier_guard": {
            "one_field_profile_scales_gauges_and_supports_fixed_first": True,
            "fixed_m_N_then_finite_stage_then_small_scale_limit": True,
            "comparison_order_fixed_T_then_radius_then_core": True,
            "no_stage_uniform_constant_is_needed": True,
        },
        "scope": {
            "downstream_chain_complete_relative_to_five_contracts": True,
            "five_contracts_independently_reproved_in_this_checkpoint": False,
            "hashes_or_tests_treated_as_continuum_proof": False,
            "unconditional_pulse_force_exclusion_or_blowup_claimed": False,
        },
    }
