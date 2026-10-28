"""Current-head dependency closure and lower primitive frontier."""

import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
PROJECT = RESULTS.parent
START = "core_current_log_reserved_pressure_comparison_join"

OLD_FRONTIER = (
    "core_leading_background_borel_join",
    "core_same_band_gap_join",
    "core_canonical_full_path_ode_join",
    "core_fixed_profile_projector_join",
    "core_compact_slow_inverse_domain_join",
)

RECONSTRUCTIONS = (
    "core_leading_background_seminorm_reconstruction",
    "core_leading_profile_certificate",
    "core_same_band_gap_reconstruction",
    "core_full_path_ode_reconstruction",
    "core_fixed_profile_projector_reconstruction",
    "core_compact_slow_inverse_reconstruction",
)

CURRENT_REPLAYS = (
    "core_current_log_source_reconciliation",
    "core_current_log_reserved_residual_reconciliation",
    "core_current_log_reserved_corner_endpoint_join",
    "core_current_log_reserved_localization_force_join",
    "core_current_log_reserved_preterminal_energy_join",
    "core_current_log_reserved_pressure_comparison_join",
)

BOUNDARY_WITNESSES = (
    "core_leading_profile_certificate",
    "core_same_band_gap_reconstruction",
    "core_reference_scaled_normal_floor",
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


def selector_ceiling(h, constants, margins):
    """Return the gap-row ceiling, or fail if its continuum premises are absent."""
    if not math.isfinite(h) or h <= 0 or len(constants) != len(margins) or not constants:
        raise ValueError("positive h and paired nonempty rows are required")
    rows = []
    for constant, margin in zip(constants, margins):
        if not math.isfinite(constant) or constant <= 0:
            raise ValueError("every actual perturbation constant must be finite and positive")
        if not math.isfinite(margin) or margin <= 0:
            raise ValueError("every actual continuum margin must be finite and positive")
        rows.append(math.log2(constant / margin))
    return max(2, math.ceil(2 + max(rows) / (2 * h)))


def primitive_frontier():
    return {
        "Q1_exact_core_and_roots": (
            "the stated pressure-family Banach fixed point, completed global backward-pressure identification, "
            "and every finite implicit repair/matching root satisfy their exact defining equations on the same profile"
        ),
        "Q2_regional_continuum_margins": (
            "on the actual finite regional cover, f=inf F0, t=inf(-g0_theta), d=inf D0, the orientation and "
            "same-tilt cone margins, and the two factored edge coefficients are strictly positive"
        ),
        "Q3_full_frame_primitive_domain": (
            "the original background/frame primitive bounds, complete path and phase jets, positive carrier/root floors, "
            "and the physical accepted-tail domain hold for that same field with finite continuum constants"
        ),
    }


def review():
    records, graph, dependencies, sources = load_graph()
    reconstruction_links = {
        f"current_head_reaches_{name}": reaches(graph, START, name) for name in RECONSTRUCTIONS
    }
    replay_links = {
        f"{later}_reaches_{earlier}": reaches(graph, later, earlier)
        for earlier, later in zip(CURRENT_REPLAYS, CURRENT_REPLAYS[1:])
    }
    old_links = {f"current_head_reaches_historical_{name}": reaches(graph, START, name) for name in OLD_FRONTIER}
    boundary_links = {
        f"current_head_reaches_boundary_witness_{name}": reaches(graph, START, name)
        for name in BOUNDARY_WITNESSES
    }
    checks = {
        "current_head_graph_is_acyclic": cycle_free(graph),
        "all_dependency_hashes_match": bool(dependencies) and all(dependencies.values()),
        "all_source_hashes_match": bool(sources) and all(sources.values()),
        "all_six_reconstructions_are_in_current_closure": all(reconstruction_links.values()),
        "current_replay_chain_is_transitively_connected": all(replay_links.values()),
        "historical_five_contract_records_remain_in_closure": all(old_links.values()),
        "lower_boundary_witnesses_are_in_current_closure": all(boundary_links.values()),
    }
    return {
        "checks": checks,
        "graph": {
            "nodes": len(graph),
            "edges": sum(len(row) for row in graph.values()),
            "dependency_hash_checks": len(dependencies),
            "source_hash_checks": len(sources),
            "reconstruction_links": reconstruction_links,
            "replay_links": replay_links,
            "historical_frontier_links": old_links,
            "boundary_links": boundary_links,
        },
        "correction": {
            "historical_five_contract_frontier_is_live": False,
            "same_object_reconstructions_present": list(RECONSTRUCTIONS),
            "downstream_current_replays_present": list(CURRENT_REPLAYS),
        },
        "live_frontier": primitive_frontier(),
        "blocking_inequality": {
            "aggregate_margin": "m_base=min(f,t,d,omega,m_phi,c_a,b_0,lambda_min,b_min^2)>0",
            "finite_constants": "0<C_i<infinity for every value, relative-stress, chart, length and frame row",
            "selector_row": "2*h*(ell-2)>=max_i log2(C_i/r_i)",
            "logical_status": "ell is not defined by this row until every actual r_i is positive and every C_i is finite; q-smallness cannot prove those antecedents",
        },
        "quantifier_guard": {
            "one_exact_profile_roots_supports_and_gauges_fixed_first": True,
            "actual_continuum_infima_precede_selector_choice": True,
            "one_common_selector_tail_precedes_band_limit": True,
            "fixed_jet_stage_and_target_precede_small_scale_limit": True,
            "comparison_order_is_fixed_T_then_R_then_terminal_core": True,
            "no_stage_uniform_constant_is_needed": True,
        },
        "scope": {
            "former_five_contracts_reconstructed_on_current_objects": True,
            "downstream_replay_complete_relative_to_lower_frontier": True,
            "Q1_Q2_Q3_independently_reproved_here": False,
            "hashes_or_tests_treated_as_continuum_proof": False,
            "unconditional_pulse_residual_force_exclusion_or_blowup_claimed": False,
        },
    }
