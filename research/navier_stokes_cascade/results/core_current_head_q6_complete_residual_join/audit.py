#!/usr/bin/env python3
"""Freeze the current-head Q6 complete six-group residual join."""

import argparse
import hashlib
import json
from pathlib import Path

from complete_residual import DIRECT_DEPENDENCIES, review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "current-head-Q6-complete-six-group-local-residual-bound-conditionally"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_closure():
    records, graph, dependency_checks, source_checks = {}, {}, {}, {}

    def visit(path):
        path = path.resolve()
        if path in records:
            return
        record = json.loads(path.read_text())
        records[path] = record
        children = []
        for relative, digest in record.get("dependency_evidence_sha256", {}).items():
            child = (PROJECT / relative).resolve()
            dependency_checks[f"{path.parent.name}:{relative}"] = (
                child.exists() and sha(child) == digest
            )
            if child.exists():
                children.append(child)
                visit(child)
        graph[path] = children
        for relative, digest in record.get("source_hashes", {}).items():
            source = PROJECT / relative
            source_checks[f"{path.parent.name}:{relative}"] = (
                source.exists() and sha(source) == digest
            )

    for name in DIRECT_DEPENDENCIES:
        visit(HERE.parent / name / "evidence.json")
    return records, graph, dependency_checks, source_checks


def cycle_free(graph):
    colour = {}

    def visit(node):
        colour[node] = 1
        for child in graph[node]:
            if colour.get(child) == 1:
                return False
            if colour.get(child, 0) == 0 and not visit(child):
                return False
        colour[node] = 2
        return True

    return all(colour.get(node, 0) or visit(node) for node in graph)


def generate():
    result = review()
    records, graph, dependencies, sources = load_closure()
    own = [HERE / name for name in (
        "README.md", "complete_residual.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    providers = [
        HERE.parent / "core_current_head_q5_weighted_inverse_source_join/weighted_inverse.py",
        HERE.parent / "core_current_log_reserved_residual_reconciliation/residual.py",
        HERE.parent / "core_log_reserved_complete_residual_replay/replay.py",
    ]
    scope = result["scope"]
    gates = {
        "nine_exact_Q6_controls": len(result["checks"]) == 9 and all(result["checks"].values()),
        "six_direct_provider_records": len(DIRECT_DEPENDENCIES) == 6 and all(
            (HERE.parent / name / "evidence.json").exists()
            for name in DIRECT_DEPENDENCIES
        ),
        "provider_graph_is_acyclic": cycle_free(graph),
        "all_dependency_hashes_match": bool(dependencies) and all(dependencies.values()),
        "all_inherited_source_hashes_match": bool(sources) and all(sources.values()),
        "scientific_scope_preserved": (
            scope["Q1_Q5_same_object_is_directly_bound_to_residual_replay"]
            and scope["every_fixed_jet_and_target_has_one_finite_comparison_stage"]
            and scope["complete_local_residual_is_flat_conditionally"]
            and not scope["disjoint_reserves_are_added_into_one_stronger_common_margin"]
            and not scope["stage_uniform_constant_or_J_of_q_used"]
            and not scope["corner_endpoint_localization_force_energy_comparison_replayed_here"]
            and not scope["hashes_or_tests_treated_as_continuum_proof"]
            and not scope["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"]
        ),
    }
    return {
        "status": STATUS,
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": result,
        "graph": {
            "nodes": len(graph),
            "edges": sum(len(children) for children in graph.values()),
            "dependency_hash_checks": len(dependencies),
            "source_hash_checks": len(sources),
        },
        "gates": gates,
        "dependency_evidence_sha256": {
            str((HERE.parent / name / "evidence.json").relative_to(PROJECT)):
                sha(HERE.parent / name / "evidence.json")
            for name in DIRECT_DEPENDENCIES
        },
        "source_hashes": {
            str(path.relative_to(PROJECT)): sha(path) for path in own + providers
        },
        "remaining_obligations": [
            "transfer the Q6 residual to the Cartesian terminal corner and join the same-field punctured endpoint jets",
            "replay all twelve localization groups and the compact-force extension with unit viscosity and full pressure cutoff derivatives",
            "only then revalidate energy, fixed-gauge pressure comparison and any conditional exclusion consequence",
            "do not promote an unconditional pulse, residual, force, exclusion or blow-up theorem without all continuum premises",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Q6 residual interfaces, provenance, arithmetic or scope changed")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--verify-record", type=Path)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit("Use a fresh output; frozen evidence is never overwritten")
    record = generate()
    validate(record)
    if args.verify_record:
        validate(json.loads(args.verify_record.read_text()))
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    graph = record["graph"]
    print(
        "current-head Q6 residual join:", graph["nodes"], "nodes;",
        graph["edges"], "edges;", graph["source_hash_checks"], "source checks",
    )


if __name__ == "__main__":
    main()
