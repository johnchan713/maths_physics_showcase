#!/usr/bin/env python3
"""Freeze the current-head Q9 preterminal energy join."""

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from preterminal_energy import DIRECT_DEPENDENCIES, review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = "current-head-Q9-preterminal-energy-and-dissipation-joined-conditionally"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    if isinstance(value, Fraction):
        return str(value)
    return value


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
    result = encode(review())
    records, graph, dependencies, sources = load_closure()
    own = [HERE / name for name in (
        "README.md", "preterminal_energy.py", "audit.py", "test_audit.py",
        "protocol.json", "requirements.txt", ".gitignore",
    )]
    providers = [
        HERE.parent /
        "core_current_head_q8_localization_force_join/localization_force.py",
        HERE.parent /
        "core_current_log_reserved_preterminal_energy_join/preterminal_energy.py",
        HERE.parent /
        "core_log_reserved_preterminal_energy_replay/preterminal_energy.py",
    ]
    scope = result["scope"]
    gates = {
        "twelve_exact_Q9_controls":
            len(result["checks"]) == 12 and all(result["checks"].values()),
        "six_direct_provider_records":
            len(DIRECT_DEPENDENCIES) == 6 and all(
                (HERE.parent / name / "evidence.json").exists()
                for name in DIRECT_DEPENDENCIES
            ),
        "provider_graph_is_acyclic": cycle_free(graph),
        "all_dependency_hashes_match":
            bool(dependencies) and all(dependencies.values()),
        "all_inherited_source_hashes_match":
            bool(sources) and all(sources.values()),
        "scientific_scope_preserved": (
            scope[
                "exact_Q1_Q8_preterminal_energy_identity_replayed_conditionally"
            ]
            and scope[
                "separate_finite_energy_and_total_dissipation_bounds_obtained_conditionally"
            ]
            and not scope["Q8_endpoint_uniform_all_order_bound_proved_here"]
            and not scope["terminal_velocity_pressure_or_energy_trace_proved"]
            and not scope["future_velocity_or_pressure_constructed"]
            and not scope["fixed_gauge_pressure_recovery_replayed_here"]
            and not scope[
                "expanding_ball_comparison_or_competitor_exclusion_proved"
            ]
            and not scope["hashes_or_tests_treated_as_continuum_proof"]
            and not scope[
                "unconditional_pulse_residual_force_exclusion_or_blowup_claimed"
            ]
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
            "prove the Q8 punctured endpoint inequality for the actual A B and fixed-gauge pressure providers",
            "replay fixed-gauge pressure recovery with -div(f) retained for this exact finite-energy field and force",
            "then prove expanding-ball comparison for every admissible smooth same-force finite-energy competitor",
            "do not infer terminal field regularity or an unconditional exclusion or blow-up theorem",
        ],
    }


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Q9 energy, provenance, endpoint scope or limits changed")


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
        "current-head Q9 energy join:",
        graph["nodes"], "nodes;", graph["edges"], "edges;",
        graph["source_hash_checks"], "source checks",
    )


if __name__ == "__main__":
    main()
