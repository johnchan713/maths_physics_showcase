#!/usr/bin/env python3
"""Freeze the current-head conditional-frontier audit."""
import argparse
import hashlib
import json
from pathlib import Path
from current_frontier import FRONTIER, START, review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def generate():
    result = review()
    dependency = HERE.parent / START / "evidence.json"
    own = [HERE / name for name in ("README.md", "current_frontier.py", "audit.py", "test_audit.py", "protocol.json", "requirements.txt", ".gitignore")]
    scope = result["scope"]
    return {
        "status": "current-head-downstream-replays-audited-with-five-contract-continuum-frontier",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": result,
        "gates": {
            "six_graph_and_frontier_controls": len(result["checks"]) == 6 and all(result["checks"].values()),
            "one_current_head_dependency": dependency.exists(),
            "five_live_contracts_enumerated": len(result["live_frontier"]) == len(FRONTIER),
            "scientific_scope_preserved": scope["downstream_chain_complete_relative_to_five_contracts"] and not scope["five_contracts_independently_reproved_in_this_checkpoint"] and not scope["hashes_or_tests_treated_as_continuum_proof"] and not scope["unconditional_pulse_force_exclusion_or_blowup_claimed"],
        },
        "dependency_evidence_sha256": {str(dependency.relative_to(PROJECT)): sha(dependency)},
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "independently audit P1 from the defining leading profile, coefficient recursion and real Borel selector, proving every fixed seminorm is finite before scale selection",
            "then audit P2 through P5 on that same object without substituting a diagnostic profile, sampled supremum or historical first band",
            "only after those five continuum contracts are discharged may the already replayed residual, force and comparison implications be considered for unconditional promotion",
        ],
    }

def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Current-head graph, frontier, provenance or scope changed")

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--verify-record", type=Path)
    args = parser.parse_args(argv)
    if args.output.exists(): raise SystemExit("Use a fresh output; frozen evidence is never overwritten")
    record = generate(); validate(record)
    if args.verify_record: validate(json.loads(args.verify_record.read_text()))
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    graph = record["review"]["graph"]
    print("current-head frontier audit:", graph["nodes"], "nodes;", graph["edges"], "edges;", graph["source_hash_checks"], "source checks")
if __name__ == "__main__": main()
