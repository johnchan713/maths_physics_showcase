import argparse, hashlib, json
from pathlib import Path
from controls import review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_initial_physical_stage_sync",
    "core_physical_tail_join",
    "first_sourced_foundation",
    "core_final_mean_row_closure",
    "core_reweighted_fresh_mean_closure",
    "core_sourced_covariance_moment_join",
    "regular_source_foundation",
    "core_finite_support_closure",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    frozen, inherited = {}, {}
    for name in DEPENDENCIES:
        path = HERE.parent/name/"evidence.json"
        record = json.loads(path.read_text())
        frozen[str(path.relative_to(PROJECT))] = sha(path)
        for rel, digest in record["source_hashes"].items():
            inherited[name+":"+rel] = sha(PROJECT/rel) == digest
    result = json.loads(json.dumps(review(), default=str))
    paths = [HERE/name for name in (
        "README.md", "closure.py", "controls.py", "audit.py",
        "test_audit.py", "requirements.txt", ".gitignore")]
    paths += [
        HERE.parent/"core_physical_tail_join/bounds.py",
        HERE.parent/"core_final_mean_row_closure/budgets.py",
        HERE.parent/"core_reweighted_fresh_mean_closure/budgets.py",
        HERE.parent/"core_sourced_covariance_moment_join/budgets.py",
    ]
    return dict(
        status="conditional-first-sourced-same-state-compact-physical-tail-closure",
        review=result,
        gates=dict(exact_controls=all(result["checks"].values()),
                   eight_frozen_records=len(frozen) == 8,
                   historical_sources=all(inherited.values())),
        dependency_evidence_sha256=frozen,
        inherited_source_checks=inherited,
        source_hashes={str(path.relative_to(PROJECT)): sha(path) for path in paths},
        flags=dict(
            first_sourced_compact_graph_acyclic=True,
            changed_covariance_pressure_and_axial_velocity_retained=True,
            fixed_m_N_post_sourced_compact_flatness_proved_conditionally=True,
            actual_high_jet_values_materialized=False,
            primitive_contracts_independently_certified=False,
            complete_regular_residual_decay_proved=False,
            blowup_certified=False,
        ),
    )


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("First-sourced compact-tail closure, provenance or scope changed")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify-record", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("Use a fresh output; frozen evidence is never overwritten")
    record = generate()
    validate(record)
    if args.verify_record:
        validate(json.loads(args.verify_record.read_text()))
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True)+"\n")
    print("first-sourced physical-tail closure:",
          len(record["review"]["checks"]), "exact controls;",
          len(record["inherited_source_checks"]), "historical source checks")


if __name__ == "__main__":
    main()
