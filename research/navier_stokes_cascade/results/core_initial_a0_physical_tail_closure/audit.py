import argparse, hashlib, json
from pathlib import Path
from controls import review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_initial_phi_high_jet_closure",
    "core_initial_a1_request_graph",
    "core_initial_target_rate_gain",
    "core_real_refresh_obligation_join",
    "core_post_signed_mean_closure",
    "core_physical_tail_join",
    "core_compact_input_sharpness",
    "core_finite_support_closure",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    bound, checks = {}, {}
    for name in DEPENDENCIES:
        path = HERE.parent/name/"evidence.json"
        record = json.loads(path.read_text())
        bound[str(path.relative_to(PROJECT))] = sha(path)
        for rel, digest in record["source_hashes"].items():
            checks[name+":"+rel] = sha(PROJECT/rel) == digest
    result = json.loads(json.dumps(review(), default=str))
    paths = [HERE/name for name in (
        "README.md", "closure.py", "controls.py", "audit.py",
        "test_audit.py", "requirements.txt", ".gitignore")]
    paths += [
        HERE.parent/"core_initial_phi_high_jet_closure/closure.py",
        HERE.parent/"core_initial_a1_request_graph/requests.py",
        HERE.parent/"core_post_signed_mean_closure/budgets.py",
        HERE.parent/"core_physical_tail_join/bounds.py",
    ]
    return dict(
        status="conditional-first-real-A0-and-compact-physical-tail-closure",
        review=result,
        gates=dict(exact_controls=all(result["checks"].values()),
                   eight_frozen_records=len(bound) == 8,
                   historical_sources=all(checks.values())),
        dependency_evidence_sha256=bound,
        inherited_source_checks=checks,
        source_hashes={str(path.relative_to(PROJECT)): sha(path) for path in paths},
        flags=dict(
            conditional_A0_source_and_tail_implication_proved=True,
            fixed_m_N_compact_physical_tail_interface_proved_conditionally=True,
            independent_giant_numerical_jet_list_required=False,
            primitive_accepted_tail_contracts_independently_certified=False,
            complete_regular_residual_improvement_proved=False,
            blowup_certified=False,
        ),
    )


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Initial A0/physical-tail closure, provenance or scope changed")


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
    print("initial A0 physical-tail closure:", len(record["review"]["checks"]),
          "exact controls;", len(record["inherited_source_checks"]), "historical source checks")


if __name__ == "__main__":
    main()
