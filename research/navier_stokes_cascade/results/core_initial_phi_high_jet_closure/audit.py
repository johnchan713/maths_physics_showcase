import argparse, hashlib, json
from pathlib import Path
from controls import review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_initial_a1_request_graph",
    "core_axial_potential_pressure_gain",
    "core_finite_support_closure",
    "core_primary_tail_recheck",
    "core_initial_haar_row_budget",
    "core_rho_projector_signed_recompute",
    "core_initial_target_rate_gain",
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
    paths = [HERE/name for name in ("README.md", "closure.py", "controls.py", "audit.py",
                                    "test_audit.py", "requirements.txt", ".gitignore")]
    paths += [
        HERE.parent/"core_initial_a1_request_graph/requests.py",
        HERE.parent/"core_axial_potential_pressure_gain/bounds.py",
        HERE.parent/"core_symbolic_initial_haar_join/bounds.py",
        HERE.parent/"core_rho_projector_signed_recompute/bounds.py",
    ]
    return dict(
        status="conditional-first-real-Phi-high-jet-closure-from-all-fixed-order-DAG",
        review=result,
        gates=dict(exact_controls=all(result["checks"].values()),
                   seven_frozen_records=len(bound) == 7,
                   historical_sources=all(checks.values())),
        dependency_evidence_sha256=bound,
        inherited_source_checks=checks,
        source_hashes={str(path.relative_to(PROJECT)): sha(path) for path in paths},
        flags=dict(
            conditional_Phi_and_DzPhi_high_jet_implication_proved=True,
            independent_giant_numerical_jet_list_required=False,
            primitive_accepted_tail_contracts_independently_certified=False,
            radial_A0_and_complete_residual_proved=False,
            blowup_certified=False,
        ),
    )


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Initial Phi high-jet closure, provenance or scope changed")


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
    print("initial Phi high-jet closure:", len(record["review"]["checks"]),
          "exact controls;", len(record["inherited_source_checks"]), "historical source checks")


if __name__ == "__main__":
    main()
