import argparse, hashlib, json
from pathlib import Path
from controls import review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_initial_target_rate_gain",
    "core_initial_final_row_obstruction",
    "core_axial_potential_pressure_gain",
    "core_rho_projector_signed_recompute",
    "core_finite_support_closure",
    "core_real_refresh_obligation_join",
    "core_initial_haar_row_budget",
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
    paths = [HERE/name for name in ("README.md", "requests.py", "controls.py", "audit.py",
                                    "test_audit.py", "requirements.txt", ".gitignore")]
    paths += [
        HERE.parent/"core_initial_final_row_obstruction/bounds.py",
        HERE.parent/"core_axial_potential_pressure_gain/bounds.py",
        HERE.parent/"core_rho_projector_signed_recompute/bounds.py",
        HERE.parent/"core_symbolic_initial_haar_join/bounds.py",
    ]
    return dict(
        status="conditional-corrected-first-real-A1-finite-request-graph",
        review=result,
        gates=dict(exact_controls=all(result["checks"].values()),
                   seven_frozen_records=len(bound) == 7,
                   historical_sources=all(checks.values())),
        dependency_evidence_sha256=bound,
        inherited_source_checks=checks,
        source_hashes={str(path.relative_to(PROJECT)): sha(path) for path in paths},
        flags=dict(
            immediate_A1_request_graph_reduced_and_padded=True,
            qualitative_profile_and_support_closure_reused=True,
            actual_weighted_high_jet_field_inequalities_proved=False,
            radial_A0_complete_residual_tail_proved=False,
            complete_residual_and_blowup_certified=False,
        ),
    )


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Initial A1 request graph, provenance or scope changed")


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
    print("initial A1 request graph:", len(record["review"]["checks"]),
          "exact controls;", len(record["inherited_source_checks"]), "historical source checks")


if __name__ == "__main__":
    main()
