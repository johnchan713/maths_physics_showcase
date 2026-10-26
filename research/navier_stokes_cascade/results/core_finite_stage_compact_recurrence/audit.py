import argparse
import hashlib
import json
from pathlib import Path
from controls import review


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_post_signed_seed_sync",
    "core_first_sourced_physical_tail_closure",
    "core_post_signed_mean_closure",
    "core_post_signed_residual_join",
    "core_next_sourced_wave_join",
    "core_next_sourced_mean_tail_join",
    "finite_stage_induction",
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
    sources = [HERE/name for name in (
        "README.md", "recurrence.py", "controls.py", "audit.py",
        "test_audit.py", "requirements.txt", ".gitignore")]
    sources += [
        HERE.parent/"core_post_signed_seed_sync/ledger.py",
        HERE.parent/"core_post_signed_mean_closure/budgets.py",
        HERE.parent/"core_next_sourced_mean_tail_join/budgets.py",
        HERE.parent/"finite_stage_induction/induction.py",
    ]
    return dict(
        status="conditional-current-seed-finite-stage-compact-recurrence",
        review=result,
        gates=dict(exact_controls=all(result["checks"].values()),
                   eight_frozen_records=len(frozen) == 8,
                   historical_sources=all(inherited.values())),
        dependency_evidence_sha256=frozen,
        inherited_source_checks=inherited,
        source_hashes={str(path.relative_to(PROJECT)): sha(path) for path in sources},
        flags=dict(
            corrected_seed_used=True,
            finite_stage_backward_request_recurrence_supplied=True,
            distinct_stage_pole_and_support_archive_supplied=True,
            all_fixed_stage_compact_tails_flat_conditionally=True,
            actual_primitive_contracts_independently_certified=False,
            stage_uniform_constants_or_fixed_q_contraction_certified=False,
            one_field_all_targets_diagonal_proved=False,
            localization_forcing_energy_or_blowup_proved=False,
        ),
    )


def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Finite-stage compact recurrence, provenance or scope changed")


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
    print("finite-stage compact recurrence:",
          len(record["review"]["checks"]), "exact controls;",
          len(record["inherited_source_checks"]), "historical source checks")


if __name__ == "__main__":
    main()
