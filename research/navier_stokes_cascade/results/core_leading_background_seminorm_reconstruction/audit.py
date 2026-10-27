#!/usr/bin/env python3
"""Freeze the direct leading/background seminorm reconstruction."""
import argparse
import hashlib
import json
from pathlib import Path
from seminorms import review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = (
    "core_current_head_conditional_frontier",
    "background_jet_foundation",
    "background_summation",
    "core_borel_prescription_closure",
    "core_leading_background_borel_join",
)
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def generate():
    result = review()
    dependencies = {str((HERE.parent / name / "evidence.json").relative_to(PROJECT)): sha(HERE.parent / name / "evidence.json") for name in DEPENDENCIES}
    inherited = {}
    for name in DEPENDENCIES:
        record = json.loads((HERE.parent / name / "evidence.json").read_text())
        for rel, digest in record.get("source_hashes", {}).items():
            target = PROJECT / rel
            inherited[f"{name}:{rel}"] = target.exists() and sha(target) == digest
    own = [HERE / name for name in ("README.md", "seminorms.py", "audit.py", "test_audit.py", "protocol.json", "requirements.txt", ".gitignore")]
    scope = result["scope"]
    return {
        "status": "fixed-order-background-seminorm-and-real-borel-selector-reconstructed-relative-to-leading-certificate",
        "protocol": json.loads((HERE / "protocol.json").read_text()),
        "review": result,
        "gates": {
            "five_exact_controls": len(result["checks"]) == 5 and all(result["checks"].values()),
            "five_frozen_dependencies": len(dependencies) == len(DEPENDENCIES),
            "all_inherited_sources_match": bool(inherited) and all(inherited.values()),
            "scientific_scope_preserved": scope["coefficient_and_borel_part_of_P1_reconstructed"] and not scope["scale_selected_before_coefficient_seminorm"] and not scope["numerical_infinite_scale_list_required"] and not scope["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"],
        },
        "dependency_evidence_sha256": dependencies,
        "inherited_source_checks": inherited,
        "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
        "remaining_obligations": [
            "prove the finalized leading-profile certificate directly from the same core, activation, continuation, modulation, repair and exterior defining maps",
            "in particular establish every fixed holomorphic-eta collar seminorm, the exact I3 power-law patch and ordered bumps, and the compact velocity/pressure-source identities on that one profile",
            "then advance the same-band gap contract without replacing continuum infima by sampled diagnostics",
        ],
    }

def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Seminorm reconstruction, provenance or scope changed")

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--verify-record", type=Path)
    args = parser.parse_args(argv)
    if args.output.exists(): raise SystemExit("Use a fresh output; frozen evidence is never overwritten")
    record = generate(); validate(record)
    if args.verify_record: validate(json.loads(args.verify_record.read_text()))
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print("leading/background seminorm reconstruction:", len(record["review"]["checks"]), "exact controls;", len(record["inherited_source_checks"]), "source checks")
if __name__ == "__main__": main()
