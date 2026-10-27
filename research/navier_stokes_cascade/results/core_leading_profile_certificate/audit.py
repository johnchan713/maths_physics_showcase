#!/usr/bin/env python3
"""Freeze the defining-map leading-profile certificate."""
import argparse, hashlib, json
from pathlib import Path
from certificate import review

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
DEPENDENCIES = ("core_all_order_foundation", "core_activation_volterra_recheck",
 "prescribed_continuation_foundation", "core_modulation_repair_join_recheck",
 "annular_exterior_foundation", "core_global_pressure_moment_recheck",
 "core_coefficient_support_recheck", "core_leading_background_seminorm_reconstruction")
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def generate():
    result = review()
    dependencies = {str((HERE.parent/name/"evidence.json").relative_to(PROJECT)): sha(HERE.parent/name/"evidence.json") for name in DEPENDENCIES}
    inherited = {}
    for name in DEPENDENCIES:
        record = json.loads((HERE.parent/name/"evidence.json").read_text())
        for rel, digest in record.get("source_hashes", {}).items():
            target = PROJECT/rel
            inherited[f"{name}:{rel}"] = target.exists() and sha(target) == digest
    own = [HERE/name for name in ("README.md","certificate.py","audit.py","test_audit.py","protocol.json","requirements.txt",".gitignore")]
    analytical = result["analytical_join"]; scope = result["scope"]
    return {
      "status": "finalized-leading-profile-certificate-joined-relative-to-exact-core-and-root-hypotheses",
      "protocol": json.loads((HERE/"protocol.json").read_text()), "review": result,
      "gates": {
        "seven_exact_structural_controls": len(result["checks"]) == 7 and all(result["checks"].values()),
        "eight_defining_map_dependencies": len(dependencies) == len(DEPENDENCIES),
        "all_inherited_sources_match": bool(inherited) and all(inherited.values()),
        "pressure_and_analyticity_scope_preserved": analytical["full_pressure_memory_retained"] and not analytical["global_complex_eta_heat_endpoint_claimed"],
        "scientific_scope_preserved": scope["leading_leaf_of_first_background_contract_joined"] and not scope["same_band_continuum_infima_joined"] and not scope["unconditional_pulse_residual_force_exclusion_or_blowup_claimed"],
      },
      "dependency_evidence_sha256": dependencies, "inherited_source_checks": inherited,
      "source_hashes": {str(path.relative_to(PROJECT)): sha(path) for path in own},
      "remaining_obligations": [
        "audit the same-band gap contract on this completed leading profile using continuum infima rather than sampled diagnostics",
        "retain the flat-edge factorizations and produce one common selector for the core, determinant, normal and inward-cone gaps",
        "only then continue through the path ODE, projector, slow inverse, actual-pulse covariance and complete residual contracts",
      ]}

def validate(record):
    expected = generate()
    if record != expected or not all(expected["gates"].values()):
        raise ValueError("Leading-profile certificate, provenance or scope changed")

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--output",required=True,type=Path); p.add_argument("--verify-record",type=Path); a=p.parse_args(argv)
    if a.output.exists(): raise SystemExit("Use a fresh output; frozen evidence is never overwritten")
    record=generate(); validate(record)
    if a.verify_record: validate(json.loads(a.verify_record.read_text()))
    a.output.write_text(json.dumps(record,indent=2,sort_keys=True)+"\n")
    print("leading-profile certificate:",len(record["review"]["checks"]),"exact controls;",len(record["inherited_source_checks"]),"source checks")
if __name__ == "__main__": main()
