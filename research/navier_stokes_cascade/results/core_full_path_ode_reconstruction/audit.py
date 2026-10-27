#!/usr/bin/env python3
"""Freeze the canonical full-path ODE reconstruction."""
import argparse,hashlib,json
from pathlib import Path
from ode import review
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=("core_same_band_gap_reconstruction","core_phase_coefficient_budget",
 "core_tangent_propagator_closure","canonical_primary_foundation",
 "core_canonical_entrance_jet_closure","core_canonical_haar_column_jet_closure",
 "sourced_inverse_foundation","core_finite_support_closure","core_canonical_full_path_ode_join")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def generate():
 result=review();records={};inherited={}
 for name in DEPENDENCIES:
  p=HERE.parent/name/"evidence.json";r=json.loads(p.read_text());records[str(p.relative_to(PROJECT))]=sha(p)
  for rel,digest in r.get("source_hashes",{}).items():
   source=PROJECT/rel;inherited[f"{name}:{rel}"]=source.exists() and sha(source)==digest
 own=[HERE/n for n in ("README.md","ode.py","audit.py","test_audit.py","protocol.json","requirements.txt",".gitignore")];scope=result["scope"]
 return {"status":"canonical-full-path-ode-reconstructed-on-current-same-band-tail",
  "protocol":json.loads((HERE/"protocol.json").read_text()),"review":result,
  "gates":{"ten_exact_controls":len(result["checks"])==10 and all(result["checks"].values()),
   "nine_frozen_dependencies":len(records)==len(DEPENDENCIES),"all_inherited_sources_match":bool(inherited) and all(inherited.values()),
   "primary_source_distinction_retained":scope["canonical_primary_is_homogeneous"] and scope["sourced_inverse_retains_unprojected_pressure"],
   "scientific_scope_retained":not scope["raw_exit_response_declared_zero"] and not scope["fixed_projector_contract_closed"] and not scope["unconditional_pulse_residual_force_or_blowup"]},
  "dependency_evidence_sha256":records,"inherited_source_checks":inherited,
  "source_hashes":{str(p.relative_to(PROJECT)):sha(p) for p in own},
  "remaining_obligations":["audit fixed compensator/projector profiles, actual supports and every fixed jet on this same accepted domain",
   "join the compact and slow inverse maps with their weighted domains and whole-path support",
   "only then join complete covariance, fresh axial/full-pressure means and every viscosity, cutoff and interaction residual term"]}
def validate(record):
 expected=generate()
 if record!=expected or not all(expected["gates"].values()):raise ValueError("Full-path ODE reconstruction, provenance or scope changed")
def main(argv=None):
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",required=True,type=Path);p.add_argument("--verify-record",type=Path);a=p.parse_args(argv)
 if a.output.exists():raise SystemExit("Use a fresh output; frozen evidence is never overwritten")
 record=generate();validate(record)
 if a.verify_record:validate(json.loads(a.verify_record.read_text()))
 a.output.write_text(json.dumps(record,indent=2,sort_keys=True)+"\n");print("full-path ODE reconstruction:",len(record["review"]["checks"]),"exact controls;",len(record["inherited_source_checks"]),"source checks")
if __name__=="__main__":main()
