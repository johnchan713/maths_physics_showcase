#!/usr/bin/env python3
"""Freeze the same-band continuum-gap reconstruction."""
import argparse,hashlib,json
from pathlib import Path
from gaps import review
HERE=Path(__file__).resolve().parent; PROJECT=HERE.parent.parent
DEPENDENCIES=("core_leading_profile_certificate","core_leading_background_seminorm_reconstruction",
 "regional_cone_join","weighted_background_edges","core_background_finite_head_budget",
 "core_primary_tail_recheck","core_reference_scaled_normal_floor","core_same_band_gap_join")
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def generate():
 result=review(); records={}; inherited={}
 for name in DEPENDENCIES:
  path=HERE.parent/name/"evidence.json"; record=json.loads(path.read_text())
  records[str(path.relative_to(PROJECT))]=sha(path)
  for rel,digest in record.get("source_hashes",{}).items():
   source=PROJECT/rel; inherited[f"{name}:{rel}"]=source.exists() and sha(source)==digest
 own=[HERE/n for n in ("README.md","gaps.py","audit.py","test_audit.py","protocol.json","requirements.txt",".gitignore")]
 scope=result["scope"]
 return {"status":"same-band-continuum-gap-contract-reconstructed-on-finalized-leading-profile",
  "protocol":json.loads((HERE/"protocol.json").read_text()),"review":result,
  "gates":{"six_exact_controls":len(result["checks"])==6 and all(result["checks"].values()),
   "eight_frozen_dependencies":len(records)==len(DEPENDENCIES),
   "all_inherited_sources_match":bool(inherited) and all(inherited.values()),
   "edge_correction_retained":not scope["stress_magnitude_positive_at_flat_edges"],
   "scientific_scope_retained":scope["same_band_contract_reconstructed"] and not scope["full_path_ode_closed"] and not scope["unconditional_pulse_residual_force_or_blowup"]},
  "dependency_evidence_sha256":records,"inherited_source_checks":inherited,
  "source_hashes":{str(p.relative_to(PROJECT)):sha(p) for p in own},
  "remaining_obligations":[
   "audit the canonical homogeneous full-path ODE on this accepted tail with moving normal, varying growth, unprojected pressure source and every cutoff derivative",
   "then audit the fixed compensator/projector supports and the compact and slow inverse weighted domains on the same field",
   "only after those contracts may actual-pulse covariance, fresh axial/full-pressure means and the complete residual be promoted"]}
def validate(record):
 expected=generate()
 if record!=expected or not all(expected["gates"].values()): raise ValueError("Same-band reconstruction, provenance or scope changed")
def main(argv=None):
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",required=True,type=Path);p.add_argument("--verify-record",type=Path);a=p.parse_args(argv)
 if a.output.exists(): raise SystemExit("Use a fresh output; frozen evidence is never overwritten")
 record=generate();validate(record)
 if a.verify_record:validate(json.loads(a.verify_record.read_text()))
 a.output.write_text(json.dumps(record,indent=2,sort_keys=True)+"\n")
 print("same-band gap reconstruction:",len(record["review"]["checks"]),"exact controls;",len(record["inherited_source_checks"]),"source checks")
if __name__=="__main__":main()
