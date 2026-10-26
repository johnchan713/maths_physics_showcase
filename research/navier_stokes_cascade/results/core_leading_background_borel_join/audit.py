#!/usr/bin/env python3
"""Freeze the leading-profile/background/Borel contract join."""
import argparse,hashlib,json
from pathlib import Path
from borel_join import review

HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
DEPENDENCIES=(
 "core_corrected_foundation_frontier","core_fixed_point_ledger_recheck",
 "core_growth_source_recheck","core_activation_volterra_recheck",
 "prescribed_continuation_foundation","core_modulation_repair_join_recheck",
 "core_global_pressure_moment_recheck","core_coefficient_support_recheck",
 "background_jet_foundation","background_summation",
 "core_borel_prescription_closure","weighted_background_edges")
FLAGS={
 "same_completed_leading_profile_used_by_background_recursion":True,
 "leading_complex_analyticity_restricted_to_protected_collar":True,
 "global_leading_profile_has_all_fixed_smooth_jets":True,
 "global_heat_endpoint_claimed_complex_analytic":False,
 "all_positive_background_coefficients_exist_at_every_fixed_order":True,
 "same_five_row_repair_and_forward_pressure_used_at_each_order":True,
 "rule_defined_borel_sum_identified_with_actual_background":True,
 "separately_fixed_unidentified_numeric_scale_sequence_certified":False,
 "full_pressure_and_curl_cutoff_commutator_retained":True,
 "same_band_phase_normal_value_cone_contract_independently_closed_here":False,
 "canonical_full_path_ODE_or_remaining_primitive_contracts_closed_here":False,
 "unconditional_pulse_residual_force_or_blowup_proved":False}
REMAINING=[
 "independently join actual same-band phase, normal, value and cone gaps on this rule-defined background tail",
 "independently validate the canonical homogeneous full-path ODE with every cutoff and pressure source on that domain",
 "retain fixed compensator/projector supports, compact/slow inverse domains and the complete viscosity/cutoff/interaction residual"]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dependency_review():
 records,sources={},{}
 for name in DEPENDENCIES:
  p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text());records[str(p.relative_to(PROJECT))]=sha(p)
  for rel,digest in r.get('source_hashes',{}).items():
   q=PROJECT/rel;sources[name+':'+rel]=q.exists() and sha(q)==digest
 return records,sources
def generate():
 records,inherited=dependency_review();result=review();protocol=json.loads((HERE/'protocol.json').read_text())
 own=[HERE/n for n in ('README.md','borel_join.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
 gates=dict(exact_join_controls=all(result['checks'].values()),twelve_frozen_dependencies=len(records)==len(DEPENDENCIES),
  all_declared_dependency_sources_match=bool(inherited) and all(inherited.values()),
  scientific_scope_retained=FLAGS['rule_defined_borel_sum_identified_with_actual_background'] and not FLAGS['same_band_phase_normal_value_cone_contract_independently_closed_here'] and not FLAGS['unconditional_pulse_residual_force_or_blowup_proved'])
 return dict(status='same-leading-profile-and-rule-defined-background-borel-contract-joined',protocol=protocol,flags=dict(FLAGS),review=result,gates=gates,
  dependency_evidence_sha256=records,inherited_source_checks=inherited,source_hashes={str(p.relative_to(PROJECT)):sha(p) for p in own},remaining_obligations=list(REMAINING))
def validate(record):
 expected=generate()
 if record!=expected or not all(expected['gates'].values()):raise ValueError('Leading/background/Borel join, provenance or scope changed')
def main(argv=None):
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args(argv)
 if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
 r=generate();validate(r)
 if a.verify_record:validate(json.loads(a.verify_record.read_text()))
 a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
 print('leading/background/Borel join:',len(r['review']['checks']),'exact controls;',len(r['inherited_source_checks']),'dependency source checks')
if __name__=='__main__':main()
