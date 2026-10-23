#!/usr/bin/env python3
"""Freeze compact trace defects; preserve all actual-input limitations."""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from controls import review

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent.parent
DEPENDENCIES=('core_axial_potential_pressure_gain','core_compact_source_quotient',
              'core_angular_pressure_normal_form','core_next_sourced_wave_join','compact_mean_foundation','auxiliary_geometry_audit')
STATUS='compact-projection-pressure-and-physical-A1-mixed-trace-gate'
FLAGS={
    'exact_classical_piece_projection_and_primitive_jump_defects_derived':True,
    'full_current_pressure_memory_and_axial_viscosity_retained':True,
    'all_seam_and_exterior_terms_in_repeated_IBP_retained':True,
    'moving_seam_mixed_trace_chain_rule_and_cutoff_commutators_retained':True,
    'zeroth_output_physical_trace_request_is_deferred_and_explicit':True,
    'finite_regularity_comparison_is_not_canonical_C_infinity_data':True,
    'canonical_support_contract_disproved_or_actual_jump_detected':False,
    'canonical_closed_cover_traces_or_all_scale_norms_verified':False,
    'higher_output_request_graph_or_regular_source_order_improved':False,
    'complete_absolute_residual_decay_infinite_field_or_blowup_proved':False,
    'historical_frozen_evidence_fields_width_or_band_changed':False
}
REMAINING=[
    'prove canonical transported whole-cover value traces of Q Omega and Q K1, including exterior endpoints and filled gaps',
    'prove identical mixed weighted Phi traces or quantify every surviving partial_Z Phi trace at each original physical request, retaining all Fourier and higher-output constants',
    'validate all actual original high-jet/profile/support/Haar inputs on unchanged bands, then advance the signed cycle',
    'complete regular physical decay, controlled infinite field, localization/force regularity and complete velocity growth before blow-up claims'
]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def provenance():
    bound,checks={},{}
    for name in DEPENDENCIES:
        p=HERE.parent/name/'evidence.json';r=json.loads(p.read_text())
        bound[str(p.relative_to(PROJECT))]=sha(p)
        for rel,expected in r['source_hashes'].items():checks[name+':'+rel]=sha(PROJECT/rel)==expected
    return bound,checks
def sources():
    paths=[HERE/n for n in ('README.md','bounds.py','controls.py','audit.py','test_audit.py','protocol.json','requirements.txt','.gitignore')]
    paths += [HERE.parent/n/'README.md' for n in DEPENDENCIES]
    paths += [HERE.parent/'core_axial_potential_pressure_gain/bounds.py',
              HERE.parent/'core_angular_pressure_normal_form/bounds.py',
              HERE.parent/'core_compact_source_quotient/controls.py',HERE.parent/'auxiliary_geometry_audit/geometry.py']
    return {str(p.relative_to(PROJECT)):sha(p) for p in paths}
def canonical(v):
    if isinstance(v,Fraction):return str(v)
    if isinstance(v,dict):return {str(k):canonical(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [canonical(x) for x in v]
    return v
def generate():
    r=canonical(review());bound,checks=provenance()
    return dict(status=STATUS,protocol=json.loads((HERE/'protocol.json').read_text()),flags=deepcopy(FLAGS),review=r,
        gates=dict(all_exact_trace_pressure_and_budget_controls_pass=bool(r['checks']) and all(r['checks'].values()),
            six_historical_records_bound=len(bound)==6,
            all_declared_historical_sources_match=bool(checks) and all(checks.values())),
        dependency_evidence_sha256=bound,inherited_source_checks=checks,source_hashes=sources(),remaining_obligations=REMAINING[:])
def validate(r):
    expected=generate()
    if r!=expected or not all(expected['gates'].values()):raise ValueError('Trace defect, support scope, pressure memory or provenance changed')
def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args(argv)
    if a.output.exists():raise SystemExit('Use a fresh output; frozen evidence is never overwritten')
    r=generate();validate(r)
    if a.verify_record:validate(json.loads(a.verify_record.read_text()))
    with a.output.open('x') as f:f.write(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(f"{STATUS}: {len(r['review']['checks'])} exact controls; {len(r['inherited_source_checks'])} historical source checks")
if __name__=='__main__':main()
