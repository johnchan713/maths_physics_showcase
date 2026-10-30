"""Q17 provenance replay; hashes and diagnostics are not continuum proofs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import support
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
Q12 = HERE.parent/'core_current_head_q16_representative_transfer'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def generate():
    sys.path.insert(0,str(Q12))
    spec=importlib.util.spec_from_file_location('q17_parent',Q12/'audit.py')
    parent=importlib.util.module_from_spec(spec);spec.loader.exec_module(parent)
    assert parent.generate()==json.loads((Q12/'evidence.json').read_text())
    seen=set(); counts={'nodes':0,'edges':0,'source_bindings':0}
    def visit(path):
        if path in seen:return
        seen.add(path); counts['nodes']+=1
        data=json.loads(path.read_text())
        for name,digest in data.get('source_hashes',{}).items():
            assert sha(PROJECT/name)==digest,name
            counts['source_bindings']+=1
        for name,digest in data.get('dependency_evidence_sha256',{}).items():
            assert sha(PROJECT/name)==digest,name
            counts['edges']+=1;visit(PROJECT/name)
    visit(Q12/'evidence.json')
    providers=['auxiliary_geometry_audit/README.md','core_primary_aggregate_atlas/README.md',
               'compact_mean_corrections/README.md','compact_mean_foundation/README.md',
               'core_current_head_q4_fixed_profile_join/README.md',
               'core_finite_support_closure/README.md']
    paths=[HERE/n for n in ('README.md','support.py','audit.py','test_audit.py')]
    paths += [HERE.parent/n for n in providers]
    samples=[{'S':S,'radial_length':'7/2','count':support.projected_count(S,support.Q(7,2)),
              'polynomial_majorant':str(support.count_constant(support.Q(7,2))*S**3)}
             for S in (1,256,10**12)]
    return {'status':'projected-mean-pressure-support-count-derived',
        'parent':'b5b4383ec413a9d7d9272c7242ce0d6ead950887','closure':counts,
        'dependency_evidence_sha256':{'results/core_current_head_q16_representative_transfer/evidence.json':sha(Q12/'evidence.json')},
        'source_hashes':{str(p.relative_to(PROJECT)):sha(p) for p in paths},
        'symbolic':support.record(),'rational_diagnostics':samples,
        'scope':{'projected_mean_pressure_count_derived':True,
                 'actual_normalized_coefficient_estimates_proved':False,
                 'mean_labels_disjoint':False,'field_or_width_reselected':False,
                 'endpoint_bound_proved':False,'unconditional_residual_or_blowup':False}}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args();r=generate()
    if a.verify_record:assert r==json.loads(a.verify_record.read_text())
    with a.output.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
    print(r['status'],r['closure'])
