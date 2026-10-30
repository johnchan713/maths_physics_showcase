"""Q29 provenance replay; hashes and diagnostics are not continuum proofs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import compact
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
Q12 = HERE.parent/'core_current_head_q28_initial_mean_composition'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def generate():
    sys.path.insert(0,str(Q12))
    spec=importlib.util.spec_from_file_location('q29_parent',Q12/'audit.py')
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
    providers=['core_fast_inverse_budget/README.md','core_signed_target_edge_budget/README.md',
        'core_initial_haar_row_budget/README.md','core_compact_slow_inverse_reconstruction/README.md',
        'core_current_head_q4_fixed_profile_join/README.md']
    paths=[HERE/n for n in ('README.md','compact.py','audit.py','test_audit.py')]
    paths += [HERE.parent/n for n in providers]
    return {'status':'conditional-explicit-weighted-shifted-compact-constants',
        'parent':'d218ebe560bc54586a7eac536b7267b981153fd1','closure':counts,
        'dependency_evidence_sha256':{'results/core_current_head_q28_initial_mean_composition/evidence.json':sha(Q12/'evidence.json')},
        'source_hashes':{str(p.relative_to(PROJECT)):sha(p) for p in paths},
        'symbolic':compact.record()}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args();r=generate()
    if a.verify_record:assert r==json.loads(a.verify_record.read_text())
    with a.output.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
    print(r['status'],r['closure'])
