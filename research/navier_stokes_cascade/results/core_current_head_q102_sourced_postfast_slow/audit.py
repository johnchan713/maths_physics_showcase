"""Q102 provenance replay; hashes and diagnostics are not continuum proofs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import q102
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = HERE.parent/'core_current_head_q101_final_rows_fast_join'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def generate():
    sys.path.insert(0,str(PARENT))
    spec=importlib.util.spec_from_file_location('q102_parent',PARENT/'audit.py')
    parent=importlib.util.module_from_spec(spec);spec.loader.exec_module(parent)
    assert parent.generate()==json.loads((PARENT/'evidence.json').read_text())
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
    visit(PARENT/'evidence.json')
    providers=('core_current_head_q99_changed_scalar_rows','core_current_head_q39_current_scalar_refresh','core_actual_scalar_gain_join','core_preserved_haar_scalar_join')
    paths=[HERE/n for n in ('README.md','q102.py','audit.py','test_audit.py')]
    paths += [HERE.parent/n/'README.md' for n in providers]
    paths += [HERE.parent/n/f for n,f in (('core_current_head_q99_changed_scalar_rows','q99.py'),('core_actual_scalar_gain_join','budgets.py'),('core_preserved_haar_scalar_join','budgets.py'))]
    return {'status':'sourced-postfast-correlations-and-fresh-slow-map',
        'parent':'5be85f46c17ba97514553f01f795b538eef763a5','closure':counts,
        'dependency_evidence_sha256':{'results/core_current_head_q101_final_rows_fast_join/evidence.json':sha(PARENT/'evidence.json')},
        'source_hashes':{str(p.relative_to(PROJECT)):sha(p) for p in paths},
        'symbolic':q102.record()}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args();r=generate()
    if a.verify_record:assert r==json.loads(a.verify_record.read_text())
    with a.output.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
    print(r['status'],r['closure'])
