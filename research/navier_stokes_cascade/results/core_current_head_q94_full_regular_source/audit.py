"""Q94 provenance replay; hashes and diagnostics are not continuum proofs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import q94
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = HERE.parent/'core_current_head_q93_current_weighted_wave'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def generate():
    sys.path.insert(0,str(PARENT))
    spec=importlib.util.spec_from_file_location('q94_parent',PARENT/'audit.py')
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
    providers=('regular_source_foundation','core_current_head_q91_weighted_mean_wave','core_current_head_q92_full_linear_envelope','core_current_head_q90_source_entrance_normal')
    paths=[HERE/n for n in ('README.md','q94.py','audit.py','test_audit.py')]
    paths += [HERE.parent/n/'README.md' for n in providers]
    paths += [HERE.parent/n/f for n,f in (('regular_source_foundation','source_checks.py'),('core_current_head_q86_final_mean_adapter','q86.py'))]
    return {'status':'same-state-complete-first-real-regular-source-join',
        'parent':'d2374fc6aa9bd55f94e2df411e681b2e4ff09973','closure':counts,
        'dependency_evidence_sha256':{'results/core_current_head_q93_current_weighted_wave/evidence.json':sha(PARENT/'evidence.json')},
        'source_hashes':{str(p.relative_to(PROJECT)):sha(p) for p in paths},
        'symbolic':q94.record()}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args();r=generate()
    if a.verify_record:assert r==json.loads(a.verify_record.read_text())
    with a.output.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
    print(r['status'],r['closure'])
