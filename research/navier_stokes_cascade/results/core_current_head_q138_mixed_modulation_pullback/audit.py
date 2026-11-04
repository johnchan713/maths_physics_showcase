"""Q138 provenance replay; hashes and diagnostics are not continuum proofs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import q138
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = HERE.parent/'core_current_head_q137_physical_repair_mixed'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def generate():
    sys.path.insert(0,str(PARENT))
    spec=importlib.util.spec_from_file_location('q138_parent',PARENT/'audit.py')
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
    paths=[HERE/n for n in ('README.md','q138.py','audit.py','test_audit.py')]
    providers=('core_current_head_q137_physical_repair_mixed','core_current_head_q131_swirl_memory_gevrey','core_current_head_q129_moving_inverse_gevrey','core_current_head_q130_mixed_loop_binding','core_current_head_q64_fast_phase_radius','core_current_head_q121_all_order_radial_tuple')
    paths += [HERE.parent/n/'README.md' for n in providers]
    return {'status':'mixed-Ny-pullback-exponential-and-physical-modulation-envelope',
        'parent':'bde05bc8d97198b7bb664b954c0f37845710ea77','closure':counts,
        'dependency_evidence_sha256':{'results/core_current_head_q137_physical_repair_mixed/evidence.json':sha(PARENT/'evidence.json')},
        'source_hashes':{str(p.relative_to(PROJECT)):sha(p) for p in paths},
        'symbolic':q138.record()}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args();r=generate()
    if a.verify_record:assert r==json.loads(a.verify_record.read_text())
    with a.output.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
    print(r['status'],r['closure'])
