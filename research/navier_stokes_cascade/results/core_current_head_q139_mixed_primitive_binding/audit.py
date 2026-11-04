"""Q139 provenance replay; hashes and diagnostics are not continuum proofs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import q139
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = HERE.parent/'core_current_head_q138_mixed_modulation_pullback'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def generate():
    sys.path.insert(0,str(PARENT))
    spec=importlib.util.spec_from_file_location('q139_parent',PARENT/'audit.py')
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
    paths=[HERE/n for n in ('README.md','q139.py','audit.py','test_audit.py')]
    providers=('core_current_head_q138_mixed_modulation_pullback','core_current_head_q130_mixed_loop_binding','core_current_head_q129_moving_inverse_gevrey','core_current_head_q64_fast_phase_radius','core_current_head_q65_signed_root_radius','core_current_head_q66_cutoff_amplitude')
    paths += [HERE.parent/n/'README.md' for n in providers]
    return {'status':'six-mixed-scalar-jets-to-fixed-origin-mixed-primitives',
        'parent':'0eb03e5b3cd6a87e1a81d9212de71ae58da3e2c0','closure':counts,
        'dependency_evidence_sha256':{'results/core_current_head_q138_mixed_modulation_pullback/evidence.json':sha(PARENT/'evidence.json')},
        'source_hashes':{str(p.relative_to(PROJECT)):sha(p) for p in paths},
        'symbolic':q139.record()}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args();r=generate()
    if a.verify_record:assert r==json.loads(a.verify_record.read_text())
    with a.output.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
    print(r['status'],r['closure'])
