"""Q133 provenance replay; hashes and diagnostics are not continuum proofs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import q133
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = HERE.parent/'core_current_head_q132_exact_repair_gevrey'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def generate():
    sys.path.insert(0,str(PARENT))
    spec=importlib.util.spec_from_file_location('q133_parent',PARENT/'audit.py')
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
    paths=[HERE/n for n in ('README.md','q133.py','audit.py','test_audit.py')]
    paths += [HERE.parent/n/f for n,f in [('core_current_head_q132_exact_repair_gevrey','README.md'),('core_current_head_q132_exact_repair_gevrey','q132.py'),('core_current_head_q131_swirl_memory_gevrey','README.md'),('core_current_head_q127_annular_radius_axis','README.md'),('core_current_head_q127_annular_radius_axis','q127.py'),('core_current_head_q126_exterior_axial_radius','README.md'),('core_current_head_q126_exterior_axial_radius','q126.py'),('core_current_head_q63_repair_anisotropy','README.md'),('core_current_head_q124_full_pressure_discrepancies','README.md'),('core_current_head_q125_uniform_angular_first_seed','README.md')]]
    return {'status':'final-leading-angular-assembly-and-norm-order-boundary',
        'parent':'efea38379f0d75edfe3b66f4a565e7196da0987f','closure':counts,
        'dependency_evidence_sha256':{'results/core_current_head_q132_exact_repair_gevrey/evidence.json':sha(PARENT/'evidence.json')},
        'source_hashes':{str(p.relative_to(PROJECT)):sha(p) for p in paths},
        'symbolic':q133.record()}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args();r=generate()
    if a.verify_record:assert r==json.loads(a.verify_record.read_text())
    with a.output.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
    print(r['status'],r['closure'])
