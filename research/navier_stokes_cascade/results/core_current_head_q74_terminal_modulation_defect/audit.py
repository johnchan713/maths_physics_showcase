"""Q74 provenance replay; hashes and diagnostics are not continuum proofs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import q74
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = HERE.parent/'core_current_head_q73_annular_target_all_orders'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def generate():
    sys.path.insert(0,str(PARENT))
    spec=importlib.util.spec_from_file_location('q74_parent',PARENT/'audit.py')
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
    providers = [name+'/README.md' for name in (
        'stress_realization_audit','core_current_head_q63_repair_anisotropy','core_current_head_q64_fast_phase_radius','repair_state_bounds','prescribed_continuation_foundation','core_all_order_foundation','activated_collar_foundation','native_sixth_angular','core_current_head_q63_repair_anisotropy',
        'modulation_repair_all_order','loop_modulation_bounds','stress_realization_audit',
        'explicit_variance_cutoff','core_modulation_repair_join_recheck')]
    paths=[HERE/n for n in ('README.md','q74.py','audit.py','test_audit.py')]
    paths += [HERE.parent/n for n in providers+['axial_stress_audit/diagnostics.py']]
    return {'status':'quadratic-terminal-modulation-and-repair-defect',
        'parent':'b935ee7091c65f4ad5eceb3d95c389630ec253df','closure':counts,
        'dependency_evidence_sha256':{'results/core_current_head_q73_annular_target_all_orders/evidence.json':sha(PARENT/'evidence.json')},
        'source_hashes':{str(p.relative_to(PROJECT)):sha(p) for p in paths},
        'symbolic':q74.record()}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args();r=generate()
    if a.verify_record:assert r==json.loads(a.verify_record.read_text())
    with a.output.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
    print(r['status'],r['closure'])
