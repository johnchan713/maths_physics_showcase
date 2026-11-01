"""Q88 provenance replay; hashes and diagnostics are not continuum proofs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import q88
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = HERE.parent/'core_current_head_q87_final_angular_adapter'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def generate():
    sys.path.insert(0,str(PARENT))
    spec=importlib.util.spec_from_file_location('q88_parent',PARENT/'audit.py')
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
    providers = ('core_current_head_q84_nonlinear_radial_envelope','core_current_head_q85_covariance_defect_adapter','core_current_head_q86_final_mean_adapter','core_current_head_q43_compact_chart','core_symbolic_fresh_mean_join')
    paths=[HERE/n for n in ('README.md','q88.py','audit.py','test_audit.py')]
    paths += [HERE.parent/n/'README.md' for n in providers]
    paths.append(HERE.parent/'core_current_head_q43_compact_chart/chart.py')
    paths += [HERE.parent/'core_current_head_q87_final_angular_adapter/test_audit.py',HERE.parent/'core_current_head_q84_nonlinear_radial_envelope/q84.py']
    paths += [HERE.parent/'core_fast_row_gain_join/budgets.py',HERE.parent/'core_symbolic_initial_haar_join/bounds.py']
    paths += [HERE.parent/n/f for n,f in (('core_background_finite_head_budget','budgets.py'),('core_current_head_q47_full_pressure_correction','full_rows.py'))]
    return {'status':'same-hull-native-nonlinear-radial-composition',
        'parent':'58bbffdde9087aef859429849d955942ef867802','closure':counts,
        'dependency_evidence_sha256':{'results/core_current_head_q87_final_angular_adapter/evidence.json':sha(PARENT/'evidence.json')},
        'source_hashes':{str(p.relative_to(PROJECT)):sha(p) for p in paths},
        'symbolic':q88.record()}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args();r=generate()
    if a.verify_record:assert r==json.loads(a.verify_record.read_text())
    with a.output.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
    print(r['status'],r['closure'])
