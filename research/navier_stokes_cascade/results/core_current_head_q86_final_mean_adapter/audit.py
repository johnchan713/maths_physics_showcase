"""Q86 provenance replay; hashes and diagnostics are not continuum proofs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import q86
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = HERE.parent/'core_current_head_q85_covariance_defect_adapter'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def generate():
    sys.path.insert(0,str(PARENT))
    spec=importlib.util.spec_from_file_location('q86_parent',PARENT/'audit.py')
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
    providers = ('core_current_head_q47_full_pressure_correction','core_current_head_q39_current_scalar_refresh','core_current_head_q41_postslow_reset','core_current_head_q52_double_viscosity','core_fast_row_gain_join')
    paths=[HERE/n for n in ('README.md','q86.py','audit.py','test_audit.py')]
    paths += [HERE.parent/n/'README.md' for n in providers]
    paths += [HERE.parent/'core_current_head_q39_current_scalar_refresh/refresh_scalars.py', HERE.parent/'core_current_head_q41_postslow_reset/postslow.py']
    paths += [HERE.parent/n/f for n,f in (('core_background_finite_head_budget','budgets.py'),('core_current_head_q47_full_pressure_correction','full_rows.py'))]
    return {'status':'corrected-full-pressure-final-mean-adapter',
        'parent':'d35af6a5d400d40878908e1a5607db0cbb2e6487','closure':counts,
        'dependency_evidence_sha256':{'results/core_current_head_q85_covariance_defect_adapter/evidence.json':sha(PARENT/'evidence.json')},
        'source_hashes':{str(p.relative_to(PROJECT)):sha(p) for p in paths},
        'symbolic':q86.record()}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args();r=generate()
    if a.verify_record:assert r==json.loads(a.verify_record.read_text())
    with a.output.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
    print(r['status'],r['closure'])
