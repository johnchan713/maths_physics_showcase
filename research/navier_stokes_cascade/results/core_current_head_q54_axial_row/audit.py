"""Q54 provenance replay; hashes and diagnostics are not continuum proofs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import axial
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = HERE.parent/'core_current_head_q53_primary_oscillation'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def generate():
    sys.path.insert(0,str(PARENT))
    spec=importlib.util.spec_from_file_location('q54_parent',PARENT/'audit.py')
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
        'core_axial_potential_pressure_gain', 'core_fast_row_gain_join',
        'core_current_head_q47_full_pressure_correction', 'core_current_head_q48_original_axial_tail',
        'core_current_head_q49_initial_flux_tail','core_initial_final_row_obstruction',
        'core_initial_target_rate_gain','core_current_head_q51_radial_curl_target',
        'core_current_head_q52_double_viscosity')]
    paths=[HERE/n for n in ('README.md','axial.py','audit.py','test_audit.py')]
    paths += [HERE.parent/n for n in providers]
    return {'status':'conditional-axial-viscosity-and-postsigned-tail',
        'parent':'d7913c10901da7d8e4aeceaac81cc059ba4d69c5','closure':counts,
        'dependency_evidence_sha256':{'results/core_current_head_q53_primary_oscillation/evidence.json':sha(PARENT/'evidence.json')},
        'source_hashes':{str(p.relative_to(PROJECT)):sha(p) for p in paths},
        'symbolic':axial.record()}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args();r=generate()
    if a.verify_record:assert r==json.loads(a.verify_record.read_text())
    with a.output.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
    print(r['status'],r['closure'])
