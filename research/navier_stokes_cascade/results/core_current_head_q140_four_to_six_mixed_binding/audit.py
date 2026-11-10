"""Q140 provenance replay; software checks are not continuum proofs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import q140

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENT = HERE.parent/'core_current_head_q139_mixed_primitive_binding'

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def generate():
    sys.path.insert(0, str(PARENT))
    spec = importlib.util.spec_from_file_location('q140_parent', PARENT/'audit.py')
    parent = importlib.util.module_from_spec(spec); spec.loader.exec_module(parent)
    assert parent.generate() == json.loads((PARENT/'evidence.json').read_text())
    seen=set(); counts={'nodes':0,'edges':0,'source_bindings':0}
    def visit(path):
        if path in seen: return
        seen.add(path); counts['nodes'] += 1
        data=json.loads(path.read_text())
        for name,digest in data.get('source_hashes',{}).items():
            assert sha(PROJECT/name) == digest, name
            counts['source_bindings'] += 1
        for name,digest in data.get('dependency_evidence_sha256',{}).items():
            assert sha(PROJECT/name) == digest, name
            counts['edges'] += 1; visit(PROJECT/name)
    visit(PARENT/'evidence.json')
    paths=[HERE/n for n in ('README.md','q140.py','audit.py','test_audit.py')]
    providers=('core_current_head_q139_mixed_primitive_binding',
               'core_current_head_q128_real_cutoff_signed_root',
               'core_current_head_q66_cutoff_amplitude',
               'core_current_head_q65_signed_root_radius',
               'core_current_head_q72_short_taper_inputs')
    paths += [HERE.parent/n/'README.md' for n in providers]
    return {'status':'four-PREMOD-mixed-fields-to-six-scalar-binding-and-short-taper-obstruction',
            'parent':'90f910a0eb4ebbbb40bad8b6e831d222c8672539',
            'closure':counts,
            'dependency_evidence_sha256':{
                'results/core_current_head_q139_mixed_primitive_binding/evidence.json':sha(PARENT/'evidence.json')},
            'source_hashes':{str(p.relative_to(PROJECT)):sha(p) for p in paths},
            'symbolic':q140.record()}

if __name__ == '__main__':
    p=argparse.ArgumentParser(); p.add_argument('--output',required=True,type=Path); p.add_argument('--verify-record',type=Path)
    a=p.parse_args(); result=generate()
    if a.verify_record: assert result == json.loads(a.verify_record.read_text())
    with a.output.open('x') as f: json.dump(result,f,indent=2,sort_keys=True); f.write('\n')
    print(result['status'], result['closure'])
