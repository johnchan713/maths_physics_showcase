"""Q13 provenance replay; hashes and diagnostics are not continuum proofs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import cutoff
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
Q12 = HERE.parent/'core_current_head_q12_symbolic_parameter'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def generate():
    sys.path.insert(0,str(Q12))
    spec=importlib.util.spec_from_file_location('q13_parent',Q12/'audit.py')
    parent=importlib.util.module_from_spec(spec);spec.loader.exec_module(parent)
    assert parent.generate()==json.loads((Q12/'evidence.json').read_text())
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
    visit(Q12/'evidence.json')
    providers=['core_complete_potential_diagonal_join/README.md',
               'core_certified_complete_potential_diagonal_join/README.md',
               'canonical_local_realization/README.md']
    paths=[HERE/n for n in ('README.md','cutoff.py','audit.py','test_audit.py')]
    paths += [HERE.parent/n for n in providers]
    return {'status':'all-order-cutoff-factor-proved-representative-bound-open',
        'parent':'18d8c52172856032ef031ad4302db0b791acbae7',
        'closure':counts,
        'dependency_evidence_sha256':{'results/core_current_head_q12_symbolic_parameter/evidence.json':sha(Q12/'evidence.json')},
        'source_hashes':{str(p.relative_to(PROJECT)):sha(p) for p in paths},
        'diagnostic_bounds':{w:str(cutoff.coefficient_bound(w,cutoff.Q(1,513))) for w in ('z','t','zz','zt','tt','zztt')},
        'scope':{'cutoff_factor_all_order_bound_proved':True,
                 'actual_representative_suprema_proved':False,
                 'numerical_cutoff_sequence_required':False,
                 'cutoff_sequence_reselected':False,
                 'endpoint_bound_proved':False,'unconditional_residual_or_blowup':False}}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args();r=generate()
    if a.verify_record:assert r==json.loads(a.verify_record.read_text())
    with a.output.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
    print(r['status'],r['closure'])
