"""Q18 provenance replay; hashes and diagnostics are not continuum proofs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import weights
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
Q12 = HERE.parent/'core_current_head_q17_mean_support_join'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def generate():
    sys.path.insert(0,str(Q12))
    spec=importlib.util.spec_from_file_location('q18_parent',Q12/'audit.py')
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
    providers=['regular_source_foundation/README.md',
               'core_log_renormalized_actual_source_induction/README.md',
               'core_log_renormalized_actual_source_induction/source_induction.py',
               'core_current_head_q5_weighted_inverse_source_join/README.md',
               'core_finite_support_closure/README.md']
    paths=[HERE/n for n in ('README.md','weights.py','audit.py','test_audit.py')]
    paths += [HERE.parent/n for n in providers]
    return {'status':'blanket-full-edge-source-implication-obstructed',
        'parent':'c5752be79a1370765eece93b1b5e91e8355b7767','closure':counts,
        'dependency_evidence_sha256':{'results/core_current_head_q17_mean_support_join/evidence.json':sha(Q12/'evidence.json')},
        'source_hashes':{str(p.relative_to(PROJECT)):sha(p) for p in paths},
        'symbolic':weights.record(),
        'diagnostics':[{'pole':p,'n':4*p+1,'log_ratio_lower':str(weights.radial_log_lower(4*p+1,p))} for p in (0,1,400009,10**30)],
        'scope':{'weight_embedding_obstruction_proved':True,
                 'actual_component_bound_proved':False,'unconditional_residual_or_blowup':False}}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args();r=generate()
    if a.verify_record:assert r==json.loads(a.verify_record.read_text())
    with a.output.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
    print(r['status'],r['closure'])
