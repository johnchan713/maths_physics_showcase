"""Q16 provenance replay; hashes and diagnostics are not continuum proofs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import transfer
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
Q12 = HERE.parent/'core_current_head_q15_flat_edge_transfer'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def generate():
    sys.path.insert(0,str(Q12))
    spec=importlib.util.spec_from_file_location('q16_parent',Q12/'audit.py')
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
    providers=['physical_pulse_transfer/README.md','core_physical_tail_join/README.md',
               'core_complete_potential_diagonal_join/diagonal.py','full_potential_realization/README.md']
    paths=[HERE/n for n in ('README.md','transfer.py','audit.py','test_audit.py')]
    paths += [HERE.parent/n for n in providers]
    samples=[]
    for h,j,m in ((transfer.Q(1,513),1,0),(transfer.Q(1,10**100),10**30,7)):
        rows=transfer.rows(h,j,m)
        assert all(r['margin']>0 for r in rows.values())
        samples.append({'h':str(h),'j':j,'m':m,'rows':{name:{k:str(v) for k,v in r.items()} for name,r in rows.items()}})
    return {'status':'conditional-complete-representative-physical-transfer',
        'parent':'0efd95500b4291224821036370c63acd78c4a3a7','closure':counts,
        'dependency_evidence_sha256':{'results/core_current_head_q15_flat_edge_transfer/evidence.json':sha(Q12/'evidence.json')},
        'source_hashes':{str(p.relative_to(PROJECT)):sha(p) for p in paths},
        'symbolic':transfer.symbolic(),'rational_diagnostics':samples,
        'scope':{'conditional_physical_transfer_proved':True,
                 'actual_complete_normalized_inputs_proved':False,
                 'mean_pressure_overlap_verified':False,'field_or_width_reselected':False,
                 'endpoint_bound_proved':False,'unconditional_residual_or_blowup':False}}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args();r=generate()
    if a.verify_record:assert r==json.loads(a.verify_record.read_text())
    with a.output.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
    print(r['status'],r['closure'])
