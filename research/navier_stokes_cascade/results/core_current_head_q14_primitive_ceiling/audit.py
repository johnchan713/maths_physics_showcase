"""Q14 provenance replay; hashes and diagnostics are not continuum proofs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import ceiling
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
Q12 = HERE.parent/'core_current_head_q13_cutoff_uniformity'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def generate():
    sys.path.insert(0,str(Q12))
    spec=importlib.util.spec_from_file_location('q14_parent',Q12/'audit.py')
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
    providers=['core_certified_finite_stage_domain_join/certified_recurrence.py',
               'core_first_real_high_jet_domain_closure/high_jet_join.py',
               'core_certified_finite_stage_domain_join/README.md']
    paths=[HERE/n for n in ('README.md','ceiling.py','audit.py','test_audit.py')]
    paths += [HERE.parent/n for n in providers]
    return {'status':'compressed-interface-to-primitive-correction',
        'parent':'f71c542e46948caebbf44d62981451ed615d8c75',
        'closure':counts,
        'dependency_evidence_sha256':{'results/core_current_head_q13_cutoff_uniformity/evidence.json':sha(Q12/'evidence.json')},
        'source_hashes':{str(p.relative_to(PROJECT)):sha(p) for p in paths},
        'primitive_offset':ceiling.OFFSET,
        'state_samples':{str(j):ceiling.expand(ceiling.q12.state_ceiling(j)) for j in (0,1,2,10**30)},
        'symbolic_root_state_primitive':'m+900047+99994*J^2+400035*J',
        'symbolic_root_combined_primitive':'max(I_J(m),I_phys)+400009',
        'maximum_common_pole_increment':2400054,
        'scope':{'fully_expanded_requests_corrected':True,
                 'representative_suprema_proved':False,'physical_field_changed':False,
                 'endpoint_bound_proved':False,'unconditional_residual_or_blowup':False}}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path)
    a=p.parse_args();r=generate()
    if a.verify_record:assert r==json.loads(a.verify_record.read_text())
    with a.output.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
    print(r['status'],r['closure'])
