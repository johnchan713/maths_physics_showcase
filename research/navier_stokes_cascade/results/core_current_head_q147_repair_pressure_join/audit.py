"""Q147 frozen provenance replay; diagnostics are not continuum proofs."""
import argparse,hashlib,json
from pathlib import Path
import q147
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent.parent
PARENT=HERE.parent/'core_current_head_q146_cap_repair_identity'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def generate():
 seen=set();counts={'nodes':0,'edges':0,'source_bindings':0}
 def visit(path):
  if path in seen:return
  seen.add(path);counts['nodes']+=1;data=json.loads(path.read_text())
  for name,digest in data.get('source_hashes',{}).items():assert sha(PROJECT/name)==digest,name;counts['source_bindings']+=1
  for name,digest in data.get('dependency_evidence_sha256',{}).items():assert sha(PROJECT/name)==digest,name;counts['edges']+=1;visit(PROJECT/name)
 visit(PARENT/'evidence.json')
 paths=[HERE/n for n in ('README.md','q147.py','audit.py','test_audit.py')]
 providers=('core_current_head_q146_cap_repair_identity','core_current_head_q74_terminal_modulation_defect',
  'core_current_head_q131_swirl_memory_gevrey','core_current_head_q132_exact_repair_gevrey',
  'core_current_head_q137_physical_repair_mixed','core_current_head_q138_mixed_modulation_pullback',
  'core_continuation_matching_recheck','matching_root_foundation','repair_state_bounds',
  'regional_cone_join','compact_global_reconciliation','core_modulation_repair_join_recheck')
 paths += [HERE.parent/n/'README.md' for n in providers]
 return {'status':'same-modulation-five-row-repair-cone-and-global-pressure-identity-joined','parent':'37aa64ffed83f55541bcb5f04eab7c806d419db1','closure':counts,
  'dependency_evidence_sha256':{'results/core_current_head_q146_cap_repair_identity/evidence.json':sha(PARENT/'evidence.json')},
  'source_hashes':{str(p.relative_to(PROJECT)):sha(p) for p in paths},'symbolic':q147.record()}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--verify-record',type=Path);a=p.parse_args();r=generate()
 if a.verify_record:assert r==json.loads(a.verify_record.read_text())
 with a.output.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
 print(r['status'],r['closure'])
