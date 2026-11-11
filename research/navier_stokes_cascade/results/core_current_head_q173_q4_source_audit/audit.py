#!/usr/bin/env python3
import argparse, hashlib, json
from fractions import Fraction
from pathlib import Path
import q173

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
PARENTS = tuple(HERE.parent / name / "evidence.json" for name in q173.PROVIDERS)
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def encode(v):
    if isinstance(v,dict): return {k:encode(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)): return [encode(x) for x in v]
    if isinstance(v,Fraction): return str(v)
    return v
def generate():
    seen=set(); counts={"nodes":0,"edges":0,"source_bindings":0}
    def visit(path):
        path=path.resolve()
        if path in seen: return
        seen.add(path); counts["nodes"]+=1
        data=json.loads(path.read_text())
        for name,digest in data.get("source_hashes",{}).items():
            assert sha(PROJECT/name)==digest,name; counts["source_bindings"]+=1
        for name,digest in data.get("dependency_evidence_sha256",{}).items():
            child=PROJECT/name; assert sha(child)==digest,name
            counts["edges"]+=1; visit(child)
    for parent in PARENTS: visit(parent)
    own=[HERE/name for name in ("README.md","q173.py","audit.py","test_audit.py")]
    readmes=[p.parent/"README.md" for p in PARENTS]
    return {
        "status":q173.record()["status"],
        "parent":"e8f403601c85d2942f4cd5bb2d9c127aa5a35f7e",
        "closure":counts,
        "dependency_evidence_sha256":{str(p.relative_to(PROJECT)):sha(p) for p in PARENTS},
        "source_hashes":{str(p.relative_to(PROJECT)):sha(p) for p in own+readmes},
        "review":encode(q173.review()),
        "symbolic":encode(q173.record()),
    }
if __name__=="__main__":
    p=argparse.ArgumentParser(); p.add_argument("--output",required=True,type=Path); p.add_argument("--verify-record",type=Path); a=p.parse_args()
    record=generate()
    if a.verify_record: assert record==json.loads(a.verify_record.read_text())
    with a.output.open("x") as f: json.dump(record,f,indent=2,sort_keys=True); f.write("\n")
    print(record["status"],record["closure"])
