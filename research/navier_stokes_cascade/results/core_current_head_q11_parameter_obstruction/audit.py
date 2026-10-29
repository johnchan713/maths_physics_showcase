"""Exact parameter-interface obstruction; hashes establish provenance only."""
import ast
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import argparse

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
RESULTS = HERE.parent
PARENT = 'd1aea3715e41ae2db7440190d6a8f4a69f5ee733'

def cross_exponent(h, stage, m=0, kappa=Q(1, 100000)):
    h, kappa = Q(h), Q(kappa)
    if not 0 < h < Q(1,100) or not 0 < kappa < Q(1,6):
        raise ValueError('positive h and stage gain required')
    if any(type(v) is not int or v < 0 for v in (stage, m)):
        raise ValueError('nonnegative integer stage and order required')
    return h * (Q(1,2)-3*kappa) * (stage+1)/2 - (15+4*h+4*m)

def stage_floor(h, m=0, target=0, kappa=Q(1,100000)):
    cross_exponent(h, 0, m, kappa)
    if type(target) is not int or target < 0:
        raise ValueError('nonnegative integer target required')
    x=2*(target+16+4*m+4*Q(h))/(Q(h)*(Q(1,2)-3*Q(kappa)))-1
    return -(-x.numerator//x.denominator)

def closure():
    seen = set(); edges = 0; sources = 0
    def visit(path):
        nonlocal edges, sources
        if path in seen: return
        seen.add(path)
        record=json.loads(path.read_text())
        for key,digest in record.get('source_hashes',{}).items():
            assert hashlib.sha256((PROJECT/key).read_bytes()).hexdigest()==digest, key
            sources+=1
        for key,digest in record.get('dependency_evidence_sha256',{}).items():
            child=PROJECT/key
            assert hashlib.sha256(child.read_bytes()).hexdigest()==digest, key
            edges+=1; visit(child)
    visit(RESULTS/'core_current_head_q10_pressure_comparison_join/evidence.json')
    return dict(nodes=len(seen), edges=edges, source_bindings=sources)

def review():
    q1=json.loads((RESULTS/'core_current_head_q1_exact_map_join/evidence.json').read_text())
    root=json.loads((RESULTS/'root_foundation/evidence.json').read_text())
    sample=json.loads((RESULTS/'core_current_head_q6_complete_residual_join/evidence.json').read_text())['review']['samples'][0]
    q7=ast.parse((RESULTS/'core_current_head_q7_corner_endpoint_join/corner_endpoint.py').read_text())
    literal=next(n.value for n in q7.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='H' for t in n.targets))
    assert isinstance(literal,ast.Call) and literal.func.id=='Q'
    selected=Q(*(ast.literal_eval(x) for x in literal.args))
    q6=ast.parse((RESULTS/'core_current_head_q6_complete_residual_join/complete_residual.py').read_text())
    calls=[n for n in ast.walk(q6) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and isinstance(n.func.value,ast.Name) and n.func.value.id=='current_residual' and n.func.attr=='comparison']
    assert len(calls)==1 and len(calls[0].args)==5 and not calls[0].keywords
    provider=ast.parse((RESULTS/'core_current_log_reserved_residual_reconciliation/residual.py').read_text())
    fn=next(n for n in provider.body if isinstance(n,ast.FunctionDef) and n.name=='comparison')
    defaults=dict(zip([a.arg for a in fn.args.args][-len(fn.args.defaults):],fn.args.defaults))
    default_h=Q(*(ast.literal_eval(x) for x in defaults['h'].args))
    assert 'h=lambda^2' in q1['protocol']['identity']
    assert 'lambda=exp(-4T)' in root['protocol']['actual_parameters']
    assert root['protocol']['actual_parameters'].startswith('T=exp(64)+10>=64')
    upper=Q(1,513) # exp(512)>1+512 and T>=64
    J=sample['log_reserved_stage']
    diagnostic=cross_exponent(selected,J)
    upper_cross=cross_exponent(upper,J)
    assert default_h==selected==Q(1,200)>upper
    assert diagnostic==Q(sample['six_group_exponents']['old_new'])
    assert upper_cross < 0
    return dict(status='same-object-parameter-interface-obstructed',parent=PARENT,
        constructed_h='exp(-8*(exp(64)+10))',constructed_h_strict_upper=str(upper),
        Q6_default_h=str(default_h),Q7_literal_h=str(selected),stage=J,
        diagnostic_cross_exponent=str(diagnostic),actual_cross_strict_upper=str(upper_cross),
        required_cross_stage='J >= ceil(2*(N+16+4*m+4*h)/(h*(1/2-3*kappa))-1)',
        same_parameter_join_valid=False,continuum_theorem_disproved=False,
        endpoint_bound_proved=False,closure=closure())

def generate():
    result=review()
    paths=[HERE/'README.md',HERE/'audit.py',HERE/'test_audit.py']
    result['source_hashes']={str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-record',type=Path)
    a=p.parse_args();record=generate()
    if a.verify_record: assert record==json.loads(a.verify_record.read_text())
    with a.output.open('x') as f: json.dump(record,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(record,indent=2))
