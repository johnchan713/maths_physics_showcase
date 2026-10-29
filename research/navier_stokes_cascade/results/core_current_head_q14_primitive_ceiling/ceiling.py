"""Expand Q12's compressed interface counts before requesting primitive domains."""
from pathlib import Path
import importlib.util
from fractions import Fraction as Q
HERE=Path(__file__).resolve().parent

def load(name, path):
    spec=importlib.util.spec_from_file_location(name,HERE.parent/path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

q12=load('q14_parameter','core_current_head_q12_symbolic_parameter/parameter.py')
first=load('q14_first_real','core_first_real_high_jet_domain_closure/high_jet_join.py')
OFFSET=400009

def expand(interface):
    q12.order(interface)
    result=first.normalized_join(interface)['combined_graph_ceiling']
    assert result==interface+OFFSET
    return result

def request(h,stage,m=0,target=0):
    old=q12.compact_request(h,stage,m,target)
    state=old['state_primitive_ceiling']
    physical=old['enlarged_physical_primitive_ceiling']
    return {'state_interface':state,'enlarged_physical_interface':physical,
            'state_primitive':expand(state),'enlarged_physical_primitive':expand(physical),
            'combined_interface':max(state,physical),'combined_primitive':expand(max(state,physical)),
            'raw_after_extra':old['raw_after_extra'],'base_strict_margins':old['base_strict_margins']}

def pole(primitive,original_pole=1,retained=30):
    q12.order(primitive)
    p,r=Q(original_pole),Q(retained)
    if min(p,r)<0:raise ValueError('nonnegative poles required')
    return max(r,6*p,6*primitive+9,6*primitive+12*p+6)

def archive(requests,original_pole=1,retained=30):
    """Preserve all distinct requested domains; no monotonicity of norm constants."""
    entries=[];distinct=[]
    for value in requests:
        retained=pole(value,original_pole,retained)
        entries.append(retained)
        if retained not in distinct:distinct.append(retained)
    return {'by_request':entries,'distinct':distinct,
            'operators':[{'KT1':6*q,'KA1':6*q,'pressure':(6*q,12*q),
                          'moment':(6*q,12*q),'oriented':6*q} for q in distinct]}
