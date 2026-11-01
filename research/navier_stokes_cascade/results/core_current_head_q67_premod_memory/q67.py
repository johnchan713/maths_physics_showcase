"""Normalized forward-memory bookkeeping, not sampled physical enclosures."""
from fractions import Fraction as Q

def primitive_bound(i,j,density,trace,length):
    if i<0 or j<0 or length<0: raise ValueError("invalid index or length")
    return trace[j]+length*density[0,j] if i==0 else Q(density[i-1,j],i)

def inputs(m,l):
    if not 0<=l<=m: raise ValueError("invalid orders")
    return {'field_total':m+1,'field_eta':l+1,'trace_eta':l+1}

def record():
    return {'state_input_orders':inputs(99998,2),
            'actual_regional_inputs_enclosed':False,'incoming_trace_required':True}
