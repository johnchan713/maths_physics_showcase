#!/usr/bin/env python3
"""Exact controls for the leading/coefficient/Borel continuum join."""
from fractions import Fraction as F


NODES=("core","activation","continued_leading","repaired_leading",
       "coefficient_recursion","borel_prescription","summed_background")
EDGES=(("core","activation"),("activation","continued_leading"),
       ("continued_leading","repaired_leading"),
       ("repaired_leading","coefficient_recursion"),
       ("coefficient_recursion","borel_prescription"),
       ("borel_prescription","summed_background"))


def topological_order(nodes=NODES,edges=EDGES):
    incoming={n:0 for n in nodes};out={n:[] for n in nodes}
    if len(incoming)!=len(nodes):raise ValueError("Duplicate node")
    for a,b in edges:
        if a not in incoming or b not in incoming:raise ValueError("Unknown node")
        incoming[b]+=1;out[a].append(b)
    ready=sorted(n for n,d in incoming.items() if d==0);order=[]
    while ready:
        a=ready.pop(0);order.append(a)
        for b in sorted(out[a]):
            incoming[b]-=1
            if incoming[b]==0:ready.append(b);ready.sort()
    if len(order)!=len(nodes):raise ValueError("Cycle")
    return order


def stage_margins(n,h,ell,log_constant,previous,boundary,b):
    """Check a supplied real-log Borel stage; never select a scale."""
    if type(n) is not int or n<1:raise ValueError("Positive integer order required")
    h,ell,log_constant,previous,boundary,b=map(F,(h,ell,log_constant,previous,boundary,b))
    if min(h,ell,log_constant,previous,boundary,b)<0 or h==0 or ell==0:
        raise ValueError("Positive gain/log2 and nonnegative logarithms required")
    g=n*h
    return {
        "normalized":g*b-log_constant-n*ell,
        "separation":b-previous-ell,
        "boundary":b-boundary-ell,
    }


def normalized_increment_bound(n,h,ell,log_constant,b):
    """C q^(nh) on q<=exp(-b), represented in exact log variables."""
    n=int(n);h,ell,log_constant,b=map(F,(h,ell,log_constant,b))
    if n<1 or min(h,ell,log_constant,b)<0 or h==0:raise ValueError("Invalid bound data")
    margin=n*h*b-log_constant-n*ell
    return {"log_margin_over_2_minus_n":margin,"bounded_by_2_minus_n":margin>=0}


def geometric_tail(start,t):
    """Sum_{n>=start} 2^-n t^n, for 0<=t<=1."""
    start=int(start);t=F(t)
    if start<1 or not 0<=t<=1:raise ValueError("Invalid geometric tail")
    return (t/F(2))**start/(1-t/F(2))


def support_count(q,scales):
    """Number of nonzero cutoff candidates a_n q<1 for supplied scales."""
    q=F(q);scales=tuple(map(F,scales))
    if q<=0 or any(a<=0 for a in scales):raise ValueError("Positive data required")
    if any(scales[i+1] < 2*scales[i] for i in range(len(scales)-1)):
        raise ValueError("Scale separation missing")
    return sum(a*q<1 for a in scales)


def review():
    order=topological_order();cycle=False
    try:topological_order(edges=EDGES+(("summed_background","core"),))
    except ValueError:cycle=True
    margins=stage_margins(3,F(1,5),F(2,3),4,7,5,20)
    tail=geometric_tail(4,1)
    return {
        "route":order,
        "controls":{
            "stage_margins":{k:str(v) for k,v in margins.items()},
            "worst_geometric_tail":str(tail),
            "support_count":support_count(F(1,20),(1,2,4,8,16,32)),
        },
        "checks":{
            "leading_to_summed_background_route_is_acyclic":len(order)==len(NODES),
            "downstream_to_core_cycle_is_rejected":cycle,
            "supplied_stage_has_strict_rule_margins":all(v>0 for v in margins.values()),
            "all_fixed_jet_tail_is_geometric":tail==F(1,8),
            "separated_scales_give_finite_support_at_fixed_q":support_count(F(1,20),(1,2,4,8,16,32))==5,
            "normalized_stage_bound_retains_nh_gain":normalized_increment_bound(3,F(1,5),F(2,3),4,20)["bounded_by_2_minus_n"],
        },
    }
