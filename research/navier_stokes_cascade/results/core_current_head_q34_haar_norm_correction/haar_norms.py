"""Complete the raw-U covariance to through-order R-Haar norm interface."""
from pathlib import Path
import importlib.util
from fractions import Fraction as Q
from math import factorial
HERE=Path(__file__).resolve().parent
def load(name,path):
    sp=importlib.util.spec_from_file_location(name,HERE.parent/path)
    m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
common=load('q34_common','core_current_head_q33_common_covariance_join/covariance.py')
reverse=load('q34_reverse','core_weighted_reverse_chart/bounds.py')
P=common.haar.P

def count(n):
    common.normal.stress.order(n)
    return sum((Q(5**j,factorial(j)) for j in range(n+1)),Q(0))

def r_inputs(n,primary,stress,*,Rmin,Rmax,kappa=Q(1,100000)):
    # Q33 wave covariance is raw in U; background is evaluated directly in R.
    raw=common.common_inputs(n,primary,stress,dict(Rmin=Rmin,Rmax=Rmax),kappa)
    chart=reverse.ReverseChart(Rmin,Rmax)
    t=P(chart.transfer(n).c)
    higher=count(n)*raw['higher']
    return dict(W=t*raw['W'],Dcurl=t*raw['Dcurl'],higher=higher,
                D=t*raw['Dcurl']+higher,pole=raw['pole'],orders=raw['orders'],
                covariance_transfer=t,background_count=count(n))

def initial_scalars(n,inputs,moments,edge):
    """Original zero-mean primary identities, after whole-hull edge absorption.

    inputs(j) is the R-chart through-order weighted output above. edge is
    ONE constant polynomial dominating the common finite pole for n+1.
    """
    common.normal.stress.order(n)
    if set(moments)!=set((-1,0,1,2)):raise ValueError('Complete whole-hull moment constants required')
    if not isinstance(edge,P):raise TypeError('Positive common-pole absorption polynomial required')
    M={k:common.normal.stress.nonnegative(v) for k,v in moments.items()}
    w=edge*inputs(n)['W'];wn=edge*inputs(n+1)['W']
    return dict(P=2*M[-1]*w+(n+1)*M[0]*wn,
                Jtheta=M[2]*w,Jz=2*M[1]*w+Q(n+1,2)*M[2]*wn)

def record():
    return dict(result='R-chart-through-order-Haar-interface-corrected',
        correction='Q33 raw polynomial interface alone does not certify through-order R norms',
        raw_to_normalized_count=[str(count(n)) for n in range(5)],
        covariance_reverse_transfer=True,background_reverse_transfer=False,
        background_direct_R_count=True,common_weight_preserved=True,
        scalar_identities_initial_only=True,full_source_residual_or_blowup=False)
