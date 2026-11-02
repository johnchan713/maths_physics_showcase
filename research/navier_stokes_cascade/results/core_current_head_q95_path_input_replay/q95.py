"""Explicit finite whole-path source inputs; no continuum certificates."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q94_full_regular_source/q94.py'
s=importlib.util.spec_from_file_location('q95_source',p);source=importlib.util.module_from_spec(s);s.loader.exec_module(source)
full=source.full;EP=source.EP;convert=source.convert

def quotient(raw,log_jets,log_poles):
    if not raw or len(raw)!=len(log_jets) or len(raw)!=len(log_poles):raise ValueError('Matching nonempty finite jet tables required')
    if any(type(p) is not int or p<0 for p in log_poles):raise ValueError('Exact nonnegative integer poles required')
    a=list(map(convert,raw));ell=list(map(convert,log_jets));H=[EP(1)];hp=[0]
    for j in range(1,len(a)):
        H.append(sum((comb(j-1,i-1)*ell[i]*H[j-i] for i in range(1,j+1)),EP()))
        hp.append(max(log_poles[i]+hp[j-i] for i in range(1,j+1)))
    G=[sum((comb(j,i)*H[i]*a[j-i] for i in range(j+1)),EP()) for j in range(len(a))]
    return dict(jets=G,Bell=H,extra_pole=max(hp),Bell_poles=hp)

def product(a,b):
    if len(a)!=len(b) or not a:raise ValueError('Same finite jet ceiling required')
    a,b=map(lambda z:list(map(convert,z)),(a,b))
    return [sum((comb(j,i)*a[i]*b[j-i] for i in range(j+1)),EP()) for j in range(len(a))]

class PathInputs:
    def __init__(self,current_source):
        if not isinstance(current_source,source.CurrentSource):raise TypeError('Same-state Q94 explicit source required')
        self.source=current_source;self.state=current_source.state
    def normal_row(self,j):
        w=self.source.waves
        native=lambda k:w.convolution(w.inverse,lambda i:w.jet('normal',i),k)
        norm=sum((Q(5**k,factorial(k))*native(k) for k in range(j+1)),EP())
        factor=2*convert(self.state.primary.chart.transfer(j))*16*188**j*factorial(j)
        return 2**j*factorial(j)*factor*source.waves.rescale(norm)
    def build(self,R,linear_table,log_table,log_poles=None):
        full.order(R)
        if len(linear_table)!=R+1 or len(log_table)!=R+1:raise ValueError('Fix the ceiling before supplying all provider tables')
        poles=[0]*(R+1) if log_poles is None else list(log_poles)
        rows=[self.source.coefficient(j,linear_table[j]) for j in range(R+1)]
        raw=[2**j*factorial(j)*row['F'] for j,row in enumerate(rows)]
        q=quotient(raw,log_table,poles);base=max(row['common_pole'] for row in rows)
        b=[self.normal_row(j) for j in range(R+1)];forcing=product(b,q['jets'])
        return dict(raw_path=raw,normalized=q['jets'],normal_row=b,source_normal_pressure=forcing,
                    six_variable_through=sum((Q(6**j,factorial(j))*v for j,v in enumerate(q['jets'])),EP()),
                    source_pole=base,common_pole=base+q['extra_pole'],log_extra_pole=q['extra_pole'],
                    source_ceiling=R,source_factory_reserve=R+2,
                    endpoint_flatness_from_norm_alone=False,actual_native_contracts_certified=False,
                    complete_residual_improvement=False)

def record():
    return dict(status='explicit-first-real-whole-path-quotient-normal-input-replay',
                source_order='1-kappa',normal_pressure_source_order='3/2-kappa',
                source_derivative_shift=0,path_coordinate_count=6,
                pole_rule='max_j 2*p_(j+2) + max_j Bell_pole_j',
                Gaussian_minimum_division=False,endpoint_flatness_from_norm_alone=False,
                inverse_response_constructed=False,actual_native_contracts_certified=False,
                full_residual_or_blowup=False)
