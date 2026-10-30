"""Exact norm adapters, not numerical certificates of the underlying fields."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
background=load('q32_background',HERE.parent/'core_background_finite_head_budget/budgets.py')
stress=load('q32_stress',HERE.parent/'core_current_head_q31_stress_input_budget/stress_inputs.py')
roots=load('q32_roots',HERE.parent/'core_primary_root_budget/budgets.py')

def physical_cost(m,Rmin,Rmax):
    stress.order(m)
    c=background.chart_constants(m,Rmin,Rmax)
    exponential=sum((2*c['log'])**k/Q(factorial(k)) for k in range(m+1))
    return 16*factorial(m)*exponential*c['composition']

def radial_cost(m,d,Umin,Rmax):
    stress.order(m)
    d,Umin,Rmax=map(stress.nonnegative,(d,Umin,Rmax))
    if min(d,Umin,Rmax)<=0:raise ValueError('Positive fixed radial chart data required')
    falling=Q(1);increment=Q(0)
    for j in range(1,m+1):
        falling*=1/d-(j-1)
        increment+=abs(falling)*Rmax/Umin**j/Q(factorial(j))
    increment+=4 if m else 0 # four identity coordinates; padding if unused
    return factorial(m)*sum(increment**k/Q(factorial(k)) for k in range(m+1))

def transfer(m,budget,Rmin,Rmax,radial=None):
    stress.order(m)
    if budget['pole']!=3*m+3:raise ValueError('Q31 order/pole mismatch')
    c=physical_cost(m,Rmin,Rmax)
    if radial is not None:c*=radial_cost(m,**radial)
    target=c*stress.nonnegative(budget['target'])
    return dict(target_raw=target,higher_stress_raw=4*c*stress.nonnegative(budget['correction']),
                root_target=roots.P(Q(5**m,factorial(m))*target),pole=3*m+3,
                higher_epsilon_order=3,full_epsilon_order=1)

def record():
    b=transfer(2,stress.fixture(2),1,2,dict(d=2,Umin=1,Rmax=2))
    return dict(result='physical-and-native-stress-normalization',
                target_raw=str(b['target_raw']),higher_stress_raw=str(b['higher_stress_raw']),
                leading_target='s^(-1-h)*T0',normalized_stress='epsilon*s^(-1-h)*Tsum',
                higher_epsilon_order=3,pole=b['pole'],derivatives_of_scale_factor_retained=True,
                pressure_memory_removed=False,all_source_premises_closed=False,
                complete_residual_or_blowup=False)
