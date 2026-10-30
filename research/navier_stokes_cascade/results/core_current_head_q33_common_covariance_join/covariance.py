"""Conditional canonical root/weighted covariance composition; no field solver."""
from pathlib import Path
import importlib.util
from fractions import Fraction as Q
HERE=Path(__file__).resolve().parent
def load(name,path):
    sp=importlib.util.spec_from_file_location(name,HERE.parent/path)
    m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
normal=load('q33_normal','core_current_head_q32_stress_normalization/normalization.py')
columns=load('q33_columns','core_canonical_haar_column_jet_closure/bounds.py')
haar=load('q33_haar','core_initial_haar_row_budget/budgets.py')

def canonical_root(template,column,cutoff,stress,geometry):
    """Rebind norm inputs only; preserve canonical value data and all curl jets."""
    joined=columns.root_join(template,column,cutoff)
    roots=columns.entrance.inverse.roots
    target=lambda j: roots.P(normal.transfer(j,stress(j),**geometry)['root_target'].c)
    return roots.RootBudget(H=joined.H,target=target,wave=joined.wave,
        curl=joined.curl,loghalf=joined.loghalf,K0=joined.K0,beta=joined.beta,
        R=joined.R,mass=joined.mass,width=joined.width,trace=joined.trace)

class PrimaryProxy:
    def __init__(self,primary):self.primary=primary
    def amplitude(self,j):return haar.P(self.primary.amplitude(j).c)
    def complete_curl(self,j):return haar.P(self.primary.complete_curl(j).c)

def common_inputs(n,primary,stress,geometry,kappa=Q(1,100000)):
    normal.stress.order(n)
    kappa=Q(kappa)
    if not 0<kappa<Q(1,4):raise ValueError('Retained initial-row kappa range required')
    cov=haar.weighted_covariance(n,PrimaryProxy(primary))
    higher=haar.P(normal.transfer(n,stress(n),**geometry)['higher_stress_raw'])
    return dict(W=cov['W'],D=cov['Dcurl']+higher,Dcurl=cov['Dcurl'],
                higher=higher,pole=12*n+18,orders={'W':Q(1),'D':2-2*kappa,'higher':Q(3)})

def record():
    return dict(result='conditional-common-covariance-and-stress-defect-inputs',
        wave_factor='16*188^n*n!; S_i<=25*S/16',
        weight_preserved=True,full_six_entry_covariance=True,
        selected_defect='Dcurl+4*Cnative_n*S_n',
        pole='12*n+18',background_evaluated_in_common_chart=True,
        duplicated_stress_or_deck_count=False,independent_W_D_callbacks_removed=True,
        canonical_value_and_ODE_inputs_independently_verified=False,
        full_source_or_residual_or_blowup=False)
