"""Raw jet majorants for multiplication/division by a positive envelope."""
from fractions import Fraction as Q
from math import comb

def nonnegative(values):
    result=[]
    for x in values:
        if isinstance(x,(float,bool)) or Q(x)<0:raise ValueError('exact nonnegative data required')
        result.append(Q(x))
    if not result:raise ValueError('nonempty data required')
    return result

def bell(log_jets):
    ell=nonnegative(log_jets);H=[Q(1)]
    # ell[0] is unused: no bound on |log P| or 1/min(P) is needed.
    for j in range(1,len(ell)):
        H.append(sum(comb(j-1,i-1)*ell[i]*H[j-i] for i in range(1,j+1)))
    return H

def convert(weighted_jets,log_jets):
    f=nonnegative(weighted_jets);H=bell(log_jets)
    if len(f)!=len(H):raise ValueError('matching finite jet orders required')
    return [sum(comb(j,i)*H[i]*f[j-i] for i in range(j+1)) for j in range(len(f))]

def log_envelope_jets(length,p_jets):
    p=nonnegative(p_jets);length=nonnegative([length])[0]
    return [Q(0)]+[length*p[j]/2+p[j-1] for j in range(1,len(p))]

def record():
    return {'recurrence':'H_j=sum_i=1^j binom(j-1,i-1) ell_i H_(j-i)',
            'quotient_bound':'Q_j=sum_i=0^j binom(j,i) H_i F_(j-i)',
            'log_jet_bound':'ell_j=(L/2)*p_j+p_(j-1)',
            'zeroth_log_or_inverse_Gaussian_required':False,
            'extra_source_derivative_orders':0,
            'conditional_whole_path_jet_equivalence':True,
            'physical_trace_bound_implies_lifted_bound':False,
            'actual_all_stage_source_enclosures_proved':False,
            'unconditional_residual_or_blowup':False}
