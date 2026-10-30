"""Finite-head/tail upper estimates; no change to accepted bands or roots."""
from fractions import Fraction as Q

def exact(x):
    if isinstance(x,(float,bool)) or not isinstance(x,(int,Q)):raise TypeError('Exact rational fixture required')
    return Q(x)

def threshold(h,radius,Lreq):
    h,r,L=map(exact,(h,radius,Lreq))
    if not 0<h<Q(1,100) or not 0<r<=Q(1,12) or L<1:raise ValueError('Invalid fixed data')
    x=max(Q(16),360*(4/h)**5,L/(2*r))
    return -(-x.numerator//x.denominator)

def head_recipe(kmax,Rmax,L,K,Nprime,d,p):
    k,R,L,K,N,d,p=map(exact,(kmax,Rmax,L,K,Nprime,d,p))
    if min(k,R,L)<=0 or min(K,N,d,p)<0:raise ValueError('Positive geometry and nonnegative norms required')
    eta=1/(k*R);B=2*K+N/eta+d
    return {'normal_floor':eta,'normal_reciprocal':eta**-2,
            'log_propagator_constant':L*(B+p),'damping_floor':eta**2}

def combine(head,tail_log,tail_reciprocal,tail_damping):
    log,rec,damp=map(exact,(tail_log,tail_reciprocal,tail_damping))
    if log<0 or min(rec,damp)<=0:raise ValueError('Invalid tail bounds')
    if head is not None:
        log=max(log,head['log_propagator_constant']);rec=max(rec,head['normal_reciprocal']);damp=min(damp,head['damping_floor'])
    return {'log_T':log,'normal_reciprocal':rec,'damping_floor':damp}

def record():
    return {'status':'finite-head-tail-uniform-upper-propagator',
        'analysis_threshold':'ceil(max(16,360(4/h)^5,Lreq/(2r0)))',
        'finite_head_normal_floor':'1/(khead Rmax)',
        'uniform_T':'max(Ttail,Thead)',
        'original_selector_changed':False,'Q22_dominance_needed_for_upper_bound':False,
        'root_or_cone_lower_bound_proved':False,'actual_source_bound_proved':False,
        'unconditional_residual_or_blowup':False}
