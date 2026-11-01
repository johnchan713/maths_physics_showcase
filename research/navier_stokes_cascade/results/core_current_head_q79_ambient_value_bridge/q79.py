"""Exact manufactured budgets; actual giant constants remain symbolic recipes."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
background=load('q79_background',HERE.parent/'core_background_finite_head_budget/budgets.py')
frame=load('q79_frame',HERE.parent/'frame_column_foundation/frame_budget.py')

def exact(v):
    if isinstance(v,(bool,float)) or not isinstance(v,(int,Q)):raise TypeError('Exact rational bound required')
    return Q(v)

def low_inputs(bg,Dbox):
    D=exact(Dbox)
    if D<0:raise ValueError('Nonnegative box radius required')
    one=bg.bound(1);two=bg.leading_bound(2)
    return dict(Db=3*D*two+one['error_FG'],Mz=one['full_FG'],
                leading_two=two,error_one=one['error_FG'],full_one=one['full_FG'])

def head_budget(khead,r,H,L,J,qtheta,qz,qrad,lam,U):
    k,r,H,L,J,qt,qz,qr,lam,U=map(exact,(khead,r,H,L,J,qtheta,qz,qrad,lam,U))
    if min(k,r,H,L,lam)<=0 or H<r or min(J,qt,qz,qr,U)<0:raise ValueError('Invalid bounds')
    nu=1/(k*H);M=2*(qt+qz)*J;K=(H+5)*J
    D=lam*(1+9*U*U/4)
    return dict(floor=nu,N=qr+qt/r+qz+L*M,M=M,K=K,Dref=D,
                log_ambient=L*(K+M/nu+D))

def tail_budget(p,J):
    b=frame.budgets(p);J=exact(J)
    if J<0:raise ValueError('Nonnegative full first jet required')
    nu=p['bmin']/2;N=p['bmax']*(1+3*p['u']/2)+b['Dr']+b['Dt']
    M=b['Drprime']+b['Dtprime'];K=(p['Rmax']+5)*J
    Qlift=K/nu+M/(nu*nu)
    kframe=(1+b['S0'])*max(1,p['cmax'])/min(1,p['cmin'])
    return dict(floor=nu,N=N,M=M,K=K,Qlift=Qlift,Lstar=b['normal_acceptance'],
                T_prefactor=kframe,T_exponent=b['CE']+b['Cd'],
                # At(S)=Tt*(1+S*linear_factor)+normal_term, not exp(S*...).
                linear_factor=Qlift*N,normal_term=N/nu)

def record():
    trace=set();bg=background.fixture(trace);low=low_inputs(bg,Q(1,10))
    p=frame.inputs(Db=low['Db'],Mz=low['Mz'])
    return {'status':'conditional-explicit-ambient-value-bridge',
        'requests':[list(t) for t in sorted(trace)],
        'fixture_low':{k:str(v) for k,v in low.items()},
        'fixture_tail':{k:str(v) for k,v in tail_budget(p,low['full_one']).items()},
        'fixture_head':{k:str(v) for k,v in head_budget(9,1,2,64,3,2,1,1,2,1).items()},
        'khead_recipe':'1+2^ceil(h ell0/2)',
        'normal_reciprocal':'max(nu_h^-2,nu_t^-2)',
        'ambient_polynomial':'Th+Tt+Nt/nu_t+S Tt Qt Nt',
        'actual_finalized_coefficients_certified':False,
        'original_selector_changed':False,'column_root_floor_proved':False,
        'unconditional_actual_pulse_bound':False,'complete_residual_or_blowup':False}
