"""Exact sufficient cone certificates; inputs are enclosures, not samples."""
from fractions import Fraction as Q

def exact(x):
    if isinstance(x,(float,bool)) or not isinstance(x,(int,Q)):raise TypeError('Exact inputs required')
    return Q(x)

def vector(x):
    if len(x)!=2:raise ValueError('Two flux components required')
    return tuple(map(exact,x))

def det(a,b):return a[0]*b[1]-a[1]*b[0]

def enclosure(plus,minus,target,errors,orientation=1):
    p,m,t=map(vector,(plus,minus,target));ep,em,et=map(exact,errors)
    if min(ep,em,et)<0 or orientation not in (-1,1):raise ValueError('Invalid error/orientation')
    up,um,ut=(sum(map(abs,x)) for x in (p,m,t))
    E=(ep*um+em*up+ep*em,et*um+em*ut+et*em,ep*ut+et*up+ep*et)
    refs=(det(p,m),det(t,m),det(p,t))
    return {'margins':tuple(orientation*x-e for x,e in zip(refs,E)),
            'error_bounds':E,'column_upper':max(up+ep,um+em)}

def root_constants(records,rho_min):
    rho=exact(rho_min)
    if rho<=0 or not records:raise ValueError('Positive stress floor and nonempty certified cover required')
    d=min(r['margins'][0] for r in records)
    a=min(min(r['margins'][1:]) for r in records)
    M=max(r['column_upper'] for r in records)
    if min(d,a,M)<=0:raise ValueError('Three strict cone margins required')
    return {'beta':rho*a/M**2,'K0':2*M/d,'d':d,'a':a,'M':M}

def record():
    return {'status':'finite-head-canonical-cone-gate',
        'required_margins':['o det(c_plus,c_minus)>0','o det(tau,c_minus)>0','o det(c_plus,tau)>0'],
        'beta_head':'rho_min a/M^2','K_head':'2M/d',
        'unbounded_target_magnitude_not_compactified':True,
        'actual_head_margins_proved':False,'selector_changed':False,
        'unconditional_residual_or_blowup':False}
