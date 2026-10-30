"""Exact rational controls for a real symbolic weighted operator estimate."""
from fractions import Fraction as Q

def constants(m,*,f,g,k,chi,one_minus_chi,chi_U,jacobian,length,edge,kernel,interior):
    if type(m) is not int or m<0:raise ValueError('Nonnegative derivative order required')
    values=(f,g,k,chi,one_minus_chi,chi_U,jacobian,length,edge,kernel,interior)
    if any(type(x) not in (int,Q) for x in values):raise TypeError('Exact rational diagnostic inputs required')
    if min(values)<0 or min(jacobian,length,edge,kernel,interior)<=0:raise ValueError('Positive geometry/kernel constants required')
    f,g,k,chi,one_minus_chi,chi_U,J,L,E,C,Z=map(Q,values)
    return {'KT':2**(3*m)*g*(chi+one_minus_chi)*f*J*C,
            'KA':2**(3*m)*k*chi_U*f*J*L*E*Z}

def collar(length,gap):
    length,gap=Q(length),Q(gap)
    if not 0<gap<=length/2:raise ValueError('Actual two-sided plateau gap required')
    return min(Q(1),length/4,gap)

def archive(poles,provider):
    if any(type(p) is not int or p<0 for p in poles):raise ValueError('Nonnegative poles required')
    rows={p:provider(p) for p in sorted(set(poles))}
    if not rows:raise ValueError('Nonempty archive required')
    return rows,{k:sum((v[k] for v in rows.values()),Q(0)) for k in ('KT','KA')}

def record():
    return {'result':'conditional-explicit-weighted-shifted-compact-constants',
        'KT':'2^(3m) g_m (chi_m+(1-chi)_m) f_m Jmax C_p',
        'KA':'2^(3m) k_m (chi_U)_m f_m Jmax Lx E_p Z_gap',
        'shift_preserved':True,'extra_M_power':0,'input_derivative_loss':0,
        'cutoff_derivatives_required':'m+1','weight_floor_on_actual_support_only':True,
        'archive_uses_sum_not_largest_pole':True,
        'slow_response_and_full_source_closed':False,'unconditional_residual_or_blowup':False}
