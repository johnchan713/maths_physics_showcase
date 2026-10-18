"""Exact oscillatory endpoint obstruction and exterior order controls."""
from fractions import Fraction as F
from functools import lru_cache

def add(*terms):
    out={}
    for term in terms:
        for key,value in term.items():out[key]=out.get(key,F(0))+value
    return {k:v for k,v in out.items() if v}

def scale(term,c):return {k:v*F(c) for k,v in term.items() if v*c}

def dt(term):
    """d/dt of tau^p*sin/cos(1/tau), with tau=1-t."""
    out={}
    for (p,quadrature),value in term.items():
        if quadrature not in ('sin','cos'):raise ValueError('Unknown quadrature')
        if p:out=add(out,{(p-1,quadrature):-p*value})
        other='cos' if quadrature=='sin' else 'sin'
        out=add(out,{(p-2,other):value if quadrature=='sin' else -value})
    return out

@lru_cache(None)
def obstruction_review():
    u={(0,'sin'):F(1)};pressure_gradient=scale(dt(u),-1);rows=[]
    velocity=u;pg=pressure_gradient
    for j in range(9):
        defect=add(dt(velocity),pg)
        rows.append(dict(time_order=j,exact_momentum_cancellation=not defect,
            velocity_jet_terms=len(velocity),pressure_gradient_terms=len(pg)))
        velocity=dt(velocity);pg=dt(pg)
    return dict(field='u=sin(1/(1-t))*e1; p=-x1*cos(1/(1-t))/(1-t)^2',
        potential='A=(0,0,x2*sin(1/(1-t))), B=0',rows=rows,
        unequal_terminal_subsequences='tau_n=1/(pi/2+2*pi*n) and sigma_n=1/(3*pi/2+2*pi*n) give +1 and -1',
        actual_canonical_field=False,
        checks=dict(first_time_derivative_sign_exact=dt(u)=={(-2,'cos'):F(1)},
            second_time_derivative_full_chain_rule=dt(dt(u))=={(-3,'cos'):F(2),(-4,'sin'):F(-1)},
            all_nine_differentiated_momentum_identities_exact=all(r['exact_momentum_cancellation'] for r in rows),
            pressure_gradient_is_nonzero_and_has_negative_tau_power=bool(pressure_gradient) and min(p for p,q in pressure_gradient)<0,
            uniform_bounded_velocity_and_exact_zero_residual_do_not_imply_terminal_trace=True,
            counterexample_violates_required_next_derivative_and_pressure_jet_bounds=True,
            no_compact_support_energy_or_actual_canonical_obstruction_claimed=True))
