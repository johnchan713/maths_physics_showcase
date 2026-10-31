"""First real-prefix angular row: keep the signed leading operator."""
from fractions import Fraction as Q

def ledger(k=Q(1,100000)):
    k=Q(k)
    if not 0<k<Q(1,10): raise ValueError('Small positive construction exponent required')
    B=1-k
    # Time, base transport, old products and at most one graph in viscosity.
    def other(h):return min(h+1,h+1-k,B+h+1-k,2*h+1-k)
    h0=B; h1=h0+1-2*k
    first_error=min(1+1-2*k,other(h0),2-2*k,2-k)
    last_error=min(first_error+1-2*k,other(h1),3-2*k,3-2*k)
    return dict(original_fast=h0,new_fast=h1,new_fast_error=first_error,
                final_leading=h1+1-2*k,final_error=last_error)

def inverse_constant():
    # Dyadic inverse L^5 Nt^-2: Parseval on the eight pure derivatives.
    return 2*(2**2)*2*(12**5)*2*(2**8)

def necessary_order(final_order,k=Q(1,100000)):
    r=ledger(k)
    return min(Q(final_order),r['final_error'])-2+5*k

def record():
    return dict(rates={key:str(value) for key,value in ledger().items()},
                inverse_aux_order=8,inverse_constant=inverse_constant(),
                principal_source='original complete W_rtheta before first fast',
                physical_lower_bound_or_complete_residual=False)
