"""Conservative phase inverse disk and normalized pullback bounds."""
from fractions import Fraction as Q
from math import factorial


def disk(m,derivative):
    if m<=0 or derivative<=0:
        raise ValueError("positive floor and derivative bound required")
    rho=min(Q(1,2),m/(4*derivative))
    return rho,m*rho/4


def pullback(r,k,N,delta,H):
    if N<1 or delta<=0 or min(r,k)<0:
        raise ValueError("invalid frozen frequency or orders")
    return Q(factorial(r)*factorial(k),N)*sum((2*N/delta)**j*H(r-j,k) for j in range(r+1))


def record():
    return {"disk_control":[str(v) for v in disk(Q(1,8),Q(3))],
            "slow_analyticity_assumed":False,"physical_inputs_enclosed":False}
