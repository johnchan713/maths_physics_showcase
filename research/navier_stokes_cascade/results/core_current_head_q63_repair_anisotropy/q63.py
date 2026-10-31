"""Explicit fixed bump and anisotropic exact-repair majorants."""
from math import factorial


def bump(r):
    if r < 0:
        raise ValueError("negative derivative")
    return 162*32600**(r+1)*factorial(r+1)**3


def repair(beta,q,z):
    out=[2*beta*z[0]]
    for k in range(1,len(z)):
        out.append(2*beta*(z[k]+q*sum(out[j]*out[k-j] for j in range(1,k))))
    return out


def record():
    return {"bump_first":[bump(i) for i in range(3)],
            "profile_total_offset":99998,"profile_eta_offset":2,
            "physical_bounds_verified":False}
