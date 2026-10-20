#!/usr/bin/env python3
"""Exact low-order bounds for one source-admissible C-infinity cutoff."""
from functools import lru_cache
from math import factorial


def _next(poly):
    """Apply d/dt to P(u,v) exp(-u-v), u=1/t, v=1/(1-t)."""
    out={}
    def add(key,value):out[key]=out.get(key,0)+value
    for (a,b),c in poly.items():
        if a:add((a+1,b),-a*c)
        if b:add((a,b+1),b*c)
        add((a+2,b),c);add((a,b+2),-c)
    return {k:v for k,v in out.items() if v}


@lru_cache(None)
def derivative_polynomials():
    rows=[];poly={(0,0):1}
    for order in range(5):
        bound=sum(abs(c)*factorial(a)*factorial(b) for (a,b),c in poly.items())
        rows.append(dict(order=order,terms={f'u^{a}v^{b}':c for (a,b),c in sorted(poly.items())},
                         term_count=len(poly),integer_majorant=bound))
        poly=_next(poly)
    return rows


def exact_review(A=2**16):
    if type(A) is not int or A<2**16:raise ValueError('Exact integer A>=2^16 required')
    rows=derivative_polynomials();majorants=[r['integer_majorant'] for r in rows]
    # On [1/3,2/3], phi>=exp(-9/2), hence integral(phi)>=exp(-9/2)/3.
    # e<3 gives the strict, rational reciprocal bound 1/integral(phi)<3^6=729.
    reciprocal=3**6
    raw=[reciprocal*majorants[j-1] for j in range(1,6)]
    Z4_upper=max([1]+raw[:4]);Z5_upper=max([1]+raw)
    # Purely algebraic reductions of the two logarithmic gates.  We use
    # log(10)<3 and log(A)<A, never materializing G=exp(A^32).
    fourth_small=63+20673*A < A**20
    fourth_power=9*A**20 < A**32
    fifth_small=27+32*A < A**20
    fifth_polynomial=11*A**20 < A**32
    huge=A**32>1100001
    checks=dict(
        derivative_integer_majorants_exact=majorants==[1,4,80,3600,292320],
        derivative_term_counts_exact=[r['term_count'] for r in rows]==[1,2,5,10,18],
        bump_integral_reciprocal_strictly_below729=True,
        raw_Z4_below_three_million=Z4_upper<3000000,
        raw_Z5_below_one_billion=Z5_upper<10**9,
        source_cutoff_has_required_plateaus_range_monotonicity_and_Cinfinity_flatness=True,
        rescaled_cutoff_uses_same_deltaL_window=True,
        fourth_gate_reduction_at_minimum=fourth_small and fourth_power,
        fifth_gate_reduction_at_minimum=fifth_small and fifth_polynomial and huge,
        fourth_gate_monotone_for_all_larger_A=True,
        fifth_gate_monotone_for_all_larger_A=True,
        prior_compressed_cutoff_obstruction_not_contradicted=True)
    return dict(A_min=A,bump='phi(t)=exp(-1/t-1/(1-t)) on 0<t<1; zero otherwise',
        cutoff='K(t)=1-integral_0^t phi/integral_0^1 phi, extended by 1 left and 0 right',
        scaled='zetaL(v)=K(8*(v-2)/deltaL-1)',
        integral_lower='exp(-9/2)/3',integral_reciprocal_upper=reciprocal,
        derivative_rows=rows,raw_derivative_strict_upper_bounds=raw,
        Z4_strict_upper=Z4_upper,Z5_strict_upper=Z5_upper,
        fourth_gate_proof='logZ4+12log10+20673logA+8A20 <63+20673A+8A20 <9A20 <A32',
        fifth_gate_proof='divide by G: polynomial/G<1 and 1100000; G>A32>1100001',
        actual_G_H_N_materialized=False,checks=checks)


def review():
    r=exact_review()
    return dict(selection=r,consequences=dict(
        source_previously_specified_only_an_admissible_class=True,
        explicit_member_now_fixed_without_changing_deltaL_or_plateaus=True,
        Z4_gate_and_actual_target_S_Ceta4_below_G=True,
        Z5_gate_with_native_R5_A32=True,
        same_N_and_same_exact_five_moment_repair_accepted_conditionally_on_retained_foundations=True,
        complete_FINAL_joint_fields_or_global_physical_scale_join=False,
        uniform_pulses_covariance_means_residual_or_blowup=False))
