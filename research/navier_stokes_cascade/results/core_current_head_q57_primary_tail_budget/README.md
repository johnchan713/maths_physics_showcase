# Q57: primary-amplitude budget for the first-real axial tail

[Q56](../core_current_head_q56_axial_pressure_moment/README.md) isolated
one signed moment with full pressure memory. Here is a sufficient
order-3−5κ estimate on that **actual moment** in terms of native
principal amplitudes, slow multipliers and only low jets of the curl.
It is conditional on their uniform polynomial enclosures. The large
counts below are for the additional moment calculation, not the
entire lower-input construction.

## Charge the complete curl without oscillatory integration

In the coarsest frozen band let P_ab=sum A_ia A_ib, with A_i the real
transferred principal amplitude including the root, temporal/spatial
cutoff and slow partition. In normalized units
W_ab=ε(P_ab+sum k_i^(-2) B_ia B_ib)/2 exactly. The real quadrature and
support separation are essential. The [aggregate atlas](../core_primary_aggregate_atlas/README.md)
bounds k_i^(-2)≤ε and the label-l² jets of ε^κ B by b_j(S). Hence the
curl part of every W entry through order j is at most
2^(j−1) b_j² ε^(2−2κ). Dz adds ε. The direct whole-line bound on the
positive hull (and the four-extra-jet temporal inverse for the two
shears) makes the entire curl contribution O_(3−2κ), at output n
using b through n+5. No smallness of the differentiated carrier is
assumed: its index is frozen. This estimate retains the full curl;
it does not assert its diagonal differences vanish.

## Native integrals and the joint divisor estimate

Write R_U=dR/dU. In dU the two diagonal inputs are
Pzz−Prr, Pθθ−Prr with slow weights u=R_U R, R_U w/R respectively.
The shear inputs Prz, Prθ have u=R_U 2R G_R, R_U(a'−2a/R).
Here a=2Vw/R and w=cχ−R²/2 as in Q56; prime means ∂R. Put r=0
for diagonal and r=1 for shear. Each principal term is

    (ε²/2) B₁ ∫ [ ∂Z(u (c_i Nt)^(-r) QP) ](Uout+x,y+Mx v_r) dx.

The Z derivative acts on the entire product, **before** multiplying
by B₁. Thus the derivative of the pressure cutoff moment cχ is still
present. The temporal inverse acts only in y, commutes with this
integration and all these slow factors, and has no zero-mode input.

With N_n the raw through-n common native norm, estimate the integral
on a dyadic shell K≤|k|₁<2K. Write l=v_r·k,t=v_t·k and
b=√2−1. The exact identity |lt|=b |(k₁+k₂)²−2k₂²|≥b for k≠0
implies |l|^(-s)|t|^(-r)≤b^(-r)|l|^(-(s−r)).
Integrate by parts s times in x and apply the [Q48 divisor spacing](../core_current_head_q48_original_axial_tail/README.md)
to s−r≥1. The inverse-symbol square-sum norm costs at most
2 b^(-r)(12K)^(s−r), omitting favorable 2π factors. By Parseval
split at the larger coordinate, q pure y derivatives give a factor
2(2/K)^q for the coefficient square-sum norm. Cauchy–Schwarz then
bounds the shell by a constant times

    c_i^(-r) M^(-s) K^(s−r−q).

This argument includes all shear modes, including those with small
temporal divisor. A uniform extra factor K^(-1) for Nt^(-1) ALONE
would be false.

Let p=99995 and take s=p+1 for K≤√M and s=p−1 otherwise, with
q=2p−r−s. As in [Q55](../core_current_head_q55_axial_tail_threshold/README.md)
the factor is M^(-p)(K/√M)^(±2), so its dyadic sum costs at most
4M^(-p). The frozen temporal floor c_i^(-r−ρp) has polynomial scale
cost after substituting M=c_i^ρ ε^(-κ). Thus ε² M^(-p)c_i^(-r)
has the target power 3−5κ. Absolute convergence follows from these
summable estimates on smooth zero-extended inputs.

## Explicit sufficient input bounds

For each of the two choices (s,q) and each weight u, bound raw
derivatives of u and the label-l² derivatives of A on the whole
integration hull. For u the needed index set is the downward closure
of α+e_Z+s e_U, |α|≤n. For A it is the downward closure of
α+e_Z+s e_U+q e_yi, |α|≤n, i=1,2. Denote the respective maxima by
D and H. Leibniz and label Cauchy–Schwarz bound derivatives of P by
2^(n+1+s+q) H². No label count or deck multiplicity is required.
For example, with L the common native hull length, a safe raw bound
for each principal term is

    ε² M^(-p)c_i^(-r) L N_n(B₁)
      × 2^(4n+6p+16) 12^(p+1) b^(-r) max_(s,q)(D H²).

Indeed outer product differentiation costs ≤2^n; splitting derivatives
between u and P costs ≤2^(n+1+s); taking a diagonal difference costs
two. The stated constant dominates these factors, the preceding Fourier
constants and both dyadic sums. Q removes only the zero coefficient.
Sum over all four weights. Compact endpoint traces, the fixed-label
chart transfer and the transition weight floor are as in Q48/Q55;
transfer the input and output norms, not derivatives of the weights
used to define those norms. A high input edge pole may be removed by
the finite fixed-profile unweighting of the atlas. The branch envelope
8·188^j πA_j((25/16)S) at j=n+1+2p−r suffices for H on its native
per-label hypotheses.

| Input for moment at n | Sufficient raw order | Restriction |
|---|---:|---|
| Diagonal principal amplitudes | n+199991 | Only n+1 derivatives in Z,T combined |
| Shear principal amplitudes | n+199990 | Temporal inverse already included |
| Four assembled slow weights u | n+99997 | Only n+1 in Z,T combined; includes R_U |
| Curl amplitudes ε^κ B | n+5 | No high-order Fourier estimate |
| Outer multiplier B₁ | n | After the Z derivative |

The slow-weight count refers to the assembled weights; V and G occur
with a radial derivative, whereas cχ depends only on Z,T and needs only
n+1 slow derivatives. These are sufficient, not minimal physical
requirements. Together with Q56's lower-order remainder these estimates
yield the axial tail rate on the stated hypotheses. Next verify the
principal root/ODE and base-weight envelopes at these specific indices,
or bound the summed near-resonant coefficients directly. Neither an
all-scale numerical enclosure, the distinct radial A₀ tail, nor a
complete residual or blow-up theorem follows here.
