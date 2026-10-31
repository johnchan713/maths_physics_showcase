# Q55: sharp native derivative threshold for the shifted axial tail

[Q54](../core_current_head_q54_axial_row/README.md) reduces the regular
axial row at order 3−5κ to the tail T=A₁ Dz QΦ on the post-signed,
pre-new-fast state. Here we determine what can be obtained from a
single isotropic native derivative norm, without using the equations
of Φ. The earlier sufficient derivative count is valid but one above
the sharp threshold for that class. This is a conditional estimate
for the physical tail and a sharpness result for the **operator**.

## Upper bound without a logarithmic loss

Use Q48's raw through-order sup norm N in (U,Z,Tslow,y), and the
[frozen operator](../core_fast_inverse_budget/README.md) with input
multiplier f₁, output multiplier b=k₁χ_U and radial support length L.
Let H be smooth, zero-Haar and compactly zero-extended in the positive
hull, M≥1, n≥0, and m≥4 integer. Then

    N_n(A₁ H) ≤ C(n,m) L N_n(b) N_(n+s_hi)(f₁) M^(−m/2) N_(n+m)(H),
    s_hi=floor(m/2)+1,  s_lo=ceil(m/2)−1.

Here for example C=2^(2n+3m+8)12^m suffices. The asserted power and
derivative count, not the numerical constant, are sharp.

Differentiate the exact whole-line integral at fixed shift a, so no
output derivative hits M a v. Take the dyadic shells K≤|k|₁<2K,
K=1,2,4,... . On one shell integrate by parts s times in the radial
Fourier phase, where s is one of the two integers above. The traces
vanish; the factor is at most M^(−s)|v·k|^(−s). As in the
[Q48 divisor estimate](../core_current_head_q48_original_axial_tail/README.md)
the inverse divisor square sum is ≤4(12K)^(2s), since s≥1.
On the differentiated amplitude use m−s additional pure torus
derivatives. Parseval split by the larger coordinate bounds its
coefficient square sum **norm** by 2(2/K)^(m−s) times the largest of
these derivative norms. The torus derivatives miss f₁. Cauchy–Schwarz,
integration over the support and the product rule therefore give the
shell factor

    M^(−s) K^(2s−m).

At most n+s derivatives hit f₁ and n+m hit H. The two product expansions
cost at most 2^(2n+s), in addition to the displayed Fourier constants.
For K≤√M take s_hi; for K>√M take s_lo. Dividing by M^(−m/2) leaves
(K/√M)^(2s−m). Its exponent is +2 and −2 for even m, +1 and −1 for
odd m. Each geometric sum costs at most two, independent of where
√M falls between shells. This proves the bound and justifies absolute
summation. Using s=m/2 everywhere for even m would not give this sum.

## A compact input that saturates the scale

Fix a positive native hull, a nonzero output multiplier at (U₀,Z₀)
and a smooth nonnegative interior bump a with integral one. Set v=(1,1−√2).
For A+B√2=(1+√2)^j take k=(A−B,B); ω=v·k=A−√2 B=±1/(A+√2 B).
Thus K=|k|₁=A grows with bounded successive ratios and |ω| is comparable
to 1/K. For each large M choose K comparable to √M, so M|ω|=O(K).
Consider the real part of

    H_M = K^(−n−m) a(U)/f₁(U) exp(2πi[k·y−Mω U]).

The reciprocal is used only on the interior bump. All raw derivatives
through n+m are uniformly bounded, since radial and torus frequencies
are O(K). The shifted phase is constant in the integration variable:
the output equals b K^(−n−m) exp(2πi[k·y−Mω U]). Taking n torus
derivatives in a largest coordinate and the sup over y gives a lower
bound c K^(−m) at the chosen slow point. Thus no positive power beyond
m/2 is available uniformly on this class, even for a fixed smooth
cutoff. Fixed polynomial losses in the band index cannot bridge an
exponential scale gap. These are test inputs, not generated solutions.

For the potential version take the same bump and phase times
sin(K(Z−Z₀)), with amplitude ε K^(−n−d) instead, d≥1. Its N_(n+d)
norm is O(ε), yet at Z₀ the tail of Dz QΦ is bounded below in N_n
by c ε² K^(1−d) (Dz=ε∂Z). All functions have zero Haar.
Slow cutoffs equal to one near the test point may be inserted. This
accounts for the extra Z derivative in a *single* potential norm;
it does not imply a lower bound for the actual Φ.

## Consequence at the current exponent

Here κ=1/100000, S=ℓ², M=c_i^ρ ε^(−κ) with inverse c_i at most Tg S.
The extra c_i power in the estimate costs a polynomial for fixed m.
As in Q48, convert normalized R jets to raw U before the bound, remove
the input weight using its finite high-pole supremum, then restore
the output weight on the cutoff transition and convert back. A bound
on the unweighted native operator does not itself supply these field
or profile constants. For H=Dz QΦ of order two, the conditional rate is
2+mκ/2. The following are minimal integer counts **for this generic
single-norm route**, for output order n:

| Goal | H through | Φ through | Resulting exponent |
|---|---:|---:|---|
| Strictly above 3−6κ | n+199989 | n+199990 | 3−11κ/2 |
| Reach 3−5κ | n+199990 | n+199991 | 3−5κ |

No stage of the pressure is reset: Φ=−2G_R Ψ+γ²+Wzz+pm on Q54's
specified state and Dz acts before A₁. Next bound the near-resonant
radial Fourier integrals of this **particular sum** (including pm), or
supply its uniform weighted high derivatives. The example prevents a
large derivative saving based on compactness and Haar cancellation alone.
It does not address the separate A₀ pressure tail, prove all-scale
pulse estimates, improve the complete residual, or establish blow-up.
