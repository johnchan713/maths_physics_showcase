# Direct leading/background seminorm reconstruction

Parent: `bd83faf6521305b8cd8d9c02c58dd7bf11108588`.
Status: `fixed-order-background-seminorm-and-real-borel-selector-reconstructed-relative-to-leading-certificate`.

The parent identifies five continuum contracts still imported by the
downstream residual, force, energy and comparison chain. Its first contract
combines two logically different assertions: regularity of the finalized
leading profile, and construction plus Borel summation of every positive-order
background coefficient. This checkpoint separates them and reconstructs the
second assertion directly.

The result is analytical. Exact controls record the combinatorics and
quantifiers, but do not evaluate the actual high-order fields or prove a
continuum supremum by testing.

## 1. Fixed coefficient before Borel scale selection

For positive order `n`, with `xi=sqrt(X)`, the zero-axis unknown is

\[
 W_n=(F_n,U_n,M_n/X-U_n,\Pi_n,\partial_\xi F_n,\partial_\xi U_n)
\]

and satisfies

\[
 \partial_\xi W_n+\xi^{-1}D_0W_n
 =A_{0,n}W_n+A_1\partial_\eta W_n+f_n,\qquad W_n(0)=0. \tag{1}
\]

At stage `n`, every lower coefficient is already globally finalized. Hence
`A_(0,n)`, `A1` and `f_n` are fixed before the new Borel scale `a_n` is
chosen. The regular formulas for `M_k/X`, `V_k/X` and the pressure source
remove every apparent axis division. In particular, the order-two pressure
source retains the shifted axial viscosity of `V0`.

The derivative block and all its `eta` derivatives map the first four
coordinates into the last two and annihilate the last two. Diagonal Volterra
kernels preserve those subspaces. Consequently two adjacent derivative
letters vanish, and a word of length `k` has at most

\[
 p_k=\lceil k/2\rceil \tag{2}
\]

surviving Cauchy losses. On a smaller parameter strip of loss `Delta`, with
finite coefficient/source bound `C_n` and `z=2C_na`, the ordered radial
simplex gives

\[
 {z^{k+1}\over(k+1)!}
 \max\{1,p_k/\Delta\}^{p_k}. \tag{3}
\]

For `k=2j`, use `j^j<=e^j j!` and
`(2j+1)!>=(j!)^2`; for `k=2j+1`, use the corresponding `(j+1)` bounds.
Thus (3) is summable for every finite `C_n`. No bound uniform in `n` and no
common complex strip for all orders is required. The same majorant proves
uniqueness, parity at the axis, and finiteness of every fixed coefficient
jet on the common radial collar.

## 2. Global repair does not introduce a circular scale choice

After the collar solve, only `F_n,U_n` are cut off. Forward reconstruction
restores `M_n,V_n,Pi_n`; pressure is not truncated on the reserved mean patch.
The current five moment rows are linear because current-order quadratic terms
occur at order `2n`. On the fixed ordered bumps, the two-row `U` determinant
is strictly positive by monotonicity of its two weights. The three exponential
`E` weights form a Chebyshev system, so their determinant has constant nonzero
sign on ordered triples. Integration against positive disjoint unit-mass bumps
therefore gives an invertible matrix. Its inverse and the finite incoming
discrepancy produce finite smooth bump coefficients.

This establishes, before selecting `a_n`, finite seminorms `M_(n,m)` for every
fixed `m`. For a fixed smooth cutoff `chi`, Euler differentiation of
`chi(a_nq)` becomes `(sigma partial_sigma)^j chi(sigma)`, `sigma=a_nq`;
its norm is independent of `a_n`. Including the radial curl commutator gives

\[
 C_{n,m}=\max\{1,2K_{m+1}(1+2nh)^mM_{n,m}\}<\infty. \tag{4}
\]

The extra cutoff derivative in `K_(m+1)` is essential for the radial entry.
Complete forward pressure, every `Pi_n`, and the commutator remain in (4).

## 3. The real Borel selector closes after (4)

Put `b_n=log a_n`, `ell=log 2`, and `v_(n,m)=log C_(n,m)`. After (4) is
known, choose the real number

\[
 b_n\ge1+\max\left\{b_{n-1}+\ell,\ B+\ell,
 {v_{n,m}+n\ell\over nh}:0\le m\le n,
 \text{the finite physical rows}\right\}. \tag{5}
\]

Every entry in this maximum is already finite. Therefore (5) is an inductive
definition, not a circular appeal to the completed infinite field. It gives
`a_n>2a_(n-1)` and, throughout `q<=a_n^-1`,

\[
 C_{n,m}q^{nh}\le2^{-n}e^{-nh}<2^{-n},\qquad m\le n. \tag{6}
\]

The coefficient increment has the additional `q^(nh)` factor, so (6) gives
the required diagonal gain. Geometric separation makes the cutoff sum locally
finite at every fixed `q>0`; the remaining tail is geometric in every fixed
jet. No infinite numerical scale list or stage-uniform constant is needed.

This reconstructs the coefficient-inverse and Borel-selector part of the
parent's first contract while retaining complete pressure memory.

## 4. Precise remaining leading-profile leaf

The derivation still needs one finalized leading profile, before any positive
coefficient is constructed. For every fixed radial order `r`, it requires a
strip width `delta_r>0` such that

\[
 \sup_{0\le X\le a^2,\ \operatorname{dist}(\eta,[-1,1])<\delta_r}
 \left|\partial_X^{\le r}\partial_\eta^{\le r}
 (F_0,U_0,\Pi_0,M_0/X,V_0/X)\right|<\infty. \tag{7}
\]

It also requires, for the same completed profile, the exact `I3` power-law
patch and five ordered bumps used by the moment inverse, plus compact velocity
and pressure-source identities after the exterior join. The strip may shrink
with `r`; a common all-order analytic width is unnecessary.

Existing core, activation, continuation, modulation, repair and exterior
checkpoints supply a route to (7), but this checkpoint does not rederive (7)
from those defining maps. Equation (7) and the stated patch identities are now
the precise first obstruction, rather than a vague requirement for all
background seminorms or an infinite scale sequence. Taking `q` smaller cannot
repair a failure of (7).

No unconditional pulse, complete residual, admissible force, competitor
exclusion or blow-up theorem is claimed. The next obligation is to prove (7)
and the patch/support identities for the same finalized leading profile, then
advance the same-band gap audit without sampled continuum infima.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_leading_background_seminorm_reconstruction/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_leading_background_seminorm_reconstruction/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_leading_background_seminorm_reconstruction/evidence.json \
  --output fresh.json
```
