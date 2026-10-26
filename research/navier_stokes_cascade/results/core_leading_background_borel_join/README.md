# Same leading profile and rule-defined background Borel join

Parent: `f3f2068e84261a7a0cd3ce4bf488871d56a35e88`.
Status: `same-leading-profile-and-rule-defined-background-borel-contract-joined`.

The parent isolated five primitive contracts beneath the corrected diagonal.
This checkpoint audits the first one: the completed leading profile, the
all-order background coefficient recursion, and the exact-real Borel
prescription refer to one construction. It joins their continuum estimates
without selecting a new scale sequence or replacing pressure by a leading
centrifugal term.

A scope correction is necessary. The leading profile is holomorphic in the
angular variable on a protected collar through the coefficient-construction
region. Globally, including the terminal heat tail, it has every fixed smooth
one-sided endpoint jet; it is not asserted complex analytic there. The
background recursion needs the local collar and global fixed-order smoothness,
so no global heat analyticity premise is required.

## 1. The same completed leading input reaches the coefficient recursion

The independently rechecked chain constructs one profile in this order:
the coefficient-space core fixed point, positive source and native exit,
flat activation collar, prescribed continuation, same matching roots,
finite-frequency modulation with its moving phase, exact five-row repair,
and the heat exterior. The global moment recheck retains

\[
 M_\infty=R=J_\infty=S_\infty=0,
 \qquad \Pi_{axis}+C_{p,\infty}=0.                         \tag{1}
\]

Thus the leading input to `background_jet_foundation` has its required
regular axis extensions, protected angular collar, reserved I3 coefficient
patch, disjoint I4 mean patch, compact positive-order support region and
same-gauge backward pressure. No leading root, support, width or pressure
datum is reselected here.

The word “analytic” in the primitive frontier is therefore read precisely:
holomorphic angular extension on the protected finite collar, and all fixed
smooth jets on every remaining closed piece. The heat endpoint may have zero
Taylor radius. This correction preserves every estimate actually used below.

## 2. Every fixed background coefficient exists before Borel selection

At positive order `n`, with `lambda_n=2nh`, the zero-axis unknown

\[
 W_n=(F_n,U_n,M_n/X-U_n,\Pi_n,\partial_\xi F_n,
             \partial_\xi U_n),\qquad \xi=\sqrt X,
\]

satisfies

\[
 \partial_\xi W_n+\xi^{-1}{\rm diag}(0,0,2,0,3,1)W_n
 =A_{0,n}W_n+A_1\partial_\eta W_n+f_n.                    \tag{2}
\]

The source contains all lower products, radial viscosity and the shifted
axial-viscosity contribution. Apparent axis quotients are interpreted by
their integral extensions. The derivative matrix maps the first four
coordinates to the last two and annihilates the last two, including after
angular differentiation. Hence adjacent derivative factors in the Volterra
series vanish. A length-`k` word has at most `ceil(k/2)` Cauchy losses, while
the ordered radial simplex contributes `1/(k+1)!`. The resulting series is
summable for every finite coefficient bound. This constructs every fixed
`n` on the same radial collar; neither a uniform complex width nor a bound
uniform in `n` is needed.

After the radial cutoff, the same five I3 rows are linear in the current
coefficient. Their normalized matrix is independent of `n` and `eta`; the
ordered U weights have positive determinant, while the three distinct E
rates form a Chebyshev system. The established inverse bounds therefore
apply at every order. Forward reconstruction restores `M_n,V_n,Pi_n` and
all five moments, retaining the full pressure source. Positive-order
velocities have one common compact support; pressure is never cut off on the
zero-velocity I4 patch.

This induction is completed before choosing any Borel scale. Consequently
each finite Euler seminorm `C_(n,m)` of the complete coefficient tuple,
including the radial curl commutator, exists as an exact continuum supremum.

## 3. The exact-real prescription identifies the summed background

Let `b_n=log a_n`, `g_n=nh`, `ell=log 2`, and let `v_(n,m)=log C_(n,m)`.
The retained prescription chooses a real `b_n` above the finite maximum

\[
 1+\max\left\{b_{n-1}+\ell,\ B+\ell,
       {v_{n,m}+n\ell\over g_n}:0\le m\le n,
       \hbox{all finite physical rows}\right\}.           \tag{3}
\]

Every entry is finite before this choice. Induction on `n` therefore defines
the whole sequence and gives `a_n>2a_{n-1}`. For every `m<=n`, throughout
the support `q<=a_n^{-1}`,

\[
 C_{n,m}q^{nh}\le2^{-n}e^{-nh}<2^{-n}.                    \tag{4}
\]

At fixed `q>0`, scale separation makes the cutoff sum locally finite. For
any fixed jet, the remaining all-smaller-`q` tail is dominated by

\[
 \sum_{n\ge N}2^{-n}t^n={ (t/2)^N\over1-t/2}
 \le2^{1-N},\qquad 0\le t=q^h\le1.                       \tag{5}
\]

The actual rule-defined background is therefore the smooth sum of the same
streamfunction, azimuthal velocity, pressure and stress coefficients.
Cutting the streamfunction produces the essential radial curl-commutator
term; it is included in `C_(n,m)`. The pressure sum retains every `Pi_n`,
including `Pi1/Pi2` memory. Equations (3)--(5) concern exact real suprema and
do not identify the construction with a diagnostic finite `mpmath` list or
an unrelated historical numerical sequence.

This closes the first primitive contract for the rule-defined background,
relative to the already independently rechecked leading construction. It
does not verify the next contract: phase, normal, value and cone gaps must
still be joined on this same accepted tail and with the actual leading first
jets. The canonical full-path ODE, fixed compensator/projector supports,
compact/slow inverse domains and complete viscosity/cutoff/interaction
residual remain later obligations. No unconditional pulse bound, force,
residual improvement or blow-up result follows.

Reproduce with

```sh
python3 research/navier_stokes_cascade/results/core_leading_background_borel_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_leading_background_borel_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_leading_background_borel_join/evidence.json \
  --output /tmp/core_leading_background_borel_join_fresh.json
```
