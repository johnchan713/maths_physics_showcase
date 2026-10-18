# Complete harmonic sources: a count-independent bound and a necessary loss

Parent: `0c1399dca354f460b80c2b22ef80f63ea285fc95`.
Status: `complete-finite-harmonic-closure-with-explicit-mode-derivative-budget`.
The parent hosted run was still in progress, with no failed step reported,
when inspected on 2026-10-17. Historical evidence is preserved.

The current branch already has finite covariance/mean cycles and a conditional
diagonal realization. Repeating those cycles would not discharge their
dependencies. This checkpoint instead quantifies a previously implicit
dependence on harmonic number in the generated COMPLETE sources. It gives
weighted absolute-sum estimates whose constants do not depend on the number
of modes or their largest index, conditional on the reviewed coefficient,
frame, entrance-support and inverse bounds. It also proves an obstruction:
complete incompressibility removes the large carrier k from same-phase
transport, but does **not** remove the ratio of harmonic indices. One angular
coefficient derivative is necessary in a generic uniform transport bound.

These are analytical estimates with exact manufactured algebra controls.
They neither evaluate the actual giant profile nor prove summability of an
infinite nonlinear sequence. No frequency, original root, tilt, width C^-20,
support, pressure gauge, viscosity or operation order is changed. The retained
exact-core/profile and remaining complete construction dependencies persist.

## 1. A norm that keeps the harmonic cost visible

In one fixed label write the complete real wave as
`u=sum_(m!=0) a_m exp(i k m Phi)`, with `a_-m=conj(a_m)`.
Different labels have separated CLOSED enlarged supports. Use the current
common lifted path and the same accepted bands, `epsilon=Q^h`, `S=ell^2`,
`1/4<=q/Q<=4`, `4Q<qstar`. Discrete labels are frozen before differentiation.

For s>=0 and a fixed coefficient jet order j define

\[
 N_{s,j}(a)=\sum_m(1+|m|)^s\|a_m\|_j.                 \tag{1}
\]

Here `||.||_j` is the sum of suprema of ordinary coefficient derivatives
through order j on the common path/chart. Tensor/vector norms use component
absolute sums. For the inverse use the DISTINCT envelope seminorm
`||a||_(j,W)=sum_(|I|<=j) sup |D^I a|/[sqrt(zeta) P d_edge^(-b_j)]`
on the open supported shell and WHOLE lifted path, with the retained
finite inverse-edge exponents b_j (separate distances at the two edges).
Its weighted sum is denoted `N^W_(s,j)`. Additional polynomial S factors
and increased finite edge/jet losses appear in the operator constant and
output seminorm. This is not an unweighted supremum estimate for forcing.
For a finite collection choose common maximum degrees. For infinitely
many modes, such majorants UNIFORM in m and convergence of the weighted
sums are extra hypotheses, not consequences of smoothness at each mode.

For coefficient products, Leibniz gives a j-dependent constant. Since
`1+|m+n| <= (1+|m|)(1+|n|)`, convolution satisfies

\[
 N_{s,j}(ab)\le C_j N_{s,j}(a)N_{s,j}(b),             \tag{2}
\]

including the zero output. This follows by the triangle inequality and
interchanging finite nonnegative sums; there is no factor equal to a mode
count. It also holds for absolutely convergent inputs with the stated jets.

## 2. The full source inverse gains two indices; pressure gains one

The [leading/frame floors](../leading_frame_foundation/README.md) now supply
the nonzero normal and specified left-inverse floors used in the
[explicit sourced inverse](../sourced_inverse_foundation/README.md).
Its high-mode estimate `1/(m^2-1)` is at most `4/(3m^2)` for |m|>=2.
The two fundamental modes have the original polynomial-length bound.
Absorbing this finite low-mode bound into C_j S^b, that inverse gives

\[
 N^W_{s+2,j}(t)\le C_j S^{b_j} N^W_{s,J(j)}(f),\qquad
 k N^W_{s+3,j}(\mathcal C)\le C_j S^{b_j}N^W_{s,J(j)}(f),
\quad k N^W_{s+1,j}(\pi)\le C_j S^{b_j}N^W_{s,J(j)}(f).
                                                               \tag{3}
\]

Fixed background/edge constants are included; J(j) is a finite sufficient
input jet count, independent of m. Ordinary coefficient, pulse and mixed
jets require the full smooth entrance compatibility already supplied by
the regular supported-source grammar. Without it even y'(0)=1 for a
constant normalized forcing blocks the claimed pulse-jet smoothing.
The pressure formula is the FULL

\[
 \pi_m={i(n\cdot Kt_m-n'\cdot t_m+n\cdot f_m)\over km|n|^2},
 \qquad \mathcal C_m={i(n\times t_m)\over km|n|^2}.    \tag{4}
\]

For a pure normal smooth supported source, t=0 and pi=i f_normal/(km|n|).
Thus pressure cannot inherit the extra two indices of the potential.
This remains true with finite smooth endpoint support.

Cut the SAME potential and pressure before taking the full curl. Its
amplitude is `a_m=psi*t_m+curl_*(psi*C_m)`. The second summand gains
three indices, carries k^-1, and costs one displayed graph derivative.
Consequently the COMPLETE cut velocity retains the two-index bound (3).
An uncut exit is not declared zero. Ordinary physical derivatives of the
carrier still cost k|m|; (3) is a coefficient estimate.

One can also make the m dependence of the reviewed REGULAR linear
remainder explicit. Its largest harmonic powers are: `m^2` on the curl
remainder, `km` on the principal velocity in eikonal/mixed-viscous terms,
and no m on the amplitude gradient of pressure. Equations (3)--(4) give
one-index smoothing in every case. With the current epsilon weights,

\[
 N^W_{s+1,j}(L_{reg}\delta u)\le
 C_j S^{b_j}\epsilon^{1/2-\kappa}N^W_{s,J(j)}(f).     \tag{5}
\]

This is the same sufficient epsilon gain as the regular-source foundation,
now with explicit harmonic control. Terms with two graph derivatives
have epsilon^(1-2kappa), which is stronger since kappa<1/2. Both viscosity
and pressure are retained. Exact cutoff commutators and (1-psi)f stay in
the separate additive Gaussian ledger, never in the next inverse input.
Actual cutoff velocities remain in every nonlinear product.

The full-path envelope hypothesis is necessary: even the scalar model
`t'=t+f`, `t(0)=0`, `P(v)=exp(v-L)<=1` has `t(L)=exp(L)-1` for f=1.
Thus an unweighted forcing supremum of one cannot give a polynomial-length
bound. An early smooth supported positive f has the same exponential
amplification. Its envelope-normalized supremum is correspondingly large.
The accepted regular sources supply the stronger hypothesis; arbitrary
physically flat additive tails are not substituted for them.

## 3. Incompressibility removes k, with a necessary remaining index loss

For a COMPLETE same-label amplitude,
`div_* a_m+i k m n.a_m=0`. Hence its transport of b_n is exactly

\[
 a_{m,r}D_rb_n+a_{m,z}D_zb_n
 +R^{-1}(-a_{m,theta}b_{n,theta},a_{m,theta}b_{n,r},0)
 -{n\over m}(div_*a_m)b_n.                            \tag{6}
\]

No output index is divided out, so n=-m and zero output cause no singularity.
Use |n/m|<=1+|n|, (2), and the bounded radial interval. A sufficient bound is

\[
 N_{s,j}(B(a,b))\le C_j\epsilon^{-\kappa}
 N_{s,j+1}(a)N_{s+1,j+1}(b).                          \tag{7}
\]

The epsilon^-kappa charges the displayed graph derivative once. The
symmetric interaction has (7) plus the exchanged term. It keeps P^2<=P,
the full edge weight, smooth endpoints and the closed label support.
Different phases cannot use (6); their derivative products vanish by the
actual separated supports. Mean-wave interaction retains its separate
half-frequency loss and is not bounded by (7).

The extra harmonic derivative in (7) is necessary for a generic bound
uniform in the largest mode. In Cartesian coordinates with common phase z
and k=1 take complete curls

\[
 a=(f(x),0,if'(x))e^{iz}
   =curl\big(if(x)e_y e^{iz}\big),\qquad
 b=e_y e^{iNz}.
\]

Both divergences vanish, but `(a.grad)b=-N f'(x)e_y e^{i(N+1)z}` and
`(b.grad)a=0`. At a point with f'=1, input coefficient jets are bounded
independently of N and transport grows as N. Choose f smooth compactly
supported and equal to x on a small interval. A compact vector-potential
cutoff for b can equal one on the same local box; it preserves the local
identity. Adding conjugates gives real complete curls, with nonzero modes
N+1 and N-1 of size proportional to N. This is a local obstruction to a
generic coefficient bound, not a selected actual-pulse counterexample.
The old finite-stage statement that m ratios enter finite constants remains
valid. It must not be read as a harmonic-uniform no-loss estimate.

## 4. Complete covariance, pressure memory and what this closes

Angular averaging keeps every matched pair:
`W_ij=sum_m a_(m,i) a_(-m,j)`. In particular all six symmetric entries,
the high-harmonic squares, axial covariance and mixed terms remain.
Equation (2) implies, without a mode-count cost,

\[
 \|C(U+V)-C(U)\|_j\le
 C_j\,[2N_{0,j}(U)+N_{0,j}(V)]N_{0,j}(V).             \tag{8}
\]

The same bound holds after auxiliary Haar averaging. A signed difference
need not be positive; the new square is retained. Changing waves while
holding current means fixed still requires CURRENT pressure reconstruction:

`Delta gr=-D1 Delta Wrr-Dz Delta Wrz+Delta Wtheta theta/R`,
`Delta pm=T0(Delta gr-rho Delta P)`, `Delta P=int bar(Delta gr)`.
The exact pressure moment and axial flux contain
`int R bar(Delta pm)=-1/2 int R^2 bar(Delta gr)+c_rho Delta P` and
the WHOLE derivative `partial_Z(Jz+c_rho P)`. The two velocity moments,
old Haar means and actual A0/A1 velocities remain. The same full background
pressure with Pi1 and Pi2 cancels only in subtraction of that same base
equation. It is not replaced by leading pressure.

Thus the existing finite source grammar and full six-entry covariance
algebra now have explicit count-independent majorants (2)--(8), with
the necessary angular loss exposed. Every finite stage has all the needed
higher finite moments by its finite support; no new first band follows.
For an infinite mode family or a uniform stage estimate one still needs
the right-hand weighted sums, common endpoint/edge majorants and constants.
The sourced inverse's two-index smoothing alone does not prove this for
the nonlinear iteration; the old/new terms can spend an index and the
mean-wave frequency factor remains. The existing diagonal realization
uses finite stages and does not require a new infinite Fourier inverse.

Next independently join the current leading/root/inverse, generated targets
and fresh mean maps into a dependency-explicit complete residual estimate,
retaining the exact core and common all-scale domain. The inherited core
and complete construction obligations must be discharged before any
unconditional residual, admissible-force exclusion or blow-up claim.

## 5. Reproduction and limits of the controls

`mode_bounds.py` checks rational convolution/transport majorants, inverse
index factors and the linear term ledger. `mode_checks.py` independently
forms exact rational complex complete curls, checks divergence and direct
transport against (6), and verifies every covariance entry against Fourier
convolution, including negative modes and high-mode squares. The real
N-mode obstruction is exact for arbitrary integer N>=2. Fixtures are local
polynomials, not actual smooth global profiles or continuum certificates.
The audit binds the frozen current dependency sources without promoting
their conditional flags; separate fresh dependency replays are reported
with publication validation. No historical evidence is overwritten.

```sh
python3 research/navier_stokes_cascade/results/harmonic_closure_foundation/test_audit.py -v
python3 research/navier_stokes_cascade/results/harmonic_closure_foundation/audit.py \
  --verify-record research/navier_stokes_cascade/results/harmonic_closure_foundation/evidence.json \
  --output /tmp/harmonic_closure_fresh.json
```
