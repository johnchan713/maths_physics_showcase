# Closed-corner transfer and the exact uniform-constant obstruction

Parent: `369f9a43e1cf0b2762f3c4144337d17b855ebb87`.
Status: `closed-corner-transfer-proved-and-uniform-constant-obstruction-isolated`.

The parent identifies the exact core pressure with the completed leading
backward pressure on the same protected field.  The current downstream
obligation is the residual estimate needed to turn conditional localization
into one smooth compact force.  This checkpoint separates the harmless
coordinate transfer from the still-unproved uniform continuum input and
states the latter as two precise closed-chart inequalities.

The same complete-potential field, gauges, supports and cutoff sequence are
retained.  No historical record is changed, and no pathwise or finite-grid
diagnostic is promoted to a compact supremum.

## 1. Exact terminal chart geometry

Keep

\[
 q=\tau+z^2q^{2h},\qquad X={r^2\over2q},\qquad
 \eta={z\over q^D},\qquad T={\tau\over q},\qquad
 D={1\over2}-h.                                         \tag{1}
\]

Since `2D+2h=1`, division of (1) by `q` gives the exact closed-chart
identity

\[
                         T+\eta^2=1.                     \tag{2}
\]

Thus the terminal boundary `T=0` is precisely `|eta|=1`; they are not two
independent missing faces.  Implicit differentiation uses

\[
 L=1-2h\eta^2\ge1-2h\ge{49\over50}>0,                  \tag{3}
\]

so the coordinate map itself has no pole on the closed angular interval.

Put `y=|z|^(1/D)`.  If `q<=2tau`, then `q<=2tau`.  Otherwise
`z^2q^(2h)>q/2`, so

\[
 q^{2D}<2z^2=2y^{2D},\qquad
 q<2^{1/(2D)}y<4y,                                      \tag{4}
\]

because `2D>=49/50>1/2`.  In both cases

\[
 q\le2\tau+4|z|^{1/D}\le4(\tau+|z|^{1/D}).             \tag{5}
\]

Consequently a **uniform** physical residual estimate

\[
 \sup_{X\in K,\ |\eta|\le1,\ T=1-\eta^2}
 |D_x^\alpha\partial_t^j\mathcal R|\le C_Pq^P          \tag{6}
\]

immediately implies the localization estimate with constant `4^P C_P`.
The extra time derivative needed for Cauchy convergence is another fixed
finite request, so arbitrary `q`-flatness absorbs it by choosing `P` after
the derivative order.

## 2. What the locally finite diagonal already supplies

For a fixed jet `m` and target power `P`, the diagonal comparison chooses a
single finite stage `J`, held fixed as `q` tends to zero, so that

\[
 hB_J-K_{m+2}\ge P+1,
 \qquad {h\delta(J+1)\over2}-L_{m+2}\ge P+1.            \tag{7}
\]

Both slopes in `J` are positive.  The exact controls recompute finite stages
for representative `(m,P)` pairs; the conservative stage can be enormous,
but existence is enough.  No stage-growth estimate and no substitution
`J=J(q)` is required.

If the constants in the fixed-stage and cutoff-tail estimates are uniform
on the closed chart, equations (7), the complete residual difference, and
(5) prove every localization corner jet.  Pressure, viscosity, both
nonlinear crosses, the new square and the cutoff-curl terms remain included.

## 3. The two missing inequalities

The scale-selection argument uses finite constants before selecting each
cutoff scale.  To include the terminal poles, it needs the actual quantities

\[
 C^{\rm cl}_{J,m}:=
 \sup_{\substack{0<q\le q_0,\ X\in K,\ |\eta|\le1\\T=1-\eta^2}}
 q^{-(hB_J-K_m)}|\mathcal R(U^{[J]})|_m<\infty,          \tag{8}
\]

and, uniformly in the cutoff parameter `a` on its transition,

\[
 \widetilde C^{\rm cl}_{j,m}:=
 \sup_{\substack{0<q\le q_0,\ X\in K,\ |\eta|\le1\\
                  T=1-\eta^2,\ aq\in[1/2,1]}}
 q^{-(h\delta j-L_m)}|(V_{j,a},P_{j,a})|_m<\infty.      \tag{9}
\]

Here `K` is every fixed compact `X` interval required by localization.  The
existing finite-stage and diagonal notes state the corresponding interior or
conditional estimates.  Their frozen evidence does not independently
establish (8)--(9) through `|eta|=1`, nor compact-`X` uniformity for the
actual infinite provider chain.  The logarithmic scale rule cannot absorb an
infinite or unproved `log C_cl`.

This is a quantifier issue, not a coordinate singularity.  For example

\[
 g_P(q,\eta)={q^{P+1}\over1-\eta^2+q}                   \tag{10}
\]

has order `q^(P+1)` at every fixed interior angle but only `q^P` at the
pole.  Interior constants therefore do not imply a closed-chart constant.
Likewise

\[
 f(X,q)={X\over q}\exp[-(X/q)^2]                         \tag{11}

is flat as `q` tends to zero for every fixed `X>=0`, yet
`f(q,q)=e^-1`; it is not uniformly flat on any compact interval containing
zero.  These are analytical counterexamples to the inference, not candidate
fluid fields.

## 4. Frontier

The next proof must derive (8)--(9) from the actual fixed-profile background,
root, inverse, support and physical-transfer formulas on the full closed
chart.  Once they are finite, the already selected scale recursion absorbs
them, (7) supplies arbitrary powers, and (5) transfers the result to the
Cartesian corner with the extra time derivative.  The force, endpoint jet
gluing, energy and expanding-ball comparison can then be replayed on the same
field without changing pressure or supports.

Until those two uniform estimates are proved, compact smooth forcing and
competitor exclusion remain conditional.  No unconditional pulse bound,
complete residual improvement, admissible-force theorem or blow-up result is
claimed.  Passing the tests below checks algebra, counterexamples and frozen
provenance; it does not compute either continuum supremum.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_corner_flatness_obstruction/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_corner_flatness_obstruction/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_corner_flatness_obstruction/evidence.json \
  --output core_corner_flatness_obstruction_replay.json
```
