# Closed-corner cutoff-tail slack reconciliation

Parent: `3dcfe37baa7005b4d513be03a56988b07bbb2ff3`.
Status: `closed-corner-cutoff-tail-slack-reconciled`.

The parent correctly proves the terminal-chart geometry and identifies the
need for closed-chart, compact-`X` estimates.  It overstates one part of the
remaining obligation: its equation (9) asks for a no-logarithm cutoff-tail
constant at the full representative exponent, uniformly in the cutoff
parameter.  The inherited representative estimate contains a fixed power of
`1+|log q|`; at the full exponent that normalized supremum need not be finite.
It is also unnecessary.  The complete-potential diagonal reserves exactly
half of every positive stage gain to absorb finite constants and logarithms.

This checkpoint corrects that quantifier and proves that the cutoff operation
adds no further closed-corner obstruction once the actual representative
coefficient estimate is available.  Historical evidence remains frozen.

## 1. The full-exponent no-log supremum is too strong

For stage `j>=1`, retain

\[
 g_j=h\delta j>0,\qquad \delta={1\over2}-3\kappa,
\]

and the inherited complete representative estimate

\[
 |(V_{j,a},P_{j,a})|_m
 \le C_{j,m}(1+|\log q|)^{P_{j,m}}q^{g_j-L_m}.       \tag{1}
\]

The constant is independent of `a` on `aq in [1/2,1]`, conditional on the
closed-chart representative seminorms.  At the exact exponent used in the
parent, (1) only gives

\[
 q^{-(g_j-L_m)}|(V_{j,a},P_{j,a})|_m
 \le C_{j,m}(1+|\log q|)^{P_{j,m}}.                  \tag{2}
\]

When `P_{j,m}>0`, the right side is unbounded as `q` tends to zero.  The model
`q^(g_j-L_m)(1+|log q|)^P` satisfies (1) with equality and refutes any inference
of a no-log constant from that estimate.  It is an analytical majorant model,
not a candidate Navier--Stokes field.

Thus the parent's full-exponent cutoff-tail supremum is neither established
nor required.  Keeping the logarithm is essential.

## 2. Half of the gain uniformly absorbs every fixed logarithm

Put `s=-log q>=1` and `gamma_j=g_j/2`.  After normalizing at the exponent
`g_j/2-L_m`, equation (1) gives

\[
 q^{-(g_j/2-L_m)}|(V_{j,a},P_{j,a})|_m
 \le C_{j,m}(1+s)^{P_{j,m}}e^{-\gamma_js}.            \tag{3}
\]

The elementary inequalities

\[
 \log(1+s)\le\sqrt{s},\qquad
 P\sqrt{s}\le{\gamma_js\over2}+{P^2\over2\gamma_j}
\]

imply the explicit uniform bound

\[
 \sup_{s\ge1}(1+s)^Pe^{-\gamma_js}
 \le \exp\!\left({P^2\over2\gamma_j}\right)
 =\exp\!\left({P^2\over g_j}\right).                 \tag{4}
\]

This proves the corrected closed-tail implication

\[
 \sup_{\substack{X\in K,\ |\eta|\le1,\ T=1-\eta^2\\
                   aq\in[1/2,1],\ a\ge e}}
 q^{-(g_j/2-L_m)}|(V_{j,a},P_{j,a})|_m<\infty,        \tag{5}
\]

provided the coefficient `C_{j,m}` and log power are finite uniformly on the
closed chart.  The implicit-coordinate denominator obeys `L>=49/50`, and on
the transition each `a` is paired with `q` through `aq`; weighted cutoff
derivatives therefore retain constants independent of `a`.  No division by
`1-eta^2` or `T` is introduced.

## 3. The existing scale recursion is already the stronger form needed

The actual construction selects `a_j` after the finite stage data are fixed.
For `c>=log C`, it chooses `s_j=log a_j` above

\[
 \max\left\{s_{j-1}+1,\log(q_*^{-1})+1,1,
 {2(c+j)\over\gamma_j}+{P^2\over\gamma_j^2}\right\}. \tag{6}
\]

The same Young bound then gives on the whole interval `0<q<=1/a_j`

\[
 C(1+|\log q|)^Pq^{\gamma_j}\le2^{-j},               \tag{7}
\]

and hence

\[
 |(V_{j,a_j},P_{j,a_j})|_m
 \le2^{-j}q^{g_j/2-L_m}.                               \tag{8}
\]

Equation (8), including every cutoff transition, is precisely the tail input
used by the fixed-stage comparison.  It does not require uniformity over
arbitrary unselected cutoffs, a no-log estimate at the full exponent, or a
stage-uniform constant.

## 4. The corrected remaining obstruction

The cutoff-tail transfer is now closed relative to two primitive
closed-chart coefficient statements:

\[
 C^{\rm stage}_{J,m}:=
 \sup_{\substack{0<q\le q_0,\ X\in K,\ |\eta|\le1\\T=1-\eta^2}}
 q^{-(hB_J-K_m)}|\mathcal R(U^{[J]})|_m<\infty,        \tag{9}
\]

and, before spending the half gain,

\[
 C^{\rm rep}_{j,m}:=
 \sup_{\substack{0<q\le q_0,\ X\in K,\ |\eta|\le1\\T=1-\eta^2}}
 {q^{-(g_j-L_m)}|Z_j|_m\over(1+|\log q|)^{P_{j,m}}}<\infty. \tag{10}
\]

Here `Z_j` is the actual fixed-gauge complete increment representative before
the cutoff derivative is applied.  Existing finite-stage notes state (9)--(10)
conditionally through retained background, root, inverse, support and
physical-transfer contracts.  Frozen arithmetic and graph evidence does not
independently prove those continuum coefficient bounds at the angular poles
and uniformly on compact `X` intervals.

Once (9)--(10) are proved for the actual providers, equations (4)--(8) close
the cutoff tail, the parent's `q`-to-Cartesian transfer closes the corner, and
the existing localization, force and comparison implications can be replayed.
No unconditional compact force, complete residual theorem, admissible-force
exclusion, or blow-up result is claimed here.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_corner_tail_slack_reconciliation/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_corner_tail_slack_reconciliation/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_corner_tail_slack_reconciliation/evidence.json \
  --output core_corner_tail_slack_reconciliation_replay.json
```
