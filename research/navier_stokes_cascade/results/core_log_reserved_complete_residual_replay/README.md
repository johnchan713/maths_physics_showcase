# Log-reserved complete residual replay

Parent: `cce45fc0cc34c48a7a9c22406687dfa0a7b3900d`.
Status: `log-reserved-complete-six-group-local-residual-replayed-conditionally`.

The parent proves the logarithmically renormalized weighted source domain for
the actual seed, one sourced/generic cycle and every fixed finite stage.  It
also records that one extra compact integration absorbs the finite logarithmic
degree of an `A0` or `A1` source.  This checkpoint performs the next residual
replay.

There are three independent logarithmic losses.  The extra `A0` or `A1`
integration can pay only its own compact-source loss.  It cannot also pay a
logarithm in the regular residual or in a cutoff tail.  This checkpoint keeps
three disjoint exponent accounts, retains complete pressure memory, and joins
the fixed-stage estimate to the existing locally finite diagonal.  The result
is an analytical implication of the frozen weighted primitive and operator
contracts.  Exact-arithmetic tests freeze the budgets and provenance; they are
not a continuum proof of those contracts or a numerical construction of the
large comparison state.

## 1. The three reserve accounts

Write

\[
 \delta=\frac12-3\kappa>0,\qquad
 B_j=B_0+j\delta,\qquad \epsilon=q^h .                  \tag{1}
\]

For a fixed physical jet `m`, the inherited finite-state loss is

\[
 K_m=2(\tfrac12+h)+\frac32+2m .                         \tag{2}
\]

Every coefficient and actual source used at a fixed stage has some finite
logarithmic degree.  For `s=|log q|`, any reserve `sigma>0` gives

\[
 (1+s)^Pq^\sigma
 \le \exp\!\left(\frac{P^2}{4\sigma}\right).           \tag{3}
\]

The constant in (3) may depend on the fixed stage, jet and target.  No
stage-uniform log degree or constant is assumed.

The reserves are booked as follows.

1. Each compact map is requested with one additional integration.  For each
   of `A0` and `A1` this separately supplies `q^(h kappa)`.  Half absorbs that
   map's actual-source logarithm and half remains strict:

   \[
    q^{h\kappa}(1+|\log\epsilon|)^{P_{e,J,m}}
       \le C_{e,J,m}q^{h\kappa/2},\qquad e=0,1.          \tag{4}
   \]

2. The regular residual receives its own reserve.  First choose a stage
   `J0` satisfying the old pure-exponent threshold, then use `J=J0+1`.
   The exact gain is `h delta`.  Half pays the fixed finite log degree of the
   regular residual and half remains:

   \[
   (1+|\log q|)^{P^R_{J,m}}q^{hB_J-K_m}
       \le C_{J,m}q^{hB_{J_0}-K_m+h\delta/2}.            \tag{5}
   \]

3. The cutoff construction keeps its existing account.  At representative
   stage `j` it spends half of `g_j=h delta j` on finite coefficient,
   derivative and cutoff logarithms and retains `g_j/2`.  This is the account
   used to choose the scale `a_j` and produce the factor `2^-j`; it is not
   reused in (4) or (5).

Thus a log at zero reserve remains a genuine obstruction, while every actual
fixed-stage logarithm is absorbed by a named positive reserve.

## 2. Fixed-stage complete residual

At stage `J`, the exact physical split remains

\[
 \mathcal R_{\rm phys}
 =\mathcal R_{\rm regular,current}
  +(-A_0Q^\circ g_{r,\rm final},0,+A_1Q^\circ E_{z,\rm post})
  +T_{\rm Gaussian,old}+T_{\psi,\rm signed}+T_{\rm base,ext}. \tag{6}
\]

The regular term in (6) includes the inherited current source, complete
signed linear term, both old/signed wave crosses, signed square, old-mean
interaction, new fast feedback and new slow feedback.  The pressure is the
fresh fixed-gauge pressure with the moving normal and `K` terms.  All six
covariance entries, `g_d-A1 g_d`, the full axial derivative
`partial_Z(J_z+c_rho P)`, and both exceptional compact pressure routes remain
in the actual-source graph.

Equation (5) controls the regular portion.  Equation (4), applied separately
to the two compact tails at the requested physical target
`hB_J-K_m`, controls their finite log degrees with extra positive margin.
The retained Gaussian, signed cutoff-time and flat exterior terms are in the
same finite-stage ledger.  Consequently, for every fixed `m,N`, a finite
stage can be chosen so that

\[
 |\mathcal R(U^{[J]})|_m\le C_{J,m,N}q^{N+1+h\delta/2}. \tag{7}
\]

This is a fixed-stage estimate.  It does not assert convergence of the raw
uncut iteration or uniform bounds in `J`.

## 3. Join to the locally finite field

Keep the already selected cutoff sequence and field

\[
 (u_*,p_*)=(u^{[0]},p^{[0]})+\sum_{j\ge1}(V_{j,a_j},P_{j,a_j}). \tag{8}
\]

For fixed `J` below the active cutoff scale, put
`(e,pi)=(u_*,p_*)-(U^[J],P^[J])`.  The exact difference is

\[
 \mathcal R(u_*,p_*)-\mathcal R(U^{[J]},P^{[J]})
 =e_t-\Delta e+\nabla\pi
 +(U^{[J]}\!\cdot\nabla)e+(e\!\cdot\nabla)U^{[J]}
 +(e\!\cdot\nabla)e .                                  \tag{9}
\]

This retains all six required groups: time, viscosity, pressure, old-new,
new-old and new-square.  The cutoff account gives, with the inherited losses
`L_(m+2)` and `H_(m+2)=L_(m+2)+1`,

\[
 \begin{array}{c|c}
 \text{group}&\text{q-exponent}\cr
 e_t,-\Delta e,\nabla\pi&g_{J+1}/2-L_{m+2}\cr
 (U^{[J]}\cdot\nabla)e,(e\cdot\nabla)U^{[J]}
    &g_{J+1}/2-L_{m+2}-H_{m+2}\cr
 (e\cdot\nabla)e&2(g_{J+1}/2-L_{m+2}).
 \end{array}                                             \tag{10}
\]

Choose `J0` by the inherited comparison inequalities for `m,N`, and set
`J=J0+1`.  Equations (5), (7), and (10) then put every term at exponent at
least `N+1`.  Hence, on the same locally finite field and for every fixed
jet and target,

\[
 \boxed{|\mathcal R(u_*,p_*)|_m=O(q^N).}                \tag{11}
\]

For the retained defaults `h=1/200`, `kappa=1/100000` and `m=N=0`, the old
coarse comparison threshold is `J0=12816`; the log-reserved comparison uses
the fixed stage `J=12817`.  This integer is an exponent-budget output, not an
evaluation of 12817 correction fields.

## 4. Scope and next obligation

Equation (11) closes the complete six-group **local residual replay
conditionally on the frozen primitive, support and operator contracts**.
It preserves the same cutoff sequence, field, supports, gauges, pressure
memory and protected inner/exterior regions.  It does not prove stage-uniform
actual pulse bounds, an unconditional residual theorem, a compact force, a
global comparison, or blow-up.

The next justified obligation is to replay the punctured closed-corner and
Cartesian flatness transfer on this corrected same field.  Localization and
compact forcing must then be replayed with viscosity and every cutoff
derivative, followed by energy, fixed-gauge pressure recovery and the
expanding-ball comparison.

Reproduce from the repository root:

```bash
python3 research/navier_stokes_cascade/results/core_log_reserved_complete_residual_replay/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_log_reserved_complete_residual_replay/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_log_reserved_complete_residual_replay/evidence.json \
  --output core_log_reserved_complete_residual_replay_fresh.json
```
