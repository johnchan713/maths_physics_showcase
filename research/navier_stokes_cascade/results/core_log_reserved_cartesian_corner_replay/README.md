# Log-reserved Cartesian-corner residual replay

Parent: `2a022f7951a72abee335bc1cac3f524fb3e01e7b`.
Status: `log-reserved-six-group-residual-transferred-to-cartesian-corner-conditionally`.

The parent closes the complete six-group local residual for the existing
locally finite field by keeping separate logarithmic reserves for the two
compact maps, the regular residual, and the cutoff tail.  This checkpoint
replays that corrected estimate on the closed terminal chart and transfers
it to Cartesian corner flatness.

The replay corrects one inherited statement without changing historical
evidence.  The scale-compatible corner checkpoint wrote a no-log,
full-exponent fixed-stage residual supremum.  The later actual-source audit
shows that finite logarithmic factors need not be bounded at zero reserve.
The correct closed-chart input is log-weighted, and the parent's extra
recurrence stage pays that log.  The cutoff tail continues to use its own
reserved half-gain.

These are analytical implications of the frozen primitive and operator
contracts.  The exact-arithmetic tests check exponents, chart identities and
dependency provenance; they do not prove the continuum estimates by
computation.

## 1. Correct closed-chart input

On the terminal chart,

\[
 T+\eta^2=1,\qquad |\eta|\le1,
 \qquad L=1-2h\eta^2\ge1-2h=\frac{99}{100}.             \tag{1}
\]

For every fixed stage `J`, jet `m`, compact normalized radial interval `K`,
and finite source graph, the corrected actual-source domain gives a finite
degree `P_(J,m,K)` and

\[
 \sup_{K,|\eta|\le1}
 q^{-(hB_J-K_m)}(1+|\log q|)^{-P_{J,m,K}}
 |\mathcal R_{\rm regular}(U^{[J]})|_m<\infty .         \tag{2}
\]

The supremum includes `eta=+/-1`.  It uses the uniform background gaps,
normalized support separation `3/8012`, strict dimensionless root
contractions, fixed normalized support hulls, and upper weighted edge
seminorms.  It does not use a positive physical support-distance floor or a
positive pointwise edge-weight floor.

The earlier unweighted expression

\[
 \sup q^{-(hB_J-K_m)}|\mathcal R(U^{[J]})|_m<\infty     \tag{3}
\]

at the same exponent is not inferred from (2).  A factor
`1+|log q|` already disproves that inference.  Equation (3) is therefore not
used in this replay.

## 2. Spend the regular reserve only once

Let

\[
 \delta=\frac12-3\kappa>0,
 \qquad B_{J+1}-B_J=\delta .                             \tag{4}
\]

For a requested Cartesian power `N`, choose the inherited comparison stage
`J0` and use `J=J0+1`.  The exact extra regular gain is `h delta`.  Spending
half in the logarithmic Young bound gives

\[
 (1+|\log q|)^{P_{J,m,K}}q^{hB_J-K_m}
 \le C_{J,m,K}q^{hB_{J_0}-K_m+h\delta/2}.               \tag{5}
\]

Each of `A0` and `A1` separately spends half of one extra compact
integration on its own source log.  The cutoff sequence separately spends
half of every representative gain `g_j=h delta j`.  No reserve is counted
twice.

The exact diagonal difference retains

\[
 e_t-\Delta e+\nabla\pi
 +(U^{[J]}\!\cdot\nabla)e+(e\!\cdot\nabla)U^{[J]}
 +(e\!\cdot\nabla)e .                                   \tag{6}
\]

Thus time, viscosity, pressure, old-new, new-old and new-square all remain.
The fixed gauge retains the moving-normal and `K` pressure, all six covariance
entries, `g_d-A_1g_d`, the full `partial_Z(J_z+c_rho P)`, and both exceptional
compact pressure tails.

For every fixed `m,N` the parent's comparison gives, uniformly on the closed
chart and fixed `K`,

\[
 |\mathcal R(u_*,p_*)|_m\le C_{m,N,K}q^{N+1}.           \tag{7}
\]

For `h=1/200`, `kappa=1/100000`, and `m=N=0`, the same log-reserved fixed
comparison stage is `J=12817`.  The minimum retained exponent is
`20030773/20000000>1`; it comes from the nonlinear tail crosses, while the
fixed-state regular reserve is much larger at this deliberately coarse
comparison stage.

## 3. Cartesian corner transfer

Put

\[
 \rho=\tau+|z|^{1/D},\qquad D=\frac12-h .               \tag{8}
\]

The unchanged implicit-coordinate calculation gives

\[
 q\le4\rho .                                             \tag{9}
\]

Therefore (7) yields

\[
 |\mathcal R(u_*,p_*)|_m
 \le 4^{N+1}C_{m,N,K}\rho^{N+1}.                        \tag{10}
\]

Because `N` is arbitrary and the comparison stage is fixed only after the
jet and target are chosen, every fixed residual jet is flat at the Cartesian
terminal corner.  There is no substitution `J=J(q)`, no raw-series
convergence assumption, and no stage-uniform coefficient bound.

This conclusion concerns the residual jets.  It does not reprove the
punctured one-sided field and pressure jets, which remain in their historical
checkpoint and must be joined to (10) before localization.

## 4. Scope and next obligation

Equations (2), (5), (7), and (10) establish conditional Cartesian corner
flatness for every fixed jet of the complete residual on the same field,
cutoff sequence, supports and gauges.  No unconditional pulse bound,
localization theorem, compact force, comparison, competitor exclusion, or
blow-up result is asserted.

The next obligation is to join the inherited punctured one-sided potential,
direct-swirl and fixed-gauge pressure jets to this corrected residual
flatness.  Then localization must be replayed with the extra time derivative,
unit viscosity, every cutoff derivative, both nonlinear crosses and the
correction square before force, energy and comparison can be advanced.

Reproduce from the repository root:

```bash
python3 research/navier_stokes_cascade/results/core_log_reserved_cartesian_corner_replay/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_log_reserved_cartesian_corner_replay/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_log_reserved_cartesian_corner_replay/evidence.json \
  --output core_log_reserved_cartesian_corner_replay_fresh.json
```
