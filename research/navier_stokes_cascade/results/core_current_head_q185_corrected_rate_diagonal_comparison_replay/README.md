# Q185: corrected-rate diagonal, localization and comparison replay

Parent: Q184, `13683a535b9d8625b1dac84c8ce6a18e9cb2b2a1`.
The live branch, Q184, Q181/Q182, the frozen provider evidence and current CI
were read first. Historical evidence, the physical profile, pressure gauge,
support hole, localization plateau, frequency and `C^-20` width are unchanged.

Q184 corrects the generated-source norm used by Q180. At the exact
logarithmically renormalized `W/M/S` level, the source recurrence keeps the
old gain. A pure-power presentation must spend `eta=kappa` to absorb every
fixed polynomial degree in `S=ell^2`, reducing the gain by `kappa`. This
checkpoint replays the downstream diagonal, residual, localization and
comparison arguments in both versions and corrects the one inherited
numerical margin that changes.

Every conclusion remains conditional on the named continuum providers.
Exact arithmetic, tests and evidence hashes verify the replay and its
provenance; they are not proofs of those providers.

## 1. Two valid routes

Keep

\[
 h={1\over200},\qquad \kappa={1\over100000}.
\]

The log-renormalized route has

\[
 B_0={3\over2}-2\kappa,\qquad
 \delta_{\log}={1\over2}-3\kappa,                 \tag{1}
\]

with the finite polynomial `S` degrees carried in the norm. The pure-power
route uses Q184's absorption inequality and therefore has

\[
 \widetilde B_0={3\over2}-3\kappa,\qquad
 \delta_{\rm pure}={1\over2}-4\kappa>0.           \tag{2}
\]

The target step supplies pure gain `1/2-2 kappa`, so it still exceeds (2) by
`2 kappa`. Thus the simultaneous finite-cycle induction still closes. The
log route is the existing `W/M/S` construction. The pure route may reselect
finite diagonal scales and is not identified byte-for-byte with Q181's
localized field or force.

## 2. Exact diagonal stage inequality

Fix a Cartesian jet, a target power `N` and its finite losses `L,H`. For
either positive gain `delta`, choose the least nonnegative `J0` satisfying

\[
 {h\delta(J_0+1)\over2}-L-H\ge N+1.                \tag{3}
\]

Explicitly,

\[
 J_0=\max\left\{0,
 \left\lceil {2(N+1+L+H)\over h\delta}\right\rceil-1
 \right\}.                                        \tag{4}
\]

Q181 uses `J=J0+1`, so the replayed cross and regular rows retain

\[
 {h\delta\over2}>0.                               \tag{5}
\]

Equation (4) is finite for every fixed request in both routes. It does not
give a stage-uniform constant, license `J=J(q)`, or prove convergence of a raw
infinite series.

For the log route, (5) is Q181's published value

\[
 {h\delta_{\log}\over2}={49997\over40000000}.      \tag{6}
\]

For the corrected pure route it is instead

\[
 {h\delta_{\rm pure}\over2}={49996\over40000000}. \tag{7}
\]

This is the required correction. Q181's value is valid in the `W/M/S` route
but is too large by `1/40000000` if its displayed source estimate is read as
an all-band pure-power norm.

The compact `A0` and `A1` logarithms use disjoint accounts. Each still gains
`h kappa`, spends half on its own logarithm, and retains

\[
 {h\kappa\over2}={1\over40000000}.                 \tag{8}
\]

They cannot pay the regular or cutoff accounts and are not added to (6) or
(7).

## 3. Complete residual, endpoint and force replay

Replacing `delta_log` by the smaller positive `delta_pure` only changes the
finite stage selected by (4) and the reserve (7). The Q181 ledger therefore
retains, in both routes:

1. time, unit viscosity and complete pressure;
2. both old/new nonlinear orientations and the new square;
3. all six covariance entries and actual axial velocity;
4. `Pi1/Pi2`, compact pressure and the whole
   `partial_Z(J_z+c_rho P)` product;
5. every potential-cutoff curl and pressure-cutoff derivative; and
6. separate `A0/A1`, Gaussian, signed cutoff-time and exterior tails.

For every fixed jet and power, (4) gives a finite stage at which all six
groups beat the request. The corner transfer and punctured endpoint estimate
remain separate. In particular, residual flatness still does not include or
regularize the terminal origin.

Localization remains the exact twelve-term expansion for

\[
 u=cu_0+\nabla c\times A,\qquad p=cp_0,             \tag{9}
\]

including viscosity, complete pressure and all nonlinear interactions. The
apparent thirteenth term vanishes because
`(grad c cross A) dot grad c=0`. The same Borel extension argument produces
one compact smooth force after choosing the actual finite norms. It does not
construct future velocity or pressure.

## 4. Comparison and protected core

Q182 uses structural properties of the rebuilt field-force pair after the
diagonal has been completed: compact support, fixed pressure gauge, finite
energy, the common force, and the protected support hole. It does not use the
numerical value of `delta`.

Consequently, on either rebuilt route, the individual pressure equations
retain `-div f`, the common force cancels once in their difference, and
momentum plus finite energy removes the harmonic pressure gradient. The
eighth-power cutoff still gives

\[
 {1\over2}E_R'+{1\over2}D_R^2
 \le \|\nabla u(t)\|_\infty E_R+C_TR^{-2}.          \tag{10}
\]

No competitor derivative or boundary-flux bound is introduced.

The reselected annular corrections still begin at `X_a=4/101` and vanish on
the core path `X_in=1/101`. The shared positive-order background is small,
not zero, so the same conditional estimate survives:

\[
 |u_\theta(x_\tau,1-\tau)|
 \ge {e_0\over2}\tau^{-101/200}.                   \tag{11}
\]

The legal order remains: fix `T_tau<1`, send the ball radius to infinity,
then send `tau` to zero. A uniform comparison constant as `T` approaches one
is not used.

## Scope and next obligation

This checkpoint closes Q184's immediate downstream replay obligation in two
ways: the existing exact log-renormalized route and a corrected pure-power
route with (7). It does not recover Q180's old rate as an exact pure-power
bound, identify the two diagonals term by term, independently prove the exact
core or finite existence/inverse/gap lemmas, construct a future solution, or
establish an unconditional pulse bound, residual improvement, admissible
force, competitor exclusion or Navier--Stokes blow-up.

The next obligation is an independent source-level audit of the earliest
retained exact-core and finite existence/inverse/gap premise, while preserving
one of the corrected rate routes above.

Reproduce with:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q185_corrected_rate_diagonal_comparison_replay/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q185_corrected_rate_diagonal_comparison_replay/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q185_corrected_rate_diagonal_comparison_replay/evidence.json \
  --output q185_replay.json
```
