# Current-head primitive frontier after the five reconstructions

Parent: `fa006e2fbbea66396bf8843cfa36a81cb5c956be`.
Status: `current-head-five-contract-frontier-corrected-to-lower-core-root-regional-frame-antecedents`.

The older current-head audit listed five live continuum contracts: the
leading/background Borel construction, same-band gaps, full-path ODE,
fixed-profile projectors, and compact/slow weighted inverses. That list is no
longer current. Six later checkpoints reconstruct the first contract's two
parts and the remaining four contracts on the same objects; the current
log-source, residual, corner, localization/force, energy and comparison
checkpoints then replay the downstream chain.

This audit verifies that statement against the full current evidence graph and
moves the frontier to the antecedents that those reconstructions still state
explicitly. It does not turn dependency hashes, controls or numerical
fixtures into continuum proofs.

## 1. What the current closure actually contains

Starting at the current pressure/comparison evidence, the transitive graph
contains all of

1. the direct fixed-order background seminorm and real Borel selector;
2. the finalized leading-profile certificate;
3. the same-band continuum-gap reconstruction;
4. the canonical full-path ODE reconstruction;
5. the fixed-profile compensator/projector reconstruction; and
6. the compact/slow weighted-inverse reconstruction.

It also contains the six current downstream replays from actual log-weighted
sources through comparison. Thus it would be incorrect to repeat the old
five-contract list as the live frontier. Historical records remain frozen and
retain their narrower scope.

## 2. The lower frontier

The reconstructed contracts are conditional on three grouped antecedents.

**Q1: exact core and exact roots.** The pressure-family Banach fixed point,
the identification of its pressure with the completed global backward
pressure, and every finite matching/repair root must satisfy the exact
defining equations on one profile. Sampled roots or small floating residuals
cannot substitute for these identities.

**Q2: regional continuum margins.** On the actual finite regional cover one
needs positive constants

\[
 f=\inf F_0,\quad t=\inf(-g_{0\theta}),\quad
 d=\inf\bigl(|g_0|^2+2F_0g_{0\theta}\bigr),
 \quad \omega,\gamma,c_a,b_0>0. \tag{1}
\]

The flat-edge stress has zero absolute infimum. The valid assertion is
`|T0| >= c zeta` together with positive factored coefficients `c_a,b_0` and
smooth limiting directions. A grid cannot prove the infima in (1).

**Q3: original full-frame primitive domain.** The original background/frame
primitive estimates, complete path and phase jets, positive carrier/root
floors, and the accepted physical tail must hold for the same field. In
particular the actual constants used in the normal estimate and tangent
propagator must be finite on the whole path, including full pressure and both
axial transport terms.

These are the live promotion premises. This checkpoint identifies and binds
their witness records; it does not independently rederive Q1--Q3.

## 3. Precise blocking inequality

After Q1 fixes the exact object, Q2 and Q3 must supply

\[
 m_{base}=\min\{f,t,d,\omega,m_\phi,c_a,b_0,
                 \lambda_{min},b_{min}^2\}>0, \qquad
 0<C_i<\infty. \tag{2}
\]

Here the finite `C_i` are the value, relative-stress, chart, length and frame
constants. Only then can the common-tail gap row select an integer `ell`:

\[
 2h(\ell-2)\ge \max_i\log_2(C_i/r_i),
 \qquad 0<r_i\le m_{base}. \tag{3}
\]

Equation (3) is the precise obstruction. If an actual `r_i` is not proved
positive or a corresponding `C_i` is not proved finite, its right-hand side
is not a finite real number. Taking `q` smaller cannot prove either missing
antecedent. Conversely, after finitely many positive margins and finite
constants are established, the ceiling in (3) is finite; no constant uniform
in coefficient stage or derivative order is required.

The executable control checks this logical domain and rejects zero, negative,
infinite and NaN margins. It does not evaluate the actual continuum values.

## 4. Quantifiers and scope

The exact profile, roots, supports and gauges are fixed first. Actual
continuum infima precede the common selector; the selector tail precedes the
band limit. For source and residual estimates, fix the physical jet, finite
stage and target before taking the small-scale limit. Comparison fixes
`T<1`, then sends the ball radius to infinity, then uses the protected core.

The current downstream implication is therefore complete only relative to
Q1--Q3. No unconditional actual-pulse bound, complete residual improvement,
admissible-force theorem, competitor exclusion or blow-up theorem follows
from this audit.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_primitive_frontier_audit/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_primitive_frontier_audit/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_primitive_frontier_audit/evidence.json \
  --output fresh.json
```
