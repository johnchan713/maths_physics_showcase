# Sharp initial target rate on the joined inverse domains

Parent: `a558aa752380e7b9aa41713287e904460dbb5ca0`.
Status: `sharper-initial-target-rate-and-strict-first-real-domain-join`.

The parent closes the compact and slow inverse domain, support and pressure-memory
contracts, but its final rate ledger reuses the obsolete downgrade

\[
  \beta=2-2\kappa\longmapsto\alpha=2-3\kappa
\]

before the linear signed-target maps.  The earlier frozen checkpoint
`core_initial_target_rate_gain` already proves that this downgrade is not needed.
This checkpoint joins that sharper rate to the parent's explicit same-profile
projectors, shifted two-edge primitive, original positive root, whole-support
extension and complete pressure memory.

The resulting first-real final row has a strict rate margin.  This is a
conditional rate correction, not a construction of the outstanding same-state
high jets and not a proof of the complete residual or blow-up.

## 1. Retain the actual Haar-row order

For the initial real-primary prefix, the frozen weighted Haar estimate gives

\[
 \|\bar E_\theta\|_n+\|\bar E_z\|_n
 \le \epsilon^\beta\zeta\delta^{-2p}P_n(S),
 \qquad \beta=2-2\kappa=\alpha+\kappa.                 \tag{1}
\]

This is the actual row with the complete primary covariance, zero-Haar fast
correlations, sharp slow fields, the actual `gd-A1 gd` axial component, compact
pressure and moving compensator.  It is special to the initial real prefix; the
full pointwise and later auxiliary-dependent rows are not promoted to (1).

The explicit same-choice projector subtracts two moment bumps of order three.
Because \(\beta<3\), its two projected sources retain order

\[
 h_e\in\mathcal M_{\min(\beta,3)}=\mathcal M_\beta.     \tag{2}
\]

Here the axial source keeps the whole derivative
\(\epsilon\partial_Z(J_z+c_\rho P)\), including
\((\partial_Zc_\rho)P\).  The logarithmic pullback, fixed projector, shifted
two-edge primitive, reciprocal edge weight and native-chart return are linear
fixed maps.  The parent's weighted-domain theorem supplies their finite
coefficient and pole losses without an epsilon loss.  Therefore

\[
 \Xi_e=-T_eh_e\in\mathcal M_\beta.                     \tag{3}
\]

The original positive-root signed inverse is likewise linear in \(\Xi\).
No new root, amplitude sign assumption, pressure truncation or unshifted
surrogate is introduced.

## 2. Complete first-real covariance rates

Using (3) in the complete phase-locked real formulas, the signed principal
velocity and complete-curl remainder have orders

\[
 \beta-\tfrac12=\tfrac32-2\kappa,
 \qquad \beta-\kappa=2-3\kappa=\alpha.                 \tag{4}
\]

The three complete covariance error groups have orders

\[
 \beta+1-2\kappa=3-4\kappa,\qquad
 2\beta-1=3-4\kappa,\qquad
 2\beta-2\kappa=4-6\kappa.                            \tag{5}
\]

Thus the real-pair remainder is order \(3-4\kappa\).  Before Haar averaging,
the selected tangential response pays one radial graph loss and has order

\[
 \beta-\kappa=\alpha.                                  \tag{6}
\]

All signed squares, curls, six covariance entries and the auxiliary-dependent
pressure remain in the response.  Equation (6) improves the stale pointwise
rate \(D=\alpha-\kappa=2-4\kappa\) by exactly \(\kappa\).

## 3. Fresh-fast and final-row join

The complete fresh-fast source is the minimum of

| Source group | Order |
|---|---:|
| original fast difference | `alpha = 2-3*kappa` |
| original slow difference | `beta = 2-2*kappa` |
| complete signed covariance response | `alpha = 2-3*kappa` |
| corrected original A1 complement | `2` |

and hence has order \(\alpha\).  The resulting fast groups are

| Fast group | Order |
|---|---:|
| slow time | `3-3*kappa` |
| base radial flux | `3-4*kappa` |
| radial viscosity | `3-5*kappa` |
| old-mean transport | `4-5*kappa` |
| changed axial pressure | `3-3*kappa` |
| new-mean square | `5-7*kappa` |

Radial viscosity is still limiting, but now lies strictly above the inherited
required order \(\nu=3-6\kappa\).  Combining it with the post-signed Haar row
of order \(3-4\kappa\), the sharp slow difference of order
\(3-3\kappa\), and the retained raw `A1` bound of order
\(3-5\kappa\), gives

\[
 E_{\mathrm{final}}\in\mathcal M_{3-5\kappa},
 \qquad (3-5\kappa)-\nu=\kappa>0.                      \tag{7}
\]

The local viscosity sharpness witness remains correct for a hypothetical
source of the stale order \(D\): such a source would indeed produce only
\(D+1-2\kappa=\nu\).  It does not block (7), because the actual source reaches
\(\alpha=D+\kappa\).  No cancellation of the nonzero selected covariance
profile is needed.

## 4. Domain, support and finite requests

The compact and slow maps are exactly those of the parent.  In particular:

- the physical characteristic shift is retained in both one-sided integrals;
- the same explicit rho-squared projector is used in the target and inverse;
- the zero-Haar temporal inverse retains the actual `-A1 gd` axial velocity;
- the fixed slow response is extended over its complete reserved potential
  hull; and
- pressure reconstruction retains old pressure, both background tails, the
  moving compensator and all covariance entries.

The parent's conservative finite request graph therefore remains sufficient:
its maximum signed current-target request is 400009 and its illustrative common
pole is 2400066.  The sharper epsilon order does not evaluate or prove those
high derivatives, and these counts are not asserted to be minimal.

## 5. Scope and next obligation

This checkpoint corrects the zero-margin conclusion in the parent.  It proves
the rational rate propagation and its compatibility with the joined inverse
domains, conditional on the frozen weighted and finite-jet premises.  The exact
controls replay eight dependency records and their declared source hashes.
Passing them is an algebraic and provenance check, not a continuum supremum.

The immediate remaining obligation is to supply the actual same-state
high-jet and whole-support inequalities for the corrected beta-order target
and retained `A1` recipe.  Then the radial `A0`, global pressure and moment
supports must be rejoined.  Later complex harmonics do not inherit the special
initial Haar rate; their complete covariance and mean corrections, all
viscosity, cutoff and interaction terms, the infinite realization and blow-up
remain open.

Reproduce with

```sh
python3 research/navier_stokes_cascade/results/core_sharp_initial_target_domain_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_sharp_initial_target_domain_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_sharp_initial_target_domain_join/evidence.json \
  --output fresh.json
```
