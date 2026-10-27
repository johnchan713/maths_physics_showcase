# Direct compact/slow provider reconstruction

Parent: `05e76e48ca83557535c30f2c0ee8994231a6a264`.
Status: `direct-compact-slow-operators-proved-with-actual-source-supremum-frontier`.

The parent reduces unconditional promotion to a fixed-stage provider estimate
for the same all-order construction.  This checkpoint reconstructs the first
two operators in that estimate from their formulas: the shifted two-edge
compact maps and the fixed five-profile slow response.  The reconstruction
keeps the physical characteristic shift, both flat edges, the actual axial
velocity, the moving pressure product and the complete pressure memory.

The operator estimates do close.  The audit also finds that the inherited
wording merges them with a different statement: uniform weighted bounds for
the actual generated sources over the infinite accepted band tail.  Smoothness
on each band does not prove that uniform supremum.  That source inequality,
given precisely in section 4, is now the first unresolved provider step.

## 1. The shifted compact map is bounded at every fixed order

In the retained characteristic coordinate `U=R^d`, let

\[
 I_-(U,y)=\int_{-\infty}^0F(U+a,y+Ma v_r)\,da,
 \qquad
 I_+(U,y)=\int_0^\infty F(U+a,y+Ma v_r)\,da,             \tag{1}
\]

and

\[
 T_eF=g_e\{(1-\chi)I_- -\chi I_+\},\qquad
 A_eF=k_e\chi_U\int_{\mathbb R}F(U+a,y+Ma v_r)\,da.     \tag{2}
\]

The shift in (1) is essential and is not removed.  A derivative in a slow,
radial or coefficient variable differentiates `F`; it does not differentiate
the fixed translation `Ma v_r`.  Auxiliary derivatives commute with torus
translation.  Thus the differentiated estimate reduces to the scalar edge
estimate for each fixed source jet.

Write

\[
 \zeta(U)=\exp\{-w^2/U^2-4/(L-U)^2\},\qquad
 \delta(U)=\min(1,U,L-U).                               \tag{3}
\]

Near `U=0`, `(1-chi)I_-` integrates only points between `0` and `U`.
Every fixed pole is absorbed because

\[
 U^{-p}e^{-c/U^2}\longrightarrow0\quad(U\downarrow0)   \tag{4}
\]

for `c>0`.  After `V=1/U`, (4) is the elementary fact that a positive
quadratic exponential dominates every fixed power.  The same argument with
`L-U` controls `chi I_+` at the outer edge.  On the compact transition where
`chi_U` is supported, `zeta` and `delta` have positive floors, so both
integrals and `A_e` have finite constants.  Splitting into these three
regions proves, for every fixed derivative request and pole,

\[
 |D^\alpha F|\le C_\alpha\epsilon^b\zeta\delta^{-p_\alpha}
 \quad\Longrightarrow\quad
 |D^\alpha(T_eF,A_eF)|
 \le C'_\alpha\epsilon^b\zeta\delta^{-p_\alpha}.       \tag{5}
\]

The constants depend on the fixed profiles, width, cutoff, order and pole.
No uniformity in derivative order is asserted.  The orientation is necessary:
`I_-` may approach a nonzero total integral as `U` approaches `L`, while the
outer weight vanishes; `I_+` has the symmetric failure at `U=0`.

For the compact remainder, the exact weighted moment condition removes the
zero Haar mode. On each nonzero Fourier mode, integration by parts `p` times
in the full-line integral in (2) gives `M^-p`; the irrational-direction
bound costs `p` Fourier derivatives, and three more make the two-dimensional
tail summable. Since `M^-1` is bounded by a fixed polynomial times
`epsilon^kappa`, this proves directly

\[
 \|A_eF\|_m\le C_{m,p}\epsilon^{b+p\kappa}P_{m,p}(S),
 \qquad J=m+2p+3.                                      \tag{6}
\]

The shift remains inside this calculation. Without the exact zero weighted
Haar moment, its zero mode survives and (6) is false for arbitrary `p`.

## 2. The zero-Haar and slow maps have finite fixed-order bounds

For a nonzero torus mode the exact quadratic-conjugate estimate remains

\[
 |v_r\!\cdot k|^{-1},|v_t\!\cdot k|^{-1}<3|k|_1.       \tag{7}
\]

Four auxiliary derivatives make the two-dimensional Fourier tail summable.
Combining (5)--(7), the actual fast block retains

\[
 g_d=-c^{-1}N_t^{-1}F_z,\qquad
 \gamma=g_d-A_1g_d,\qquad
 \beta=-\epsilon\partial_ZT_1g_d.                      \tag{8}
\]

The fixed slow response is a finite linear combination of three smooth
compact profiles and their potential on `1/4<=s<=4`.  Its only exceptional
conditioning is the already fixed finite factor `lambda^-1`.  Finite chain
and product rules therefore give, schematically,

\[
 \|(u,\gamma,\psi,\beta)_{slow}\|_m
 \le C_m\lambda^{-1}
  (\|P\|_{m+1}+\|J_z\|_{m+1}+\|J_\theta\|_{m+2}).      \tag{9}
\]

Equation (9) keeps `partial_Z(J_z+c_rho P)`, including
`(partial_Z c_rho)P`.  The response profiles have their exact five moment
rows and fixed whole-potential support.  Hence neither compact modulation nor
the slow exact repair introduces a new loss beyond its explicit finite
derivative request.

## 3. What has actually been established

Equations (5)--(9) are continuum operator estimates. Their proof uses the
actual shifted formulas, flat weights, fixed cutoffs, Fourier divisors and
fixed response profiles; it does not infer boundedness from tests.  They show
that once a matching-stage source belongs to the stated weighted space, the
compact and slow outputs do too.  Complete covariance, `g_d-A_1g_d`,
`J_z+c_rho P`, old pressure and both `Pi1/Pi2` memories remain inputs.

The exact controls freeze the edge orientation, reject a zero Gaussian
coefficient, track the slow derivative requests, and retain all six residual
groups.  They are reproducibility checks for the proof and do not establish
the source hypothesis of (5).

## 4. Exact remaining source inequality

For every fixed stage `j`, physical jet `m` and target power `N`, let
`R_j(m,N)` be the finite raw request supplied by the certified backward
recurrence.  Before (5) can be applied, the actual matching-stage generated
rows must satisfy

\[
 \boxed{
 \sup_{\ell\ge\ell_*}\epsilon_\ell^{-b_j}
 \max_{|\alpha|\le R_j(m,N)}
 \sup_{\Omega_{j,\ell}}
   \zeta^{-1}\delta^{p_{j,\alpha}}
   |D^\alpha F_{j,\ell}|<\infty .}                     \tag{10}
\]

Here `F_(j,ell)` denotes the actual source tuple at band `ell`: all six
covariance entries, the axial velocity, the current pressure moments and the
whole supported row.  It must be computed from the same fixed profile, roots,
supports and gauges used by `U^[j]`.

For each individual band, smooth compact data give a finite number in (10).
That pointwise statement does not imply that their supremum over
`ell>=ell_*` is finite.  The compact/slow operator bound also cannot prove its
own input estimate. The existing source recurrence asserts (10) through
provider contracts, but this checkpoint does not find a raw-formula proof of
the band-uniform normalized bound.  This is an insufficiency of the recorded
proof, not a counterexample to the intended construction.

Until (10) is proved, the master estimate

\[
 |\mathcal R(U^{[j]})|_m
 \le C_{j,m,N}q^{hB_j-(3/2+2h+2m+\eta)}                \tag{11}
\]

cannot be promoted from the provider chain. Equation (11) must still retain
time, viscosity, complete pressure, both nonlinear cross orientations and the
new square.  No unconditional residual, compact-force, comparison or blow-up
claim is made.

The next obligation is to expand the actual matching-stage source formulas
into (10), beginning with the corrected initial real target and then one
sourced/generic cycle.  Only after those uniform normalized source bounds are
proved should (5)--(9) be inserted into the six-group residual and the force
and comparison replays.

## 5. Reproduction

```sh
python3 research/navier_stokes_cascade/results/core_direct_compact_slow_provider_reconstruction/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_direct_compact_slow_provider_reconstruction/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_direct_compact_slow_provider_reconstruction/evidence.json \
  --output fresh.json
```
