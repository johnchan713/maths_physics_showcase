# Compact and slow inverse weighted-domain join

Parent: `4cab7bfa11a69e62c7ef2b5d5b892eb2d47fa967`.
Status: `compact-slow-weighted-domain-and-support-contract-joined`.

The parent closes the fixed compensator/projector contract by selecting one
smooth compact pressure profile, choosing explicit same-profile projectors,
and recomputing every projector-dependent signed, covariance and pressure
object. This checkpoint joins those objects to the original shifted compact
maps, zero-Haar temporal inverse and fixed five-profile slow response on one
common weighted domain.

This closes the fifth primitive contract from the corrected foundation
frontier. The conclusion is finite-order and conditional on the previously
established leading/background, gap, canonical ODE and current weighted-row
contracts. It does not provide a strict residual gain. The exact remaining
zero-margin inequality is recorded in section 5.

## 1. Weighted compact maps keep the physical shift

Use the retained characteristic coordinate `U=R^d` and the actual shifted
one-sided integrals

\[
 I_-(U,y)=\int_{-\infty}^0F(U+a,y+Ma v_r)\,da,
 \qquad
 I_+(U,y)=\int_0^\infty F(U+a,y+Ma v_r)\,da.             \tag{1}
\]

The compact primitive and error are

\[
 T_e f=g_e\{(1-\chi)I_- -\chi I_+\},
 \qquad A_e f=k_e\chi_U\int_{\mathbb R}F(U+a,y+Ma v_r)\,da. \tag{2}
\]

At fixed `a` the translation in (1) is independent of every differentiated
coefficient variable. Thus fixed coefficient derivatives introduce no
power of `M`. Removing the shift would change the operator identity and is
not an admissible estimate shortcut.

On the original logarithmic shell write

\[
 \zeta=\exp(-w^2/x^2-4/(L-x)^2),\qquad
 \delta=\min(1,x,L-x),                                  \tag{3}
\]

with the unchanged width `w=C_setup^-20`. For a source jet bounded by
`epsilon^b zeta delta^-p`, the established two-edge theorem bounds both
directed integrals in (2) by the same weight and pole. Its coefficient is a
finite positive sum of three terms `c_j exp(e_j)`: the inner edge, outer
edge and interior weight-floor branches. These exponentials may be enormous,
but they remain exact symbolic constants. Their numerical materialization is
not an analytical hypothesis.

The cutoff derivative in `A_e` is supported in the fixed interior plateau,
so the output weight is recovered there. The whole shifted moment is kept;
zero Haar mean alone does not delete it. For an auxiliary-independent source
with exact weighted radial moment zero, `A_e=0` and the shifted map reduces
to the native two-edge primitive. The new rho-squared projectors give those
two exact moments for `e=1,2` on the same actual moving support.

## 2. Fast inverse and full pressure domain

For the actual irrational directions, the integer quadratic conjugate gives

\[
 |v_r\!\cdot k|^{-1}, |v_t\!\cdot k|^{-1}<3|k|_1,
 \qquad k\ne0.                                          \tag{4}
\]

Four additional torus derivatives leave a summable two-dimensional Fourier
tail. Hence the zero-Haar temporal inverse preserves the slow/radial weighted
support and costs source order `m+4`. For the current rows
`F_theta,F_z`, the exact fast update is

\[
 u=-c^{-1}N_t^{-1}F_\theta,\quad
 g_d=-c^{-1}N_t^{-1}F_z,\quad
 \psi=T_1g_d,
\]
\[
 \beta=-\epsilon\partial_Z\psi,\qquad
 \gamma=D_1\psi=g_d-A_1g_d.                             \tag{5}
\]

The actual axial velocity includes `-A_1 g_d`. It remains in all later
products and in the exact fast-time error. Tangential and axial outputs use
source order `m+4`; the radial potential uses `m+5`. After the fresh scalar
and slow steps, output order `m` reaches prefast sources through `m+8`, and
the whole axial pressure moment

\[
 \epsilon\partial_Z(J_z+c_\rho P)                       \tag{6}
\]

reaches `m+9`. Formula (6) retains
`(partial_Z c_rho)P`; freezing the compensator during differentiation gives
a different expression.

The initial current Haar-row checkpoint supplies the required weighted
source inequalities before unweighting, including the full covariance,
moving pressure and fixed-map outputs. The same-choice projector checkpoint
removes the former arbitrary projector callback. Therefore the symbolic
compact-map constants apply to the actual first-real norm chain, conditional
on the cited primitive bounds, rather than to an unweighted surrogate.

## 3. Fixed slow response on the complete reserved hull

The slow inverse uses the three fixed response profiles and their compact
potential from the original five-map. With `A=1/2+h`, its physical scale
factors are distinct:

\[
 u=a^{-1}\{s^A P F_P+s^{A-1}J_zF_z\},
 \quad \gamma=a^{-1}s^{A-3/2}J_\theta F_g,
\]
\[
 \psi=a^{-1}s^{A-1}J_\theta\Psi,
 \qquad \beta=-\epsilon\partial_Z\psi.                  \tag{7}
\]

The exact five moment rows and `D_1 beta+epsilon partial_Z gamma=0` follow
from the fixed profile identities. The axial block keeps its unavoidable
fixed `lambda^-1` conditioning. Since `lambda>0` is already fixed, this is a
finite constant, without uniformity as `lambda` tends to zero.

The response is supported on the entire convex potential hull, including
gaps between the individual bumps. That hull lies strictly inside the same
reserved annulus. The reserved weight comparison is used only on this actual
support, never on the rectangular union hull. Finite profile seminorms and
the scale-coordinate composition give all fixed mixed jets. The compact
maps can fill the auxiliary torus and radial gaps, but their smooth zero
extensions remain inside the declared whole potential hull.

## 4. Finite first-real request graph

For every fixed output order the new dependency graph is finite. At
`kappa=10^-5` and output zero, the least strict raw compact counts are

\[
 p_1=99999,\qquad p_0=199996,                            \tag{8}
\]

which give

| request | index |
|---|---:|
| `A1` raw source | 200001 |
| `A0` raw source | 399995 |
| maximum new prefast row | 400005 |
| maximum new covariance coefficient | 400008 |
| maximum signed current-target input | 400009 |

For the documented illustrative nested pole `p=1`, quotient allowance
`d=12` and old pole 30, one sufficient common pole is `q=2400066`.
These integers describe finite analytical requests; no factorial, giant
exponential or physical field at those orders is numerically evaluated.
Constants may depend on the fixed order, pole, width, tilt and profiles.
There is no order-uniform or infinite-stage conclusion.

The native projected primitives, weighted fast maps, slow response and
whole-support extension all accept these finite domains. Full pressure is
reconstructed after each fast and slow update, preserving old pressure,
both background pressure tails, the moving compensator, actual `A0/A1`
velocities and all six covariance entries. This establishes the requested
compact/slow inverse domain and support contract.

## 5. Exact obstruction after the primitive join

Write

\[
 B=1-\kappa,\quad \alpha=2-3\kappa,\quad
 D=2-4\kappa,\quad \nu=3-6\kappa.                       \tag{9}
\]

The complete first-real signed covariance response supplies the new fast
source only at order `D`; its positive-majorant source margin is zero. In
the full new-fast row difference, the radial viscosity has order

\[
 D+1-2\kappa=\nu,                                       \tag{10}
\]

again with zero margin. All other listed groups have positive margins.
If `delta_source` is an additional gain in the complete fast source and
`delta_viscosity` is an additional gain or exact cancellation in its radial
viscosity contribution, strict absorption requires precisely

\[
       \boxed{\delta_{\rm source}+\delta_{\rm viscosity}>0}. \tag{11}
\]

Neither term is proved positive by the current evidence. Finite operator
constants cannot convert equality in (10) into (11). This is the next
analytical obstruction, rather than a failure of compact-map existence.

## 6. Scope and verification

This checkpoint closes the fifth primitive contract on the separately
certified later tail. It does not identify that tail with an unspecified
historical first band, prove actual uniformity over derivative order or
stage, or establish a strict complete residual improvement. Later generated
harmonics, complete repeated covariance/mean corrections, every viscosity,
cutoff and interaction term, admissible compact forcing and blow-up remain
open.

The exact controls verify dependency order, symbolic two-edge constants,
distinct fast and slow request counts, full moving-pressure differentiation,
the finite request graph and the zero-margin equality. They are algebra and
provenance diagnostics. The continuum conclusions are the cited weighted
kernel, Fourier, finite-composition and support arguments.

Reproduce with

```sh
python3 research/navier_stokes_cascade/results/core_compact_slow_inverse_domain_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_compact_slow_inverse_domain_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_compact_slow_inverse_domain_join/evidence.json \
  --output core_compact_slow_inverse_domain_join_fresh.json
```
