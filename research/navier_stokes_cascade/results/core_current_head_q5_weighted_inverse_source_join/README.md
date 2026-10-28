# Current-head Q5 weighted inverses and log-renormalized actual sources

Parent: `40172e50fe02956cac8b89e199288e97bb3e8d26`.
Status: `current-head-Q5-weighted-inverse-and-log-source-domain-bound`.

Q1--Q3 bind the exact core/root object, its regional margins and its complete
frame/path domain. Q4 chooses one fixed smooth pressure profile and explicit
same-profile projectors, then recomputes every dependent target, covariance
and pressure object. Q5 joins those same objects to the shifted compact maps,
zero-Haar fast inverse, fixed slow response and actual finite-stage source
induction.

The audit also corrects Q4's displayed next inequality. A pure-power
normalized supremum is stronger than the established source class and can
fail in the presence of a finite logarithmic coefficient. The current route
uses a log-renormalized fixed-stage space and one extra compact integration;
it does not manufacture a positive margin in the older zero-margin ledger.

## 1. Shifted compact maps on the two-edge weighted domain

In the unchanged characteristic coordinate `U=R^d`, retain

\[
 I_-(U,y)=\int_{-\infty}^0F(U+a,y+Ma v_r)\,da,
 \qquad
 I_+(U,y)=\int_0^\infty F(U+a,y+Ma v_r)\,da,            \tag{1}
\]

and

\[
 T_e f=g_e\{(1-\chi)I_- -\chi I_+\},\qquad
 A_e f=k_e\chi_U\int_{\mathbb R}F(U+a,y+Ma v_r)\,da.  \tag{2}
\]

The physical translation `y+Ma v_r` remains inside every integral. At fixed
`a` it is independent of the differentiated coefficient variables, so those
derivatives do not introduce an artificial power of `M`. Replacing (1) by
unshifted primitives would change the operator.

With the unchanged width `w=C_setup^-20`, set

\[
 \zeta=\exp(-w^2/U^2-4/(L-U)^2),\qquad
 \delta=\min(1,U,L-U).                                 \tag{3}
\]

The established inner-edge, outer-edge and interior-floor estimates give,
for every fixed derivative request and pole,

\[
 |D^\alpha F|\le C_\alpha\epsilon^b
   (1+|\log\epsilon|)^P\zeta\delta^{-p_\alpha}
 \Longrightarrow
 |D^\alpha(T_eF,A_eF)|\le C'_\alpha\epsilon^b
   (1+|\log\epsilon|)^P\zeta\delta^{-p_\alpha}.       \tag{4}
\]

The edge orientations in (2) are necessary. The `A_e` cutoff derivative is
supported on a fixed interior plateau where the weight has a positive floor.
For an auxiliary-independent input with its exact `R^e` moment zero,
`A_e=0`; Q4 supplies those exact same-choice moments for `e=1,2`. Zero Haar
mean alone does not delete the whole shifted moment.

## 2. Fast inverse, axial velocity and complete pressure

For every nonzero torus mode the retained quadratic irrational directions
satisfy

\[
 |v_r\!\cdot k|^{-1},\ |v_t\!\cdot k|^{-1}<3|k|_1.    \tag{5}
\]

Four extra torus derivatives leave a summable two-dimensional Fourier tail.
The exact fast block is

\[
 u=-c^{-1}N_t^{-1}F_\theta,\qquad
 g_d=-c^{-1}N_t^{-1}F_z,\qquad \psi=T_1g_d,
\]

\[
 \beta=-\epsilon\partial_Z\psi,\qquad
 \gamma=D_1\psi=g_d-A_1g_d.                           \tag{6}
\]

The `-A1 gd` axial velocity is generally nonzero and remains in subsequent
products and the fast-time error. At output order `m`, the tangential and
axial outputs request source order `m+4`, the radial potential requests
`m+5`, the later prefast source reaches `m+8`, and the complete axial
pressure moment

\[
 \epsilon\partial_Z(J_z+c_\rho P)                     \tag{7}
\]

reaches `m+9`. Equation (7) retains `c_{rho,Z}P`. Pressure is reconstructed
after the fast and slow steps with old pressure, both `Pi1/Pi2` tails, all
six covariance entries and the current axial velocity.

## 3. Fixed slow response and whole reserved hull

The slow inverse keeps the original three fixed response profiles and their
compact potential. With `A=1/2+h`, its distinct physical scale powers are

\[
 u=a^{-1}\{s^A P F_P+s^{A-1}J_zF_z\},\qquad
 \gamma=a^{-1}s^{A-3/2}J_\theta F_g,
\]

\[
 \psi=a^{-1}s^{A-1}J_\theta\Psi,\qquad
 \beta=-\epsilon\partial_Z\psi.                       \tag{8}
\]

The five exact moment rows and
`D1 beta+epsilon partial_Z gamma=0` follow from the fixed profile identities.
The axial block retains its fixed finite `lambda^-1` conditioning; no
uniformity as `lambda` tends to zero is claimed. The compact potential can
fill gaps between response bumps, so its support is the complete convex
potential hull, which remains strictly inside Q4's reserved annulus.

Every fixed output order produces a finite request graph. For the recorded
zero-order first-real row at `kappa=10^-5`, the sufficient raw requests
include

| request | index |
|---|---:|
| `A1` raw source | 200001 |
| `A0` raw source | 399995 |
| maximum prefast row | 400005 |
| maximum covariance coefficient | 400008 |
| maximum signed target | 400009 |

These are finite analytical requests. No selected field is numerically
evaluated at those orders, and no constant uniform in order or stage is
asserted.

## 4. Correct source norm and strict reserve

Q4 repeated the pure-power source request

\[
 \sup_{\ell\ge\ell_*}\epsilon_\ell^{-b_j}
 \max_{|\alpha|\le R}\sup_{\Omega_{j,\ell}}
 \zeta^{-1}\delta^{p_{j,\alpha}}|D^\alpha F_{j,\ell}|<\infty. \tag{9}
\]

That is not the exact current domain. For example, an otherwise admissible
family with one factor `1+|log epsilon_l|` makes (9) diverge. The example
does not identify the canonical source; it proves only that finite log degree
does not imply a zero-reserve pure-power bound.

Fix the construction stage `j`, physical jet `m` and target power `N` first.
The finite backward recurrence then returns finite derivative ceilings,
poles, harmonics and a finite log degree. The actual source formulas are
finite compositions of sums, products, coefficient derivatives, Haar
projections, radial moments and the maps in Sections 1--3. Hence they satisfy

\[
 \boxed{
 \sup_{\ell\ge\ell_*}\epsilon_\ell^{-b_j}
 (1+|\log\epsilon_\ell|)^{-P_{j,R}}
 \max_{|\alpha|\le R}\sup_{\Omega_{j,\ell}}
 \zeta^{-1}\delta^{p_{j,\alpha}}
 |D^\alpha F_{j,\ell}|<\infty,}                       \tag{10}
\]

with finite `P_(j,R)`. This includes the actual initial source and each
fixed sourced/generic cycle, all six covariance entries, `g_d-A1 g_d`, and
the whole moving pressure product. It does not insert a desired-source
callback.

One additional compact integration gains `epsilon^kappa=q^(h kappa)` and
costs two additional source derivatives for each of `A0` and `A1`. With
`s=|log q|`, half the gain absorbs every fixed log degree:

\[
 q^{h\kappa}(1+|\log\epsilon|)^P
 \le \exp\!\left({P^2\over2h\kappa}\right)
       q^{h\kappa/2}.                                  \tag{11}
\]

The other half is a strict positive reserve. The enlarged derivative
requests remain finite after `(j,m,N)` is fixed.

## 5. The old equality is retained and reclassified

For the historical positive-majorant pure-power route,

\[
 D=2-4\kappa,\qquad \nu=3-6\kappa,
\]

and the radial viscosity order is exactly

\[
 D+1-2\kappa=\nu.                                     \tag{12}
\]

Its source and viscosity margins are both zero, so that route would need
`delta_source+delta_viscosity>0`. Finite operator constants cannot change
the equality. Equations (10)--(11) use a different norm and an enlarged
finite request, so they do not contradict or erase (12). The equality is a
failure of the obsolete pure-power proof, not the current global frontier.

## 6. Result and remaining frontier

The Q1--Q4 object, shifted two-edge maps, zero-Haar inverse, fixed slow
response and actual log-renormalized source induction now form one acyclic
same-object chain on the explicit later tail. This closes the fifth primitive
contract and its fixed-stage source input analytically. Exact controls check
the request arithmetic, support/order interfaces and provenance; they are
not continuum estimates of the selected large field.

The next obligation is the complete six-group residual replay on this exact
object, with separate compact-source, regular-residual and cutoff-tail log
accounts. It must retain time, unit viscosity, full pressure, both nonlinear
cross orientations, the new square, every potential/pressure cutoff
derivative, repeated covariance means and fresh axial velocity. No
unconditional residual, force, exclusion or blow-up result is claimed.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q5_weighted_inverse_source_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q5_weighted_inverse_source_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q5_weighted_inverse_source_join/evidence.json \
  --output fresh.json
```
