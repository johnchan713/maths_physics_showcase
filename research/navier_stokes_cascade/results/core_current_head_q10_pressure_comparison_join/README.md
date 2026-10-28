# Current-head Q10 pressure recovery and expanding-ball comparison

Parent: `45f5106e3e805007f5895cc5330af784d18d5d9c`.
Status: `current-head-Q10-pressure-and-expanding-ball-comparison-joined-conditionally`.

Q9 gives, conditional on Q8's endpoint premise and all upstream contracts, one
divergence-free compactly supported field `u`, its complete fixed-gauge
pressure `p`, a compact smooth force `f`, zero initial data, and finite
preterminal energy and total dissipation. This checkpoint recovers the
pressure difference and proves expanding-ball comparison for those exact
representatives.

The conclusion remains conditional on the frozen Q1--Q9 continuum contracts.
Executable tests reproduce exponent arithmetic, rational constants and
provenance; they are not an independent continuum proof.

## 1. Difference equation and force cancellation

Fix `T<1`. Let `v,P` be smooth on compact subsets of
`R^3 x [0,T]`, solve Navier--Stokes with the same force and zero datum, and
satisfy

\[
 \sup_{0\le t\le T}\|v(t)\|_2<\infty.             \tag{1}
\]

No global derivative, pressure-decay or dissipation hypothesis is imposed.
Set

\[
 w=v-u,\qquad \pi=P-p,
 \qquad g=v\otimes v-u\otimes u
 =w\otimes w+w\otimes u+u\otimes w.                 \tag{2}
\]

The common force cancels exactly once:

\[
 w_t+\operatorname{div}g=\Delta w-\nabla\pi,
 \qquad \operatorname{div}w=0,
 \qquad -\Delta\pi=\partial_i\partial_jg_{ij}.     \tag{3}
\]

The individual pressure equations retain Q8's force divergence:

\[
 -\Delta p=\partial_i\partial_j(u_i u_j)-\nabla\cdot f,
 \qquad
 -\Delta P=\partial_i\partial_j(v_i v_j)-\nabla\cdot f. \tag{4}
\]

Thus the force is neither dropped from the individual equations nor
subtracted twice in the difference equation.

## 2. Finite energy removes the harmonic gradient

The fixed-interval `L2` bounds make `w` uniformly `L2` and `g`
uniformly `L1`. Define

\[
 \pi_* = \sum_{i,j}R_iR_jg_{ij}.                       \tag{5}
\]

An `L1` function embeds into `H^{-s}` only for `s>3/2`, not at the
endpoint. Hence `pi_*` is uniformly `H^{-2}`. For a compactly supported
smooth time test `a`, momentum gives

\[
 \int a\nabla\pi\,dt
 =\Delta\int aw\,dt+\int a'w\,dt
  -\operatorname{div}\int ag\,dt\in H^{-3}.         \tag{6}
\]

Subtracting `grad(pi_*)` leaves a harmonic `H^{-3}` distribution. Its
Fourier transform is a weighted `L2` function supported at the single point
zero, so it vanishes. Therefore

\[
 \nabla\pi=\nabla\pi_*                              \tag{7}
\]

as a space-time distribution. No decay or temperedness of the competitor
pressure is assumed.

Finite energy is essential. Without it, `v=t e1, P=-x1` and `u=p=0`
have the same zero datum and force but a nonzero harmonic pressure gradient.
Poisson alone cannot exclude that example.

## 3. Compact mixed pressure and quadratic commutator

Split

\[
 g=w\otimes w+h,\qquad h=w\otimes u+u\otimes w.     \tag{8}
\]

The Q9 field `u` has fixed compact support, so `h` does too. Kernel
separation on `R<=|x|<=2R` gives

\[
 |\pi_h(x)|\le C H_T R^{-3},                          \tag{9}
\]

and pairing with `R^{-1}w` gives an `O(R^{-5/2})` annular flux. No
isotropy of the compact mixed moment is assumed.

Use `chi_R=phi(x/R)^8`, where `phi=1` on the unit ball and vanishes
outside the ball of radius two. The quadratic-pressure commutator kernel is

\[
 C|z|^{-3}\min\{|z|/R,1\}.                           \tag{10}
\]

At exponent `4/3`, the near and far radial powers are `-2/3` and `-2`;
the radial coefficient is exactly four, and the commutator norm is
`O(R^{-3/4})`. The eighth cutoff power is sufficient for this pairing; the
seventh is not.

With

\[
 E_R=\int\chi_R|w|^2,\qquad
 D_R^2=\int\chi_R|\nabla w|^2,\qquad
 B_R=\|\phi_R^4w\|_6,
\]

local Sobolev gives `B_R<=C(D_R+R^-1 W_T)`. The pressure and transport
fluxes reduce to

\[
 CR^{-1}B_R^{3/2},\qquad CR^{-7/4}B_R^{3/4}.           \tag{11}
\]

Young absorption leaves `R^-4` and `R^(-14/5)`. The Sobolev offsets and
compact mixed-pressure flux are `O(R^-5/2)`; the cutoff-Laplacian term is
`O(R^-2)` and is the limiting remainder.

## 4. Radius limit and conditional comparison

The localized difference identity yields

\[
 {1\over2}E_R'(t)+{1\over2}D_R(t)^2
 \le \|\nabla u(t)\|_\infty E_R(t)+C_TR^{-2},
 \qquad E_R(0)=0.                                     \tag{12}
\]

For `L_T=sup_[0,T] ||grad u||_infinity`, Gronwall gives

\[
 E_R(t)\le2C_TT e^{2TL_T}R^{-2}.                      \tag{13}
\]

Fix `T<1` first and send `R` to infinity. Every fixed spatial ball
eventually lies in the cutoff plateau, so `w=0` there. The cutoffs need not
be monotone and the constants need not be uniform as `T` approaches one.
Since `T<1` is arbitrary,

\[
 v=u\quad\hbox{on }\mathbb R^3\times[0,1).          \tag{14}
\]

Only after this radius limit may one approach the protected core.

## 5. Conditional protected-core consequence

The same-diagonal core estimate, finalized leading-profile collar, current
source chain and localization plateau give, relative to their frozen
contracts,

\[
 |u_\theta(x_\tau,1-\tau)|
 \ge {e_0\over2}\tau^{-A}.                           \tag{15}
\]

A globally smooth competitor is locally bounded near the terminal point, so
(14)--(15) exclude a globally smooth finite-energy competitor with the same
zero datum and force, conditional on every Q1--Q10 premise.

The still-open endpoint condition is unchanged:

\[
 \boxed{
 \sup_{E\times[1-\delta,1)}
 |D_x^\alpha\partial_t^j(A,B,p)|<\infty
 }                                                     \tag{16}
\]

for every fixed order on each punctured terminal compact.

## 6. Scope and next obligation

Q10 establishes conditional pressure-gradient recovery, same-force
preterminal comparison, and finite-energy competitor exclusion. It does not
remove finite energy, prove (16), construct future velocity or pressure,
independently prove every upstream continuum contract, or establish an
unconditional admissible-force, exclusion or blow-up theorem.

The next justified obligation is a current-head dependency and quantifier
audit enumerating every remaining conditional primitive and operator
contract. Each must be discharged analytically before any downstream
statement can be promoted.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q10_pressure_comparison_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q10_pressure_comparison_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q10_pressure_comparison_join/evidence.json \
  --output fresh.json
```
