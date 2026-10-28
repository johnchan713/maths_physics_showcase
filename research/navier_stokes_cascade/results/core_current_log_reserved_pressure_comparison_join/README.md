# Current log-reserved pressure recovery and comparison join

Parent: `14140ab7973092807c61503d6a9b41f50190173d`.

The parent establishes, relative to the frozen current continuum contracts,
one divergence-free compactly supported preterminal field `u`, its complete
fixed-gauge pressure `p`, a smooth compact force `f`, zero initial data, and
the correct preterminal energy bounds. This checkpoint replays pressure
recovery and expanding-ball comparison for those exact representatives.

The conclusion remains conditional on the current field, force, and
protected-core contracts in the frozen chain. The executable controls
reproduce exact exponents, rational constants, and evidence hashes. They are
not an independent continuum proof of those contracts.

## 1. Difference equation and finite energy

Fix `T<1`. Let `v,P` be smooth on compact subsets of
`R^3 x [0,T]`, solve Navier--Stokes with the same force and zero datum, and
satisfy

\[
 \sup_{0\le t\le T}\|v(t)\|_2<\infty.                   \tag{1}
\]

No global derivative, pressure-decay, or dissipation hypothesis is added.
Set

\[
 w=v-u,\qquad \pi=P-p,
 \qquad g=v\otimes v-u\otimes u
 =w\otimes w+w\otimes u+u\otimes w.                    \tag{2}
\]

The common force cancels exactly once:

\[
 w_t+\operatorname{div}g=\Delta w-\nabla\pi,
 \qquad \operatorname{div}w=0,
 \qquad -\Delta\pi=\partial_i\partial_jg_{ij}.          \tag{3}
\]

The individual pressure equations still contain the force:

\[
 -\Delta p=\partial_i\partial_j(u_i u_j)-\nabla\cdot f,
 \qquad
 -\Delta P=\partial_i\partial_j(v_i v_j)-\nabla\cdot f. \tag{4}
\]

Thus neither pressure is replaced by an unforced formula, and the force
divergence is not subtracted a second time.

## 2. Recovering the pressure gradient

The fixed-interval energy bounds make `w` uniformly `L2` and `g` uniformly
`L1`. Define

\[
 \pi_* = \sum_{i,j}R_iR_jg_{ij}.                         \tag{5}
\]

An `L1` function belongs to `H^{-s}` through its bounded Fourier transform
only for `s>3/2`; the endpoint is unavailable. Hence `pi_*` is uniformly
`H^{-2}`. For a compactly supported smooth time test `a`, momentum gives

\[
 \int a\nabla\pi\,dt
 =\Delta\int aw\,dt+\int a'w\,dt
  -\operatorname{div}\int ag\,dt\in H^{-3}.             \tag{6}
\]

Subtracting `grad pi_*` leaves a harmonic `H^{-3}` distribution. Its
Fourier transform is a weighted `L2` function supported at the single point
zero, so it vanishes. Therefore

\[
 \nabla\pi=\nabla\pi_*                                  \tag{7}
\]

as a space-time distribution. This uses momentum rather than Poisson alone
and assumes no decay or temperedness of `P`.

Finite energy is essential. Without it, `v=t e_1`, `P=-x_1` and `u=p=0`
have the same zero datum and zero force but differ by a nonzero harmonic
pressure gradient.

## 3. Compact mixed pressure and quadratic commutator

Write

\[
 g=w\otimes w+h,
 \qquad h=w\otimes u+u\otimes w.                        \tag{8}
\]

Since `u` is supported in a fixed `K subset B_L`, so is `h`. On
`R<=|x|<=2R`, for `R>=max(1,2L)`, kernel separation gives

\[
 |\pi_h(x)|\le C H_T R^{-3}.                            \tag{9}
\]

Pairing with `R^{-1}w` on the annulus gives an `O(R^{-5/2})` flux. No
isotropy of the compact mixed moment is assumed.

Use `chi_R=phi(x/R)^8`, where `phi=1` on the unit ball and vanishes outside
the ball of radius two. For the quadratic pressure, the commutator kernel is
bounded by

\[
 C|z|^{-3}\min\{|z|/R,1\}.                              \tag{10}
\]

At exponent `4/3`, the near and far radial powers are `-2/3` and `-2`, and
their coefficient is exactly four. The commutator norm is `O(R^{-3/4})`.
The eighth cutoff power is sufficient for this pairing; the seventh is not.

With

\[
 E_R=\int\chi_R|w|^2,
 \quad D_R^2=\int\chi_R|\nabla w|^2,
 \quad B_R=\|\phi_R^4w\|_6,
\]

local Sobolev gives `B_R<=C(D_R+R^-1 W_T)`. The pressure and transport
fluxes reduce to

\[
 CR^{-1}B_R^{3/2},
 \qquad CR^{-7/4}B_R^{3/4}.                              \tag{11}
\]

Young absorption leaves `R^-4` and `R^(-14/5)`; the Sobolev offsets and
compact mixed flux decay like `R^-5/2`. The coarse cutoff-Laplacian term is
`O(R^-2)` and is the limiting remainder.

## 4. Expanding-ball limit and current core

The localized difference identity yields

\[
 {1\over2}E_R'(t)+{1\over2}D_R(t)^2
 \le \|\nabla u(t)\|_\infty E_R(t)+C_TR^{-2},
 \qquad E_R(0)=0.                                       \tag{12}
\]

For `L_T=sup_[0,T] ||grad u||_infinity`, Gronwall gives

\[
 E_R(t)\le2C_TT e^{2TL_T}R^{-2}.                        \tag{13}
\]

Fix `T` first and let `R` tend to infinity. Every fixed spatial ball
eventually lies in the plateau `chi_R=1`, so (13) gives `w=0` there. The
cutoffs need not be monotone, and constants need not remain uniform as `T`
approaches one. Since `T<1` is arbitrary,

\[
 v=u\quad\hbox{on }\mathbb R^3\times[0,1).              \tag{14}
\]

Only after this radius limit may one approach the core path. The frozen
same-diagonal core bound, finalized leading-profile protected collar,
current source chain, and current localization plateau together give
conditionally

\[
 |u_\theta(x_\tau,1-\tau)|
 \ge {e_0\over2}\tau^{-A}.                              \tag{15}
\]

A globally smooth competitor is bounded on a compact neighborhood of the
terminal point, contradicting (14)--(15). Such a smooth finite-energy
same-force competitor is therefore excluded **relative to every frozen
current contract**.

## 5. Scope and next obligation

This checkpoint establishes conditional fixed-gauge pressure-gradient
recovery, exact-current preterminal comparison, and competitor exclusion.
It does not remove finite energy, construct future velocity or pressure,
independently reprove every frozen continuum contract, or establish an
unconditional admissible-force or blow-up theorem.

The next obligation is a new current-head dependency and quantifier audit
that enumerates every still-conditional primitive or operator contract. Each
such contract must be discharged analytically before any downstream
admissible-force or exclusion statement can be promoted.

Reproduce with

```sh
python3 test_audit.py -v
python3 audit.py --verify-record evidence.json --output fresh.json
```
