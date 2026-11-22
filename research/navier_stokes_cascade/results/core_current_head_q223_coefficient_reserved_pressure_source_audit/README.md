# Q223: finite coefficient, reserved support and pressure-floor source audit

Parent: Q222, `f76255c013a02198843db0a94b4602b0307484b1`.
The live branch, Q222 evidence and CI, the frozen coefficient/support records,
the compact-repair hierarchy, the rule-selected Borel schedule, the uniform
primary theorem and the reserved-coordinate record were read first.  The
finalized profile, gauge, roots, supports, selected scales, frequency and
unchanged width `C^-20` are preserved.

This checkpoint source-audits the finite positive-order coefficient solve and
joins its velocity, pressure and stress supports to the reserved mean interval.
It also corrects one overbroad phrase in Q222's proposed obligation: the sharp
positive radial pressure floor is forced by the order-one term `Pi1`.  `Pi2`
must remain in the complete pressure because its constant contains later-source
memory, but the current records do not prove that the actual `K2` is nonzero
and a radially constant `Pi2` alone has zero two-point oscillation.

All conclusions below are conditional on the named continuum providers.
Exact rational tests and frozen evidence hashes verify the algebra and source
composition; they are not a numerical proof of the continuum estimates.

## 1. Fixed positive-order coefficient solve

After every lower order is finalized, the ordinary five discrepancies are
normalized in the order `(M,I,J,S,Cp)`.  At order `n>=1`, with
`lambda=2nh>0`, the correction is

\[
 z_n=\left(M,{M-J\over\lambda},I,-S,C_p\right),\qquad
 c_n=-B^{-1}z_n.                                             \tag{1}
\]

The frozen positive-mass determinant proof gives the continuum bound
`||B^-1||<800`.  Equation (1) is a finite-dimensional linear solve, so it has
no small-target condition.  It retains the raw `lambda^-1` loss and the
reciprocal plateau-amplitude jets.  Hence it supplies finite seminorms for
each selected coefficient and derivative order, conditional on the finite
inner Volterra providers; it supplies no constant uniform in coefficient
order or derivative order.

This is distinct from the leading compact repair.  Q188's leading repair uses
the nonlinear modulation and exact repair machinery with its own cone and
smallness hypotheses.  The two statements cannot be exchanged merely because
both repair five rows.

## 2. Exact support ledger and reserved separation

The correction bumps end at logarithmic radius `51/20`, the positive-order
velocity edge is `3`, and the reserved I4 interval is `(6,11)`.  Thus

\[
 {51\over20}<3<6<11.                                      \tag{2}
\]

The first exact total makes `M_n` and all available angular jets vanish after
`Xvel`.  The flux identity then makes the whole positive-order velocity zero
there.  The potential-cutoff commutator is proportional to `M_n`, so it also
vanishes on I4, including during cutoff transitions.  The exact right support
envelopes remain

| Order | Velocity | Pressure | Stress |
|---|---|---|---|
| `1` | `Xvel` | `Xv` | `Xb` |
| `2` | `Xvel` | `Xv` | `Xv` |
| `n>=3` | `Xvel` | `Xvel` | `Xvel` |

Pressure is reconstructed from its complete source and is not a velocity
component.  Therefore zero positive-order velocity on I4 does not imply zero
positive-order pressure there.  In the unchanged gauge,

\[
 \Pi_1=K_1(\eta)-{\mu^2\over4X},\qquad
 \Pi_2=K_2(\eta),\qquad \Pi_n=0\quad(n\ge3).             \tag{3}
\]

The functions `K1` and `K2` contain the later pressure-source memory.  A local
zero of the radial pressure derivative fixes neither integration function.
The historical manufactured-profile control with nonzero order-two plateau is
a counterexample to deleting such memory; it is not a computation of the
actual finalized profile's `K2`.

## 3. Sharp order-one pressure obstruction

At fixed `(q,eta)`, all `X`-independent pieces of the reserved pressure,
including both cutoff-weighted `K1` and `K2`, may be combined into an arbitrary
constant `C`.  Write

\[
 f(X)=C-s{\mu^2\over4X},\qquad
 s=|\chi(c_1q)|q^{2h}\ge0.                              \tag{4}
\]

For `0<X_a<X_b`, the exact minimax identity is

\[
 \inf_C\max\{|f(X_a)|,|f(X_b)|\}
 ={s\mu^2\over8}\left(X_a^{-1}-X_b^{-1}\right).          \tag{5}
\]

Indeed, the triangle inequality gives half the endpoint difference as a
lower bound.  Choosing `C` to center the two endpoint values makes them equal
and opposite, attaining the bound.  Taking one-sided limits to the open I4
endpoints yields

\[
 \sup_{X\in I4}|q^{2A}p_B-\Pi_0|
 \ge { |\chi(c_1q)|q^{2h}\mu^2\over8\rho_3}
       (e^{-6}-e^{-11}).                                \tag{6}
\]

Equation (6) is an analytical two-point obstruction to deleting `Pi1` or
replacing the reserved full pressure by `Pi0`.  The `1/8` factor is sharp
even after arbitrary `X`-independent cancellation.

The correction is precise: `Pi2=K2(eta)` has the same value at both radial
points, so its two-point difference and its sharp radial minimax floor are
zero.  This does not authorize deleting `Pi2`; complete covariance and mean
equations use its value and physical derivatives.  It only means that the
available source records prove no positive lower bound for the actual `K2`.
Equation (6) is also not a lower bound for the complete Navier--Stokes
residual, where other terms and cancellations must still be estimated.

## 4. Reconciled frontier

Q204 already corrected the older all-order coefficient-scale concern: the
implemented Borel construction chooses each scale after the finite
coefficient seminorms are known.  Q149 and Q204 already provide the uniform
primary-pulse value, derivative, pressure, potential and closed-support
bounds on every used physical band, with constants allowed to depend on the
fixed derivative and profile.  Repeating the fixed coefficient estimate
cannot strengthen these results into a hierarchy-uniform raw generated-pulse
constant.

The next source obligation is therefore to audit that uniform primary theorem
at its exact interfaces, then carry (3), full axial velocity and all pressure
derivatives into the complete covariance/mean and residual join.  A common raw
constant for the uncut generated hierarchy, the complete viscosity, cutoff
and all-interaction residual, and unconditional blow-up remain unproved.

## 5. Reproduction

Run

```sh
python3 test_audit.py -v
python3 audit.py --verify-record evidence.json --output fresh.json
```

from this directory.  The standard-library tests replay the five-row target,
all three support regimes, strict reserved separation, the sharp rational
minimax identity and the constant-`Pi2` negative control.  The evidence record
recursively verifies the frozen source closure.  These checks distinguish
symbolic diagnostics from the retained continuum assumptions.
