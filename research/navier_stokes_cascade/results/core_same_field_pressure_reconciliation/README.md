# Same-field core pressure through the completed leading backward pressure

Parent: `ec9eec28c874eb626b6e5857def432a351d4e1bf`.
Status: `same-field-core-pressure-identified-with-completed-leading-backward-pressure`.

The parent joined the certified complete-potential field to the compact force
and expanding-ball comparison, but retained identification of the exact core
pressure family with the global backward pressure as its next obligation.
That wording needs an order distinction.  The exact identity is true for the
leading `q^0` pressure on the same field and the whole protected core.  It is
false for the full corrected pressure, which retains positive-order
background, wave, mean and localization terms.  This checkpoint proves the
former and records the latter as a required part of the residual.

No gauge, profile, root, support, width, cutoff or pressure term is changed.
Historical evidence stays frozen.

## 1. The axis datum is the same exact pressure family

For the unedited reference, `y=log(X/X_R)` gives `dX/X=dy`.  The axis datum
used by the coefficient-space core is therefore

\[
 \Pi_{\rm axis}(\eta)
 =-\frac12\int_{-\infty}^{\infty}E_{\rm id,sched}(y,\eta)^2\,dy
 =-P^2\int_0^1(1+\eta^2)^{-2\theta}\,d\mu(\theta).       \tag{1}
\]

Here `mu` is the exact positive pushforward measure, including the inherited
`(5/2) delta_1` atom.  It is independent of the later attachment scale; no
diagnostic quadrature replaces it.  The fixed-point and growth rechecks use
this same measure and gauge.

The five-row repair ledger retains the pressure row as a functional identity
in `eta` and gives for the completed leading swirl `E_0`

\[
 \Pi_{\rm axis}+\int_0^\infty {E_0(x,\eta)^2\over2x}\,dx=0. \tag{2}
\]

Thus the global pressure has the same axis value as the core.  Omitting the
fifth row would leave the nonzero functional defect on the left of (2); four
moment equations or a zero pressure value at one angle cannot repair it.

## 2. Equality holds as a radial function, not only at the axis

The independently rederived core fixed point uses `Y=Lambda X`,

\[
 E_{\rm core}=\sqrt{2X}\,{\varphi_*\Phi\over C},\qquad
 \Pi_{\rm core}=\Pi_{\rm axis}+\Lambda^{-1}p,
 \qquad p_Y=(\varphi_*\Phi/C)^2.                          \tag{3}
\]

Consequently the scale factors cancel exactly:

\[
 \partial_X\Pi_{\rm core}
 ={(\varphi_*\Phi)^2\over C^2}
 ={E_{\rm core}^2\over2X}.                               \tag{4}
\]

The activation, continuation, modulation and exact repair construction does
not touch the protected inner core, and it does not reset the axis pressure.
Hence on that whole common collar `E_0=E_core`.  Equations (1) and (4) show
that the core and completed leading pressures have the same value at `X=0`
and the same radial derivative.  Therefore

\[
 \Pi_{\rm core}(X,\eta)=\Pi_0(X,\eta)
 =\Pi_{\rm axis}(\eta)+\int_0^X{E_0(x,\eta)^2\over2x}\,dx. \tag{5}
\]

Using (2) in (5) yields the same-gauge backward formula

\[
 \boxed{\Pi_{\rm core}(X,\eta)=\Pi_0(X,\eta)
 =-\int_X^\infty{E_0(x,\eta)^2\over2x}\,dx}             \tag{6}
\]

throughout the protected core.  This closes the parent's pressure
identification on the same field, rather than merely matching an axis
constant from a different profile.

Every fixed `eta` derivative in (6) is legitimate.  The density is smooth at
the regular axis, smooth on finite repaired intervals, and bounded in the
unchanged heat tail by `C_k X^(-2-2h)` for each fixed `k`.  Dominated
differentiation applies.  This is fixed-order smoothness, not a uniform
all-order constant or a positive Taylor radius at the heat endpoint.

## 3. Correction: the full corrected pressure is not the core pressure

The exact complete-potential pressure has the schematic decomposition

\[
 p_{\rm full}=\Pi_0+p_{B,+}+p_{\rm waves}+p_{\rm means}+p_{\rm loc}. \tag{7}
\]

All terms in (7) are retained by the finite-stage, diagonal, localization and
comparison checkpoints.  They cannot be erased by (6).  For example, on the
reserved I4 patch the first background pressure contains

\[
 a(q,\eta)\left(K_1(\eta)-{\mu^2\over4X}\right),
 \qquad \partial_X={a(q,\eta)\mu^2\over4X^2}>0           \tag{8}
\]

whenever its positive cutoff amplitude `a` is active.  The positive-order
background velocity vanishes there, so (8) is a direct counter-control to
collapsing the complete pressure to the leading centrifugal pressure.

The correct conclusion is therefore: the local exact core pressure equals
the completed global **leading** backward pressure on the same field; the
full pressure continues to carry every later term.  This distinction also
keeps `-div f` in the localized forced pressure equation and does not alter
the double-Riesz recovery of a competitor pressure difference.

## 4. Current analytic obstruction

The five primitive accepted-tail contracts named by the older foundation
frontier were subsequently joined on this branch: leading/background Borel
identification, same-band gaps, the full-path ODE, fixed-profile projectors,
and compact/slow inverse domains.  The finite-stage recurrence and locally
finite complete-potential diagonal then retained full pressure memory and
all residual groups.

The remaining implication used by localization is still conditional.  For
each fixed spacetime jet one needs a bound on the closed endpoint chart,
including `|eta|=1` and `T=0`, with the extra time derivative generated by
the force.  Moreover that fixed-order corner flatness must hold uniformly on
each compact `X`-region.  Pointwise or pathwise decay does not imply this
uniform statement.  In schematic form the missing estimate is

\[
 \forall m,N,K\Subset[0,\infty)_X\quad
 \sup_{X\in K,\ |\eta|\le1,\ 0<T\le T_0}
 T^{-N}\,|\partial_{X,\eta,T}^{\alpha}\mathcal R(X,\eta,T)|<\infty,
 \qquad |\alpha|\le m,                                   \tag{9}
\]

with a smooth extension at `T=0` and the corresponding `T` derivative
required by localization.  The
existing locally finite diagonal supplies finite active stages away from the
corner, but its cited fixed-stage estimates have not yet been assembled into
(9) on the closed chart.  Passing finite tests cannot provide that continuum
supremum.

Until (9) is proved, compact smooth forcing, global comparison and competitor
exclusion remain conditional.  No unconditional pulse bound, complete
residual improvement, admissible-force theorem or blow-up result is claimed.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_same_field_pressure_reconciliation/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_same_field_pressure_reconciliation/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_same_field_pressure_reconciliation/evidence.json \
  --output core_same_field_pressure_reconciliation_replay.json
```
