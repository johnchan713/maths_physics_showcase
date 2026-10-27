# Same-field protected-core lower-bound reconciliation

Parent: `47a30fc47a9b55ad9d81d190559775a98ee95a5e`.
Status: `same-complete-diagonal-positive-core-bound-reconciled`.

The parent closed the scale-compatible corner-to-force implication but left
the positive-core estimate for the same complete-potential diagonal as the
next obligation.  That boundary was too pessimistic.  The exact support hole,
the Cartesian core ansatz and the rule-defined Borel background estimates
already give the required lower bound when combined carefully.

There is one necessary correction to the proposed shortcut.  Every annular
diagonal increment vanishes exactly in the protected core, but positive-order
background coefficients generally do not.  They are instead
`O(q^(2h))` after the leading velocity normalization.  This checkpoint keeps
those two mechanisms separate and retains the complete fixed-gauge pressure
and field.  Historical evidence is unchanged.

## 1. The leading physical coefficient is positive

Fix any

\[
                 0<X_{\rm in}<X_a={4\over\Lambda}.       \tag{1}
\]

At `z=0` the implicit equation gives `q=tau`, `eta=0`.  On
`r=sqrt(2 X_in tau)`, the normalized Cartesian coordinates are
`xi=sqrt(2 X_in)`, `upsilon=0`.  The exact core ansatz is

\[
 u_2=q^{-1/2}{\upsilon v_0\over2}+q^{-A}\xi F,
 \qquad F={\varphi_*\Phi\over C}.                        \tag{2}
\]

Hence the radial term vanishes on this path and the azimuthal component is

\[
 q^A u_\theta=\sqrt{2X_{\rm in}}F(X_{\rm in},0)
               =E_0(X_{\rm in},0)=:e_0.                 \tag{3}
\]

The independently rechecked infinite core has `Phi>33/125`; `varphi_*` is a
positive exponential and `C>0`.  Therefore

\[
 e_0>0,\qquad
 e_0^2\ge2X_{\rm in}
 \left({\inf\varphi_*\over C}{33\over125}\right)^2>0.  \tag{4}
\]

This is a velocity coefficient identity.  The local shear, stress-growth and
cone estimates are not used as a substitute for (3).

## 2. The positive-order background is small, not zero

The rule-defined background cuts its order-`n` velocity coefficients as

\[
 q^{2nh}\chi(a_nq)(E_n,U_n,W_n^r)
\]

with the complete radial curl commutator retained.  The exact Borel selector
gives, for every fixed normalized velocity jet, a finite-head plus geometric-
tail estimate.  At order zero on the fixed core path,

\[
 \left|q^A u_{\theta,{\rm bg}}-E_0(X_{\rm in},0)\right|
 \le C_{\rm bg}q^{2h}                                   \tag{5}
\]

for some finite `C_bg` and all sufficiently small `q`.  To see the quantifier
directly, hold a finite `J>=1`: the first `J` terms are a finite
`O(q^(2h))` sum, while the selected tail is bounded by
`2^(-J)q^(h(J+1))`, also `O(q^(2h))`.  No constant uniform in coefficient
order is assumed.

Choose

\[
 q_{\rm core}=\min\left\{q_0,
       \left({e_0\over2C_{\rm bg}}\right)^{1/(2h)}\right\}.
                                                                  \tag{6}
\]

Then the normalized background swirl is at least `e0/2` in magnitude for
`0<q<=q_core`.  Equation (5), rather than exact background vanishing, is the
missing estimate in the parent checkpoint.

## 3. Every annular diagonal correction vanishes exactly

Each fixed-stage correction potential, direct swirl and pressure has a
smooth closed support in

\[
                         X_a\le X\le X_b.                \tag{7}
\]

The complete diagonal cuts the potential before taking curl:

\[
 V_{j,a}=\operatorname{curl}(\chi(aq)A_j)
              +\chi(aq)b_je_\theta.                     \tag{8}
\]

Multiplication, differentiation and the cutoff commutator do not enlarge
the closed `X` support.  Thus at the fixed `X_in<X_a`, in a whole
neighborhood,

\[
 V_{j,a_j}=0\quad\hbox{for every }j\ge1.                \tag{9}
\]

The same applies to the annular pieces already included in the synchronized
stage-zero initialization.  Local finiteness is more than enough here: every
summand is individually zero at the protected path.

Combining (3), (5) and (9) for the same complete-potential diagonal gives

\[
 q^A u^*_\theta(x_\tau,1-\tau)=e_0+R_{\rm bg}(\tau),
 \qquad |R_{\rm bg}(\tau)|\le C_{\rm bg}\tau^{2h}.      \tag{10}
\]

Therefore

\[
 \boxed{|u^*_\theta(x_\tau,1-\tau)|
              \ge{e_0\over2}\tau^{-A}}
 \quad(0<\tau\le q_{\rm core}).                         \tag{11}
\]

This discharges equations (10)--(12) of the parent checkpoint without
claiming that all background coefficients vanish.

## 4. Localization and comparison order

For the same localization radius `r_L` and time width `tau_0`, the path is on
the unit plateau whenever

\[
 \tau\le\min\left\{{\tau_0\over2},
                    {r_L^2\over8X_{\rm in}}\right\}.    \tag{12}
\]

Thus (11) survives localization after also restricting by (12).  The
expanding-ball comparison must retain its order of limits: first fix `T<1`,
then send the comparison radius to infinity, and only afterward let
`tau` tend to zero along (11).  Relative to the already joined field, force,
energy and pressure-recovery contracts, this supplies the positive-core input
to the conditional same-force competitor exclusion.

## 5. Scope and next obligation

This checkpoint proves a same-field analytical lower bound relative to the
frozen exact-core, Borel-background and support contracts.  Its exact controls
check rational thresholds and identities; they do not numerically construct
the enormous profile or establish continuum estimates by sampling.

It does not turn the collection of conditional provider contracts into an
unconditional pulse hierarchy, admissible-force construction, or blow-up
theorem.  The next justified checkpoint is a consolidated dependency and
quantifier audit of the complete field/force/comparison theorem: verify that
one selected profile, Borel sequence, pulse sequence, potential/pressure
gauge and localization is used throughout, identify every remaining
conditional continuum premise, and replay the fixed-time/radius/core limit
only after that list is empty.

## 6. Reproduction

```sh
python3 research/navier_stokes_cascade/results/core_same_field_positive_core_reconciliation/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_same_field_positive_core_reconciliation/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_same_field_positive_core_reconciliation/evidence.json \
  --output fresh.json
```
