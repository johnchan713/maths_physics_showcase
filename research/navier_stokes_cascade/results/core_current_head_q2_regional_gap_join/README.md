# Current-head Q2 regional gaps, factored edges and same-band join

Parent: `290f6d1b20675997b9eb47799283781220e3eef1`.
Status: `current-head-Q2-regional-continuum-factored-edges-and-same-band-gap-bound`.

The parent bound Q1: the exact pressure-family core, selected root families,
same compact repair and completed leading backward pressure now belong to one
current-head object. Q2 asks whether that object has the positive regional
value, determinant, orientation and cone margins required by the primitive
promotion, including the flat annular edges.

The existing analytical providers answer Q2, but two inherited shortcuts
cannot be used. The leading stress has zero absolute magnitude infimum at its
flat edges. Also, the older `leading_cone_margin` sufficient gate for its
already chosen tilt was not certified. This checkpoint binds the valid finite
regional cover, exact edge factors and separately certified explicit
later-tail tilt to Q1, and rechecks their common perturbation radius.

## 1. Exact finite regional cover

The protected core and first activation collar retain the exact width
`C^-20`, positive core shear/growth and `K_act<C^3`. The continuation,
shape, five-moment repair, axial, intermediate and compensated exterior
providers then cover the active region by ten analytic pieces with shared
endpoints:

`Xminus`, `Xa exp(delta)`, `Xi`, `Xsep`, `XR exp(-8)`, `XR exp(-5)`,
`XR`, `e XR`, `XR exp(T+1)`, `XR exp(T+2)`, `Xrep=Xstar`.

On every closed nonflat piece their displayed inequalities give finite
constants

\[
 f\le F_0\le\bar F,\qquad -g_{0\theta}\ge t,\qquad
 |g_0|\le\bar G,\qquad
 D_0=|g_0|^2+2F_0g_{0\theta}\ge d>0,                 \tag{1}
\]

plus orientation and same-tilt cone margins `omega,gamma>0`. These are
continuum inequalities and compact-cover arguments from the cited providers.
The interval computations freeze fixed scalar choices and outward rounding;
they are not sampled proofs of the `X,eta` infima.

The first outer-ramp comparison must use

\[
 p_1(y)\ge1.1X_R e^{-0.6y}
  -{0.01\over0.98}X_R(e^y-1).                         \tag{2}
\]

At `y=1`, (2) is greater than `0.586 X_R`, hence greater than
`X_R/2`. The historical stronger expression that multiplied both terms by
`e^-0.6` is false and is not used. Frozen historical evidence remains
unchanged.

## 2. The flat-edge obstruction and its replacement

Flatness forces

\[
 \inf_{X_a<X<X_b}|T_0|=0.                              \tag{3}
\]

Thus a positive absolute stress floor is unavailable. The exact factors are

\[
 T_0=e^{-t_1^2/y_a^2}A_a,\quad |A_a|\ge c_a>0,
\]

\[
 T_0=e^{-4/y_b^2}y_b^{-3}(b_0,y_b^6b_z),\quad b_0>0. \tag{4}
\]

For

\[
 \zeta=\exp[-t_1^2/y_a^2-4/y_b^2],
\]

(4) gives the valid bound `|T0|>=c zeta`; the normalized directions extend
smoothly through both edges. The background provider proves the relative
estimate

\[
 D^I[(\widehat T-T_0)/|T_0|]=O_I(q^{2h}),              \tag{5}
\]

while retaining `Pi1`, `Pi2`, full axial velocity and all later pressure
memory. No absolute division by a nonexistent edge floor occurs.

## 3. One positive same-band radius

For a perturbation of the value and relative-stress data, define

\[
 L_D=4\bar G+2\bar F+3,\quad L_P=\bar F+\bar G+1,
\]

\[
 \bar D=\bar G^2+2\bar F\bar G,
\quad L_c={2L_D\over ft}+{2\bar D L_P\over f^2t^2},
\quad L_H=1+{2\over t},
\]

\[
 L_\phi=3u^2L_H+(1+u^2)(4L_c+3\bar c^2L_H),
\quad m_\phi=u^2\omega^2\gamma(2-\gamma)>0.          \tag{6}
\]

Then

\[
 \rho_{gap}=\min\left\{{1\over2},{f\over2},{t\over2},
 {d\over2L_D},{\omega\over2L_H},{m_\phi\over2L_\phi}\right\}>0. \tag{7}
\]

Equation (7) preserves the value, determinant, spectral, orientation and
same-tilt cone rows. The separately certified
`explicit_tilt_pulse_tail -> core_primary_tail_recheck` route supplies the
accepted later tail used by the same-band join. The uncertified historical
existing-tilt gate is excluded rather than silently promoted.

## 4. Result and remaining frontier

Q2 is now bound to the exact Q1 object: the ten-region continuum estimates,
correct outer-ramp comparison, weighted flat-edge factors, explicit later
tilt and positive same-band perturbation radius form one acyclic provider
chain. This is an analytical join with exact interface checks. It does not
independently reprove every provider, and passing tests or interval
diagnostics are not a continuum proof.

Q3 remains: audit the original full-frame primitive and moving-normal path
ODE, complete phase jets and accepted physical domain on this same tail.
Fixed projector/compensator supports, compact and slow inverse domains,
complete-wave covariance and fresh axial/full-pressure means, and the full
viscosity/cutoff/interaction residual remain downstream. No unconditional
actual-pulse bound, residual improvement, admissible force, exclusion or
blow-up theorem follows here.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_current_head_q2_regional_gap_join/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_current_head_q2_regional_gap_join/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_current_head_q2_regional_gap_join/evidence.json \
  --output fresh.json
```
