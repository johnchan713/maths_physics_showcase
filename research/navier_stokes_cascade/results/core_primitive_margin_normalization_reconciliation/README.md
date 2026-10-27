# Scale-compatible primitive margin reconciliation

Parent: `498e2f994e9acc9af070b05f3ffa4ea28da44598`.
Status: `scale-compatible-primitive-margin-certificate-reconciled`.

The parent reduces the closed-chart coefficient bounds to seven primitive
domains and summarizes their remaining nondegeneracy by one scalar minimum
`Delta_prim`. That summary is too strong: it mixes differently normalized
objects. Background gaps have a common positive floor on the accepted tail.
Global physical support distances, absolute root radii and the compact edge
weight do not. Their proofs use normalized separation, dimensionless
contraction slack and weighted seminorms.

This checkpoint corrects that premise and joins the corrected certificate to
the five primitive contracts already established on this branch. Their fixed
profiles, complete pressure memory, width, gauges and accepted tail are kept.

## 1. Uniform background part

Let the same-band contract supply `F0>=f`, `-g0_theta>=t`, `D0>=d` and
`|g0|<=Gbar`. Inside the certified perturbation radius,

\[
 F\ge f/2,\quad |g|\ge -g_\theta\ge t/2,\quad
 -N_\theta={-g_\theta\over |g|}\ge {t\over2(\bar G+1)}, \tag{1}
\]
\[
 \lambda^2\ge {ftd\over4(\bar G+1)^2}>0,
 \qquad 1-r_{\rm cone}\ge\delta_{\rm cone}>0.          \tag{2}
\]

The leading/background Borel join supplies every fixed closed-chart
seminorm before the tail is selected. Thus these are genuine uniform floors,
including at the angular poles and with complete `Pi1/Pi2` pressure memory.

## 2. Support separation is normalized

The constructive geometry uses 1,250 colors, prime `P=2003` and radius at
most `1/6057072`. Its exact center-congruence remainder is

\[
 {1\over P}-{1\over4P}={3\over8012}>0.                 \tag{3}
\]

This proves closed disjointness in normalized torus coordinates. It is not a
uniform physical Euclidean separation: after physical scale `s`, the distance
is proportional to `3s/8012`, whose infimum is zero. Hence the parent cannot
require one positive `Delta_support` interpreted in physical coordinates.
The finite graph needs closed disjointness, fixed normalized radial hulls and
flat zero extensions, all supplied by the existing support contracts.

## 3. Root inputs use relative balls

For the annular root, `delta=(143/500)epsilon`, `beta=610`, and the quadratic
bound is 50. Therefore

\[
 r_{ann}=2(610){143\over500}\epsilon,
 \qquad
 \vartheta_{ann}=8(610)^2(50){143\over500}\epsilon<1.  \tag{4}
\]

The accepted construction has `epsilon=10^-16`. The absolute ball shrinks
with its input scale; the dimensionless contraction slack is what matters.
The pointwise inverse is at most 1220 and fixed root jets follow from the same
root identity. Compact repair similarly has

\[
 \vartheta_{rep}\le2^{-219}<1,
 \qquad 2000-1600=400>0                                \tag{5}
\]

in its normalized uniqueness ball. Its absolute radius is proportional to
`D/N`; `lambda^-1<=A` is retained rather than taking `lambda=0`.

## 4. Compact inputs are weighted

The compact two-edge map uses

\[
 \zeta(x)=\exp\{-w^2/x^2-4/(L-x)^2\}.                  \tag{6}
\]

Its closed-shell infimum is zero. A positive scalar
`Delta_compact_input` would reject the intended flat-edge behavior. The
proved domain is instead

\[
 |D^m f|\le C_m\epsilon^b\zeta\delta^{-p},              \tag{7}
\]

propagated by the shifted kernels, zero-Haar temporal inverse and fixed slow
response on their normalized support hulls. Each `C_m` is finite for fixed
`m`; no pointwise lower bound for `zeta` is used.

## 5. Correct certificate

The heterogeneous scalar minimum must be replaced by

\[
 \Delta_{bg}>0,\quad {3\over8012}>0,\quad
 \vartheta_{ann}<1,\quad\vartheta_{rep}<1,\quad
 C_m^{profile},C_m^{root},C_m^{compact}<\infty
 \quad(m\text{ fixed}),                                \tag{8}
\]

plus the fixed normalized support inclusions and flat extensions. The five
existing contracts supply (8): the leading/Borel field, same-band gaps,
full-path ODE, fixed profile/projectors, and compact/slow weighted domains.
Relative to (8), the parent's finite-composition argument supplies the fixed
stage residual and increment-representative closed-chart constants.

This correction does not provide constants uniform in stage, derivative
order or scale. Nor does it revalidate downstream localization, compact
force, energy, comparison or core growth after the premise correction. No
unconditional pulse hierarchy, complete residual theorem or blow-up result
is claimed.

Next, replay the logarithmic-slack and Cartesian-corner joins using (8), then
audit localization, force, energy and comparison on the same fixed-gauge
field without replacing weighted edge control by a positive pointwise floor.

Reproduce from the repository root:

```sh
python3 research/navier_stokes_cascade/results/core_primitive_margin_normalization_reconciliation/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_primitive_margin_normalization_reconciliation/audit.py \
  --verify-record research/navier_stokes_cascade/results/core_primitive_margin_normalization_reconciliation/evidence.json \
  --output core_primitive_margin_normalization_reconciliation_replay.json
```
